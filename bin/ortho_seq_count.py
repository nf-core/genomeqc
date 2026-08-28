#!/usr/bin/env python3

# Written by Chris Wyatt, released under the MIT license.
# Counts, per species, how many sequences (chromosomes/scaffolds/contigs) carry a
# gene orthologous to another species, from an OrthoFinder orthogroups table and
# per-species GFF files.

import argparse
import os
import re
import sys
from collections import defaultdict

import pandas as pd


def load_gff_files(gff_files):
    found = []
    for item in gff_files:
        if os.path.exists(item) and item.endswith(".gff"):
            found.append(item)
    return found


def load_orthogroups(orthogroups_file):
    if not os.path.exists(orthogroups_file):
        print(f"[ERROR] {orthogroups_file} not found")
        sys.exit(1)

    try:
        orthogroups_df = pd.read_csv(orthogroups_file, sep='\t')
        print(f"[INFO] Loaded {len(orthogroups_df)} orthogroups")
        print(f"[INFO] Orthogroup columns: {list(orthogroups_df.columns)}")
    except Exception as e:
        print(f"[ERROR] Failed to load {orthogroups_file}: {e}")
        sys.exit(1)

    orthogroup_species = [col for col in orthogroups_df.columns if col != 'Orthogroup']
    print(f"[INFO] Species in orthogroups: {orthogroup_species}")
    return orthogroups_df, orthogroup_species


def species_name_from_gff(gff_file):
    species_name = os.path.basename(gff_file)
    for suffix in [".gff", "_busco.gff", "_busco_combined.gff", ".longest"]:
        species_name = species_name.replace(suffix, "")
    species_name = species_name.replace("_busco", "")
    return species_name


def match_species(species_name, orthogroup_species, debug_info):
    for ortho_species in orthogroup_species:
        if species_name == ortho_species:
            return ortho_species
        if species_name in ortho_species or ortho_species in species_name:
            print(f"[INFO] Using partial match: '{species_name}' -> '{ortho_species}'")
            debug_info.append(f"Partial match: '{species_name}' -> '{ortho_species}'")
            return ortho_species
    return None


def extract_gene_id(attributes, is_busco):
    if is_busco:
        match = re.search(r'Name=([^;]+)', attributes)
        return match.group(1) if match else None

    for pattern in [r'ID=([^;]+)', r'Name=([^;]+)', r'gene_id=([^;]+)']:
        match = re.search(pattern, attributes)
        if match:
            return match.group(1)

    for attr in attributes.split(';'):
        attr = attr.strip()
        if '=' in attr:
            key, value = attr.split('=', 1)
            if key.lower() in ['id', 'name', 'gene_id']:
                return value
    return None


def parse_gff(gff_file, final_species_name, debug_info):
    gene_to_chr = {}
    gene_count = 0
    is_busco = False

    with open(gff_file, 'r') as f:
        first_line = f.readline().strip()
        if "BUSCO" in first_line or "busco" in first_line:
            is_busco = True
            print(f"[INFO] Detected BUSCO file format for {gff_file}")
            debug_info.append(f"Detected BUSCO format for {gff_file}")

    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                # "transcript" is treated as a synonym for "mRNA": AUGUSTUS (and
                # AGAT-normalised multi-source annotations combining it with other
                # evidence, e.g. GeneMark) commonly emits "transcript" instead of
                # "mRNA" for the same SO term, and that's the feature gffread keys
                # protein/CDS sequence names on.
                feature_type = parts[2]
                if feature_type not in ['gene', 'mRNA', 'transcript', 'CDS']:
                    continue

                chromosome = parts[0]
                attributes = parts[8]

                gene_id = extract_gene_id(attributes, is_busco)
                if gene_id:
                    gene_id = gene_id.strip().strip('"').strip("'")
                    gene_id = gene_id.replace(":", "_")  # avoid mismatches with orthofinder
                    gene_to_chr[gene_id] = (final_species_name, chromosome)
                    gene_count += 1

                    if gene_count <= 5:
                        debug_info.append(f"  Gene {gene_count}: {gene_id} -> {final_species_name}:{chromosome}")
    except Exception as e:
        print(f"[WARNING] Error processing {gff_file}: {e}")
        debug_info.append(f"ERROR processing {gff_file}: {e}")

    print(f"[INFO] Mapped {gene_count} genes from {final_species_name}")
    debug_info.append(f"Mapped {gene_count} genes from {final_species_name}")
    return gene_to_chr


def build_gene_to_chr(gff_files, orthogroup_species, debug_info):
    gene_to_chr = {}
    for gff_file in gff_files:
        species_name = species_name_from_gff(gff_file)
        print(f"[INFO] Processing {gff_file} for species '{species_name}'")
        debug_info.append(f"Processing {gff_file} -> species '{species_name}'")

        matching_species = match_species(species_name, orthogroup_species, debug_info)
        if not matching_species:
            print(f"[WARNING] No matching species found in orthogroups for '{species_name}'")
            debug_info.append(f"WARNING: No match for '{species_name}' in orthogroups")
            continue

        debug_info.append(f"Final species name: '{matching_species}'")
        gene_to_chr.update(parse_gff(gff_file, matching_species, debug_info))

    return gene_to_chr


def summarise_overlap(orthogroups_df, gene_to_chr, debug_info):
    orthogroup_genes = set()
    for _, row in orthogroups_df.iterrows():
        for species, genes in row.items():
            if species == "Orthogroup":
                continue
            if pd.isna(genes) or str(genes).strip() == "":
                continue
            for gene in str(genes).split(","):
                orthogroup_genes.add(gene.strip())

    mapped_genes = set(gene_to_chr.keys())
    overlap = orthogroup_genes.intersection(mapped_genes)

    print(f"[INFO] Genes in orthogroups: {len(orthogroup_genes)}")
    print(f"[INFO] Genes mapped from GFF: {len(mapped_genes)}")
    print(f"[INFO] Overlap: {len(overlap)}")

    debug_info.append(f"Genes in orthogroups: {len(orthogroup_genes)}")
    debug_info.append(f"Genes mapped from GFF: {len(mapped_genes)}")
    debug_info.append(f"Overlap: {len(overlap)}")

    if gene_to_chr:
        print("[INFO] Sample gene mappings from GFF:")
        for i, (gene, (sp, chrom)) in enumerate(gene_to_chr.items()):
            if i < 10:
                print(f"  {gene} -> {sp}:{chrom}")
                debug_info.append(f"  GFF gene: {gene} -> {sp}:{chrom}")

    if orthogroup_genes:
        print("[INFO] Sample genes from orthogroups:")
        for gene in list(orthogroup_genes)[:10]:
            print(f"  {gene}")
            debug_info.append(f"  Orthogroup gene: {gene}")

    return orthogroup_genes


def compute_orthologous_pairs(orthogroups_df, gene_to_chr, debug_info):
    orthologous_chr_pairs = defaultdict(int)
    species_chr_syntenic = defaultdict(set)

    processed_orthogroups = 0
    for _, row in orthogroups_df.iterrows():
        gene_chroms = defaultdict(set)

        for species, genes in row.items():
            if species == "Orthogroup":
                continue
            if pd.isna(genes) or str(genes).strip() == "":
                continue
            for gene in str(genes).split(","):
                gene = gene.strip()
                if gene in gene_to_chr:
                    sp, chrom = gene_to_chr[gene]
                    gene_chroms[sp].add(chrom)

        species_list = list(gene_chroms.keys())
        if len(species_list) >= 2:
            processed_orthogroups += 1
            for i in range(len(species_list)):
                for j in range(i + 1, len(species_list)):
                    sp1, sp2 = species_list[i], species_list[j]
                    for chr1 in gene_chroms[sp1]:
                        for chr2 in gene_chroms[sp2]:
                            key = tuple(sorted([(sp1, chr1), (sp2, chr2)]))
                            orthologous_chr_pairs[key] += 1
                            species_chr_syntenic[sp1].add(chr1)
                            species_chr_syntenic[sp2].add(chr2)

    print(f"[INFO] Processed {processed_orthogroups} orthogroups with multi-species genes")
    debug_info.append(f"Processed {processed_orthogroups} orthogroups with multi-species genes")
    return orthologous_chr_pairs, species_chr_syntenic


def build_output_frames(orthologous_chr_pairs, species_chr_syntenic):
    if orthologous_chr_pairs:
        pairwise_output = pd.DataFrame([
            {"Species1": k[0][0], "Chr1": k[0][1], "Species2": k[1][0], "Chr2": k[1][1], "Orthogroup_Count": v}
            for k, v in orthologous_chr_pairs.items()
        ])
        species_chr_count = pd.DataFrame([
            {"Species": sp, "Syntenic_Chromosomes": len(chrs)}
            for sp, chrs in species_chr_syntenic.items()
        ])
    else:
        pairwise_output = pd.DataFrame(columns=["Species1", "Chr1", "Species2", "Chr2", "Orthogroup_Count"])
        species_chr_count = pd.DataFrame(columns=["Species", "Syntenic_Chromosomes"])

    return pairwise_output, species_chr_count


def main():
    ap = argparse.ArgumentParser(description="Count, per species, sequences carrying an orthologous gene.")
    ap.add_argument("--orthogroups", required=True, help="OrthoFinder orthogroups TSV (e.g. N0.tsv)")
    ap.add_argument("--gff", nargs='+', required=True, help="Per-species GFF files")
    ap.add_argument("--pairwise-out", default="pairwise_ortho_seq_count.tsv")
    ap.add_argument("--species-out", default="species_ortho_seq_count.tsv")
    ap.add_argument("--debug-out", default="debug_gene_mapping.txt")
    args = ap.parse_args()

    print("[INFO] Starting sequence orthology summary script...")

    gff_files = load_gff_files(args.gff)
    if not gff_files:
        print("[ERROR] No GFF files found in the input")
        sys.exit(1)

    print(f"[INFO] Found {len(gff_files)} GFF files:")
    for f in gff_files:
        print(f" - {f}")

    orthogroups_df, orthogroup_species = load_orthogroups(args.orthogroups)

    debug_info = []
    gene_to_chr = build_gene_to_chr(gff_files, orthogroup_species, debug_info)
    print(f"[INFO] Total mapped genes: {len(gene_to_chr)}")
    debug_info.append(f"Total mapped genes: {len(gene_to_chr)}")

    summarise_overlap(orthogroups_df, gene_to_chr, debug_info)

    orthologous_chr_pairs, species_chr_syntenic = compute_orthologous_pairs(orthogroups_df, gene_to_chr, debug_info)
    pairwise_output, species_chr_count = build_output_frames(orthologous_chr_pairs, species_chr_syntenic)

    pairwise_output.to_csv(args.pairwise_out, sep='\t', index=False)
    species_chr_count.to_csv(args.species_out, sep='\t', index=False)

    with open(args.debug_out, "w") as f:
        f.write("\n".join(debug_info))

    print("\n[INFO] Results saved:")
    print(f"- {args.pairwise_out}")
    print(f"- {args.species_out}")
    print(f"- {args.debug_out}")

    print("\n[INFO] Pairwise sequence orthology summary:")
    print(pairwise_output.head())

    print("\n[INFO] Species orthologous sequence counts:")
    print(species_chr_count)


if __name__ == "__main__":
    main()
