process ORTHOLOGOUS_CHROMOSOMES {
    tag "orthologous_chromosomes"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pip_pandas:2fd05a70c67560f2"

    input:
    path orthogroups_tsv
    path gff_files

    output:
    path "species_orthologous_chromosomes.tsv", emit: species_summary
    path "pairwise_chromosome_orthology.tsv", emit: pairwise_summary
    path "debug_gene_mapping.txt", emit: debug_info
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/bin/bash

    python3 << 'EOF'
import os
import sys
from collections import defaultdict
import re

try:
    import pandas as pd
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", "pandas"])
    import pandas as pd

print("[INFO] Starting chromosome orthology summary script...")

# Check that the expected files exist
orthogroups_file = "${orthogroups_tsv}"
if not os.path.exists(orthogroups_file):
    print(f"[ERROR] {orthogroups_file} not found")
    sys.exit(1)

# Get GFF files from input
gff_files = []
for item in "${gff_files}".split():
    if os.path.exists(item) and item.endswith(".gff"):
        gff_files.append(item)

if not gff_files:
    print("[ERROR] No GFF files found in the input")
    sys.exit(1)

print(f"[INFO] Found {len(gff_files)} GFF files:")
for f in gff_files:
    print(f" - {f}")

# Load the orthogroups file
try:
    orthogroups_df = pd.read_csv(orthogroups_file, sep='\\t')
    print(f"[INFO] Loaded {len(orthogroups_df)} orthogroups")
    print(f"[INFO] Orthogroup columns: {list(orthogroups_df.columns)}")

    # Get species names from orthogroups (excluding 'Orthogroup' column)
    orthogroup_species = [col for col in orthogroups_df.columns if col != 'Orthogroup']
    print(f"[INFO] Species in orthogroups: {orthogroup_species}")

except Exception as e:
    print(f"[ERROR] Failed to load {orthogroups_file}: {e}")
    sys.exit(1)

# Parse GFF files to create gene to chromosome mapping
gene_to_chr = {}
debug_info = []

for gff_file in gff_files:
    # Extract species name from filename - handle both BUSCO and regular naming
    species_name = os.path.basename(gff_file)

    # Remove common suffixes
    for suffix in [".gff", "_busco.gff", "_busco_combined.gff", ".longest"]:
        species_name = species_name.replace(suffix, "")

    # Remove _busco suffix specifically
    species_name = species_name.replace("_busco", "")

    print(f"[INFO] Processing {gff_file} for species '{species_name}'")
    debug_info.append(f"Processing {gff_file} -> species '{species_name}'")

    # Check if this species name matches any orthogroup column
    matching_species = None
    for ortho_species in orthogroup_species:
        if species_name == ortho_species:
            matching_species = ortho_species
            break
        # Also try partial matching (in case of slight differences)
        elif species_name in ortho_species or ortho_species in species_name:
            matching_species = ortho_species
            print(f"[INFO] Using partial match: '{species_name}' -> '{ortho_species}'")
            debug_info.append(f"Partial match: '{species_name}' -> '{ortho_species}'")
            break

    if not matching_species:
        print(f"[WARNING] No matching species found in orthogroups for '{species_name}'")
        debug_info.append(f"WARNING: No match for '{species_name}' in orthogroups")
        continue

    # Use the matching species name from orthogroups
    final_species_name = matching_species
    debug_info.append(f"Final species name: '{final_species_name}'")

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
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\\t')
                if len(parts) < 9:
                    continue

                # Look for gene features
                feature_type = parts[2]
                if feature_type not in ['gene', 'mRNA', 'CDS']:
                    continue

                chromosome = parts[0]
                attributes = parts[8]

                # Extract gene ID from attributes
                gene_id = None
                
                if is_busco:
                    # For BUSCO files, use name istead of ID
                    match = re.search(r'Name=([^;]+)', attributes)
                    if match:
                        gene_id = match.group(1)
                else:
                    # Try different attribute patterns
                    for pattern in [r'ID=([^;]+)', r'Name=([^;]+)', r'gene_id=([^;]+)']:
                        match = re.search(pattern, attributes)
                        if match:
                            gene_id = match.group(1)
                            break

                    # If no pattern matched, try splitting by semicolon
                    if not gene_id:
                        for attr in attributes.split(';'):
                            attr = attr.strip()
                            if '=' in attr:
                                key, value = attr.split('=', 1)
                                if key.lower() in ['id', 'name', 'gene_id']:
                                    gene_id = value
                                    break

                if gene_id:
                    # Clean gene ID
                    gene_id = gene_id.strip().strip('"').strip("'")
                    gene_id = gene_id.replace(":", "_")
                    gene_to_chr[gene_id] = (final_species_name, chromosome)
                    gene_count += 1

                    # Debug: log first few mappings
                    if gene_count <= 5:
                        debug_info.append(f"  Gene {gene_count}: {gene_id} -> {final_species_name}:{chromosome}")

    except Exception as e:
        print(f"[WARNING] Error processing {gff_file}: {e}")
        debug_info.append(f"ERROR processing {gff_file}: {e}")

    print(f"[INFO] Mapped {gene_count} genes from {final_species_name}")
    debug_info.append(f"Mapped {gene_count} genes from {final_species_name}")

print(f"[INFO] Total mapped genes: {len(gene_to_chr)}")
debug_info.append(f"Total mapped genes: {len(gene_to_chr)}")

# Debug: Check overlap between orthogroups and GFF genes
orthogroup_genes = set()
for _, row in orthogroups_df.iterrows():
    for species, genes in row.items():
        if species == "Orthogroup":
            continue
        if pd.isna(genes) or str(genes).strip() == "":
            continue
        for gene in str(genes).split(","):
            gene = gene.strip()
            orthogroup_genes.add(gene)

mapped_genes = set(gene_to_chr.keys())
overlap = orthogroup_genes.intersection(mapped_genes)

print(f"[INFO] Genes in orthogroups: {len(orthogroup_genes)}")
print(f"[INFO] Genes mapped from GFF: {len(mapped_genes)}")
print(f"[INFO] Overlap: {len(overlap)}")

debug_info.append(f"Genes in orthogroups: {len(orthogroup_genes)}")
debug_info.append(f"Genes mapped from GFF: {len(mapped_genes)}")
debug_info.append(f"Overlap: {len(overlap)}")

# Show sample mappings and orthogroup genes for comparison
if gene_to_chr:
    print("[INFO] Sample gene mappings from GFF:")
    for i, (gene, (sp, chr)) in enumerate(gene_to_chr.items()):
        if i < 10:
            print(f"  {gene} -> {sp}:{chr}")
            debug_info.append(f"  GFF gene: {gene} -> {sp}:{chr}")

if orthogroup_genes:
    print("[INFO] Sample genes from orthogroups:")
    for i, gene in enumerate(list(orthogroup_genes)[:10]):
        print(f"  {gene}")
        debug_info.append(f"  Orthogroup gene: {gene}")

# Track orthologous chromosome pairs
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

    # Get all species pairs and update orthologous chromosome counts
    species_list = list(gene_chroms.keys())
    if len(species_list) >= 2:
        processed_orthogroups += 1
        for i in range(len(species_list)):
            for j in range(i+1, len(species_list)):
                sp1, sp2 = species_list[i], species_list[j]
                for chr1 in gene_chroms[sp1]:
                    for chr2 in gene_chroms[sp2]:
                        key = tuple(sorted([ (sp1, chr1), (sp2, chr2) ]))
                        orthologous_chr_pairs[key] += 1
                        species_chr_syntenic[sp1].add(chr1)
                        species_chr_syntenic[sp2].add(chr2)

print(f"[INFO] Processed {processed_orthogroups} orthogroups with multi-species genes")
debug_info.append(f"Processed {processed_orthogroups} orthogroups with multi-species genes")

# Create output DataFrames
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

# Save results
pairwise_output.to_csv("pairwise_chromosome_orthology.tsv", sep='\\t', index=False)
species_chr_count.to_csv("species_orthologous_chromosomes.tsv", sep='\\t', index=False)

# Save debug info
with open("debug_gene_mapping.txt", "w") as f:
    f.write("\\n".join(debug_info))

print("\\n[INFO] Results saved:")
print("- pairwise_chromosome_orthology.tsv")
print("- species_orthologous_chromosomes.tsv")
print("- debug_gene_mapping.txt")

print("\\n[INFO] Pairwise chromosome orthology summary:")
print(pairwise_output.head())

print("\\n[INFO] Species orthologous chromosome counts:")
print(species_chr_count)
EOF

    cat <<-END_VERSIONS > versions.yml
"ECOFLOW_GENOMEQC:GENOMEQC:GENOME_AND_ANNOTATION:ORTHOLOGOUS_CHROMOSOMES":
    python: \$(python --version | sed 's/Python //g')
    pandas: \$(python -c "import pandas; print(pandas.__version__)")
END_VERSIONS
    """
}
