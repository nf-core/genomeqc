process BUSCO_TSV_TO_GFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/python:3.9--1"

    input:
    tuple val(meta), path(busco_dir)

    output:
    tuple val(meta), path("${meta.id}_busco.gff"), emit: gff
    tuple val(meta), path("${meta.id}_busco_stats.json"), emit: stats

    script:
    """
    #!/usr/bin/env python3

    import os
    import json
    import glob
    import csv
    from datetime import datetime

    def convert_busco_tsv_to_gff(busco_dir, sample_id):
        # Find the full_table.tsv file
        full_table_pattern = os.path.join(busco_dir, "**/full_table.tsv")
        full_table_files = glob.glob(full_table_pattern, recursive=True)

        if not full_table_files:
            raise FileNotFoundError("No full_table.tsv found in BUSCO output")

        full_table_file = full_table_files[0]
        output_gff = f"{sample_id}_busco.gff"
        stats_file = f"{sample_id}_busco_stats.json"

        print(f"Processing BUSCO table: {full_table_file}")

        stats = {
            "sample_id": sample_id,
            "timestamp": datetime.now().isoformat(),
            "input_file": os.path.basename(full_table_file),
            "total_buscos": 0,
            "complete_buscos": 0,
            "complete_single_copy": 0,
            "complete_duplicated": 0,
            "fragmented_buscos": 0,
            "missing_buscos": 0,
            "output_gff": output_gff
        }

        features = []

        # Use CSV reader to handle TSV properly
        with open(full_table_file, 'r') as inf:
            # Read first line to check if it's a header
            first_line = inf.readline().strip()
            inf.seek(0)  # Reset file pointer

            # Skip header if present
            if first_line.startswith('# Busco id') or first_line.startswith('Busco id'):
                next(inf)  # Skip header line
            else:
                inf.seek(0)  # Reset if no header

            reader = csv.reader(inf, delimiter='\\t')

            for row_num, row in enumerate(reader):
                if len(row) < 6:  # Need at least 6 columns
                    continue

                busco_id = row[0].strip()
                status = row[1].strip()
                sequence = row[2].strip()
                start_str = row[3].strip()
                end_str = row[4].strip()
                strand = row[5].strip() if len(row) > 5 else '.'
                score = row[6].strip() if len(row) > 6 else '.'
                length = row[7].strip() if len(row) > 7 else '.'
                description = row[9].strip() if len(row) > 9 else ''

                # Skip if this looks like a header row
                if busco_id == 'Busco id' or start_str == 'Gene Start':
                    continue

                stats["total_buscos"] += 1

                # Count by status
                if status == "Complete":
                    stats["complete_single_copy"] += 1
                    stats["complete_buscos"] += 1
                elif status == "Duplicated":
                    stats["complete_duplicated"] += 1
                    stats["complete_buscos"] += 1
                elif status == "Fragmented":
                    stats["fragmented_buscos"] += 1
                elif status == "Missing":
                    stats["missing_buscos"] += 1
                    continue  # Skip missing BUSCOs

                # Skip entries without valid coordinates
                if sequence in ["N/A", ""] or start_str in ["N/A", ""] or end_str in ["N/A", ""]:
                    continue

                # Validate coordinates are numeric
                try:
                    start_int = int(start_str)
                    end_int = int(end_str)
                except ValueError:
                    print(f"Warning: Skipping {busco_id} - invalid coordinates")
                    continue

                # Convert strand notation
                if strand in ["N/A", ""]:
                    strand = "."
                elif strand == "Plus":
                    strand = "+"
                elif strand == "Minus":
                    strand = "-"

                # CREATE COORDINATE-BASED GENE ID (like your BUSCO sequences)
                # Format: NC_000908.2:24468-26006|+
                coordinate_gene_id = f"{sequence}_{start_int}-{end_int}|{strand}"

                # Create attributes using coordinate-based ID
                attributes = f"ID={coordinate_gene_id};Name={coordinate_gene_id};BUSCO_ID={busco_id};Status={status};Type=BUSCO_gene"
                if score not in ["N/A", ""]:
                    attributes += f";Score={score}"
                if length not in ["N/A", ""]:
                    attributes += f";Length={length}"
                if description and description != "N/A":
                    clean_desc = description.replace(';', ',').replace('=', ':')
                    attributes += f";Description={clean_desc}"

                # Create GFF line with coordinate-based ID
                gff_line = [
                    sequence,
                    "BUSCO",
                    "gene",
                    str(start_int),
                    str(end_int),
                    score if score not in ["N/A", ""] else ".",
                    strand,
                    ".",
                    attributes
                ]

                features.append(gff_line)

        # Sort features
        if features:
            features.sort(key=lambda x: (x[0], int(x[3])))

        # Write GFF file
        with open(output_gff, 'w') as outf:
            outf.write("##gff-version 3\\n")
            outf.write(f"##source BUSCO analysis for {sample_id}\\n")
            outf.write(f"##date {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
            outf.write("##description BUSCO gene predictions with coordinate-based IDs\\n")

            for feature in features:
                outf.write('\\t'.join(feature) + '\\n')

        # Write stats
        with open(stats_file, 'w') as statsf:
            json.dump(stats, statsf, indent=2)

        print(f"Successfully converted {len(features)} BUSCO features to GFF")
        print(f"Stats: Complete={stats['complete_buscos']}, Fragmented={stats['fragmented_buscos']}, Missing={stats['missing_buscos']}")

        # Print some example coordinate IDs for verification
        if features:
            print("Sample coordinate-based gene IDs:")
            for i, feature in enumerate(features[:5]):
                coord_id = feature[8].split(';')[0].replace('ID=', '')
                print(f"  {coord_id}")

    # Run the conversion
    convert_busco_tsv_to_gff("${busco_dir}", "${meta.id}")
    """
}
