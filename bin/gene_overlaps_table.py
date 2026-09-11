#!/usr/bin/env python3

# Written by Fernando Duarte with AI assistance, released under the MIT license.
# Extracts gene overlap statistics from per-sample files into a combined table

import pandas as pd
import argparse

# Set up the argument parser
parser = argparse.ArgumentParser(description='Extract gene statistics from files.')
parser.add_argument('input_files', nargs='+', help='List of input files.')
parser.add_argument('output_file', help='Path to save the output TSV file.')
parser.add_argument('--include-same-strand', action='store_true', help='Include same-strand genes count.')
parser.add_argument('--include-opposite-strand', action='store_true', help='Include opposite-strand genes count.')

# Parse the arguments
args = parser.parse_args()

# Initialize an empty list to store the results
results = []

# Process each input file
for file in args.input_files:
    try:
        # Load the file into a DataFrame
        df = pd.read_csv(file, sep='\t', header=None, names=['Statistic', 'Count'])

        # Extract required statistics
        total_genes = df.loc[df['Statistic'] == 'Total number of genes', 'Count'].values[0]
        overlapping_genes = df.loc[df['Statistic'] == 'Total number of overlapping genes', 'Count'].values[0]

        # Optional statistics
        same_strand_genes = df.loc[df['Statistic'] == 'Number of genes fully contained in same strand direction', 'Count'].values[0] if args.include_same_strand else "NA"
        opposite_strand_genes = df.loc[df['Statistic'] == 'Number of genes fully contained in opposite strand direction', 'Count'].values[0] if args.include_opposite_strand else "NA"
        print(same_strand_genes)

        # Collect results in a dictionary
        entry = {
            'File': file,
            'Total_genes': total_genes,
            'Overlapping_genes': overlapping_genes,
        }
        if args.include_same_strand:
            entry['Fully_contained_same_strand_genes'] = same_strand_genes
        if args.include_opposite_strand:
            entry['Fully_contained_opposite_strand_genes'] = opposite_strand_genes

        results.append(entry)
    except Exception as e:
        print(f"Error processing {file}: {e}")
        continue

# Convert the results to a DataFrame
columns = ['File', 'Total_genes', 'Overlapping_genes']
if args.include_same_strand:
    columns.append('Fully_contained_same_strand_genes')
if args.include_opposite_strand:
    columns.append('Fully_contained_opposite_strand_genes')

result_df = pd.DataFrame(results, columns=columns)

# Write the result to the output file
result_df.to_csv(args.output_file, sep='\t', index=False)
print(f"Extraction completed successfully. Output saved to {args.output_file}.")
