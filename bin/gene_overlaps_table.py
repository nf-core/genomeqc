#!/usr/bin/python3

# Written using ChatGPT

import pandas as pd
import argparse

# Set up the argument parser
parser = argparse.ArgumentParser(description='Extract gene statistics from files.')
parser.add_argument('input_files', nargs='+', help='List of input files.')
parser.add_argument('output_file', help='Path to save the output TSV file.')
parser.add_argument('--include-sense', action='store_true', help='Include sense genes count.')
parser.add_argument('--include-antisense', action='store_true', help='Include antisense genes count.')

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
        sense_genes = df.loc[df['Statistic'] == 'Number of genes fully contained in sense direction', 'Count'].values[0] if args.include_sense else "NA"
        antisense_genes = df.loc[df['Statistic'] == 'Number of genes fully contained in antisense direction', 'Count'].values[0] if args.include_antisense else "NA"
        print(sense_genes)

        # Collect results in a dictionary
        entry = {
            'File': file,
            'Total_genes': total_genes,
            'Overlapping_genes': overlapping_genes,
        }
        if args.include_sense:
            entry['Fully_contained_sense_genes'] = sense_genes
        if args.include_antisense:
            entry['Fully_contained_antisense_genes'] = antisense_genes
        
        results.append(entry)
    except Exception as e:
        print(f"Error processing {file}: {e}")
        continue

# Convert the results to a DataFrame
columns = ['File', 'Total_genes', 'Overlapping_genes']
if args.include_sense:
    columns.append('Fully_contained_sense_genes')
if args.include_antisense:
    columns.append('Fully_contained_antisense_genes')

result_df = pd.DataFrame(results, columns=columns)

# Write the result to the output file
result_df.to_csv(args.output_file, sep='\t', index=False)
print(f"Extraction completed successfully. Output saved to {args.output_file}.")