#!/usr/bin/python3

# Written by Chris Wyatt and released under the MIT license. 
# Converts a group of busco outputs to a table to plot on a tree

import pandas as pd
import argparse

# Set up the argument parser
parser = argparse.ArgumentParser(description='Extract and merge specific columns from a table.')
parser.add_argument('input_file', type=str, help='Path to the input TSV file.')
parser.add_argument('output_file', type=str, help='Path to save the output TSV file.')

# Parse the arguments
args = parser.parse_args()

# Read the input table into a pandas DataFrame
df = pd.read_csv(args.input_file, sep='\t')

# Select the required columns
df_extracted = df[['Input_file', 'Single', 'Duplicated', 'Fragmented', 'Missing']]

# Merge the columns from 'Complete' to 'Missing' into a single column, with values separated by commas
df_extracted['busco'] = df_extracted[['Single', 'Duplicated', 'Fragmented', 'Missing']].astype(str).agg(','.join, axis=1)

# Drop the individual 'Complete' to 'Missing' columns
df_extracted = df_extracted[['Input_file', 'busco']]

# Write the header and custom line first
with open(args.output_file, 'w') as f:
    # Write the header
    f.write('species\tbusco\n')
    # Insert 'NA<tab>stacked' as the second line
    f.write('NA\tpie\n')

# Append the DataFrame content to the file without the header
df_extracted.to_csv(args.output_file, sep='\t', index=False, mode='a', header=False)

print(f"Extraction completed successfully. Output saved to {args.output_file}.")
