#!/usr/bin/python3

# Written by Chris Wyatt and released under the MIT license. Converts a group of quast outputs to a table to plot on a tree

import pandas as pd
import argparse
import os
import re

def main():
    parser = argparse.ArgumentParser(description="Combine QUAST summary files, transpose the output, and specify plot types.")
    parser.add_argument('files', nargs='+', help="List of QUAST summary files to combine.")
    parser.add_argument('-o', '--output', required=True, help="Output CSV file.")
    parser.add_argument('-col', '--columns', help="Comma-separated list of columns to include in the final output.")
    parser.add_argument('-plot_types', '--plot_types', help="Comma-separated list of plot types for the columns.")
    args = parser.parse_args()

    combined_df = None

    for file in args.files:
        # Read each QUAST file into a DataFrame
        df = pd.read_csv(file, sep='\t', index_col=0, skiprows=0)

        # Extract the filename without the extension
        species_name = os.path.basename(file)

        # Remove any extension by using regex
        species_name_cleaned = re.sub(r'\..*$', '', species_name)

        # Rename the columns using the cleaned species name
        df.columns = [species_name_cleaned]

        if combined_df is None:
            combined_df = df
        else:
            combined_df = combined_df.join(df, how='outer')  # Outer join to combine based on row index

    # Transpose the DataFrame so that rows become columns and vice versa
    transposed_df = combined_df.T

    # If the -col flag is set, filter to include only the specified columns
    if args.columns:
        selected_columns = args.columns.split(',')
        transposed_df = transposed_df[selected_columns]

    # Create the plot types line
    if args.plot_types:
        plot_types_list = args.plot_types.split(',')
        if len(plot_types_list) != len(transposed_df.columns):
            raise ValueError("Number of plot types must match the number of columns in the output.")
    else:
        plot_types_list = ['text'] * len(transposed_df.columns)  # Default to 'text' for all columns if not provided

    # Insert the "species" column and NA for the second row
    header = ['species'] + list(transposed_df.columns)
    plot_type_row = ['NA'] + plot_types_list

    # Insert the plot types as the second line in the output file
    with open(args.output, 'w') as f_out:
        # Write the header (with species column)
        f_out.write("\t".join(header) + "\n")
        
        # Write the plot types row (with NA in the species column)
        f_out.write("\t".join(plot_type_row) + "\n")

        # Write the data (species column and the rest)
        transposed_df.to_csv(f_out, sep='\t', header=False)

if __name__ == "__main__":
    main()
