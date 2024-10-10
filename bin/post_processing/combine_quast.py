import os
import pandas as pd
import argparse

def combine_files(file_list, output_filepath):
    # Initialize an empty DataFrame to store the consolidated data
    consolidated_df = pd.DataFrame()

    # Loop through each file in the provided file list
    for filepath in file_list:
        # Get the filename without the directory path
        filename = os.path.basename(filepath)
        
        # Read the file into a DataFrame, skipping the first line (header)
        df = pd.read_csv(filepath, sep="\t", skiprows=1, header=None, index_col=0)
        
        # Use the file's name (without extension) as the column header
        column_name = os.path.splitext(filename)[0]
        df.columns = [column_name]
        
        # Merge the current file's DataFrame with the consolidated DataFrame
        if consolidated_df.empty:
            consolidated_df = df
        else:
            consolidated_df = consolidated_df.join(df, how='outer')

    # Reset the index to have the metric names as the first column
    consolidated_df.reset_index(inplace=True)

    # Rename the first column to 'Metric'
    consolidated_df.rename(columns={'index': 'Metric'}, inplace=True)

    # Save the consolidated table to a CSV file
    consolidated_df.to_csv(output_filepath, index=False)

    print(f"Consolidated data saved to {output_filepath}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Combine multiple QUAST summary files into a single CSV file.")
    parser.add_argument('files', metavar='F', type=str, nargs='+', help='Input QUAST summary files (e.g., *.tsv)')
    parser.add_argument('-o', '--output', type=str, default='consolidated_table.csv', help='Output CSV file name')

    # Parse the arguments
    args = parser.parse_args()

    # Combine files and save to the output file
    combine_files(args.files, args.output)

