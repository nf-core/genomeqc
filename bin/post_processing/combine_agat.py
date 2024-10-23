import os
import pandas as pd
import argparse

def extract_mrna_section(filepath):
    """
    Extracts the mrna section from the input file and returns it as a DataFrame.
    """
    with open(filepath, 'r') as file:
        lines = file.readlines()
    
    # Find the starting line for the mrna section
    start_idx = None
    for i, line in enumerate(lines):
        if "mrna" in line.lower() and "-----" in line:
            start_idx = i + 1
            break

    if start_idx is None:
        raise ValueError(f"mrna section not found in file: {filepath}")
    
    # Collect lines until an empty line or another section starts
    mrna_data = []
    for line in lines[start_idx:]:
        if line.strip() == "" or ("-----" in line):
            break
        mrna_data.append(line.strip())
    
    # Parse the extracted data into a dictionary
    data_dict = {}
    for line in mrna_data:
        if line:
            key, value = line.rsplit(maxsplit=1)
            data_dict[key.strip()] = value.strip()

    # Convert the dictionary to a DataFrame
    df = pd.DataFrame.from_dict(data_dict, orient='index', columns=[os.path.basename(filepath).split('.')[0]])
    
    return df

def combine_files(file_list, output_filepath):
    # Initialize an empty DataFrame to store the consolidated data
    consolidated_df = pd.DataFrame()

    # Loop through each file in the provided file list
    for filepath in file_list:
        # Extract the mrna section data as a DataFrame
        df = extract_mrna_section(filepath)
        
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
    parser = argparse.ArgumentParser(description="Extract mrna section data from multiple files and combine them into a single CSV file.")
    parser.add_argument('files', metavar='F', type=str, nargs='+', help='Input files to process (e.g., *.txt)')
    parser.add_argument('-o', '--output', type=str, default='consolidated_mrna_table.csv', help='Output CSV file name')

    # Parse the arguments
    args = parser.parse_args()

    # Combine files and save to the output file
    combine_files(args.files, args.output)

