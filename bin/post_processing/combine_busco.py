import os
import pandas as pd
import argparse

def extract_busco_data(filepath):
    """
    Extracts the relevant BUSCO metrics from the input file and returns it as a DataFrame.
    """
    data_dict = {}

    with open(filepath, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if line.startswith("C:"):
            # Extract percentages (Complete, Fragmented, Missing)
            metrics = line.split(',')
            data_dict['Complete BUSCOs (C)'] = metrics[0].split(':')[1].split('%')[0].strip()
            data_dict['Fragmented BUSCOs (F)'] = metrics[1].split(':')[1].split('%')[0].strip()
            data_dict['Missing BUSCOs (M)'] = metrics[2].split(':')[1].split('%')[0].strip()
        elif "Complete BUSCOs (C)" in line:
            data_dict['Complete BUSCOs'] = line.split()[0].strip()
        elif "Complete and single-copy BUSCOs (S)" in line:
            data_dict['Complete and single-copy BUSCOs'] = line.split()[0].strip()
        elif "Complete and duplicated BUSCOs (D)" in line:
            data_dict['Complete and duplicated BUSCOs'] = line.split()[0].strip()
        elif "Fragmented BUSCOs (F)" in line:
            data_dict['Fragmented BUSCOs'] = line.split()[0].strip()
        elif "Missing BUSCOs (M)" in line:
            data_dict['Missing BUSCOs'] = line.split()[0].strip()
        elif "Total BUSCO groups searched" in line:
            data_dict['Total BUSCO groups searched'] = line.split()[0].strip()

    # Convert the dictionary to a DataFrame
    df = pd.DataFrame.from_dict(data_dict, orient='index', columns=[os.path.basename(filepath).split('.')[0]])
    
    return df

def combine_files(file_list, output_filepath):
    # Initialize an empty DataFrame to store the consolidated data
    consolidated_df = pd.DataFrame()

    # Loop through each file in the provided file list
    for filepath in file_list:
        # Extract the BUSCO data as a DataFrame
        df = extract_busco_data(filepath)
        
        # Ensure unique column names by appending the full basename
        column_prefix = os.path.basename(filepath).replace('.', '_')
        df.columns = [f"{column_prefix}_{col}" for col in df.columns]
        
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
    parser = argparse.ArgumentParser(description="Extract BUSCO data from multiple files and combine them into a single CSV file.")
    parser.add_argument('files', metavar='F', type=str, nargs='+', help='Input BUSCO summary files (e.g., *.txt)')
    parser.add_argument('-o', '--output', type=str, default='consolidated_busco_table.csv', help='Output CSV file name')

    # Parse the arguments
    args = parser.parse_args()

    # Combine files and save to the output file
    combine_files(args.files, args.output)

