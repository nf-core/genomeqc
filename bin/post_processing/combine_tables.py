import argparse

def combine_csv_files(file_list, output_filepath):
    with open(output_filepath, 'w') as outfile:
        for i, filepath in enumerate(file_list):
            with open(filepath, 'r') as infile:
                if i == 0:
                    # Write the first file entirely including its header
                    outfile.write(infile.read())
                else:
                    # Skip the header for subsequent files
                    next(infile)
                    outfile.write(infile.read())

    print(f"Combined data saved to {output_filepath}")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Combine multiple CSV files with the same header into one CSV file.")
    parser.add_argument('files', metavar='F', type=str, nargs='+', help='Input CSV files (e.g., *.csv)')
    parser.add_argument('-o', '--output', type=str, default='combined_output.csv', help='Output CSV file name')

    # Parse the arguments
    args = parser.parse_args()

    # Combine the CSV files and save to the output file
    combine_csv_files(args.files, args.output)

