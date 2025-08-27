#!/usr/bin/env python3

# This script takes the output n buscos per sequence
# table from plot_busco_ideogram.R and concatenates
# the results from different samples/species into a
# single table

# Credits to ChatGPT

import argparse
import os
import pandas as pd

def count_sequences_above_threshold(file, threshold):
    df = pd.read_csv(file)
    count = (df['Complete_BUSCOs'] > threshold).sum()
    return count

def main():
    parser = argparse.ArgumentParser(
        description="Count number of sequences with Complete_BUSCOs above a threshold."
    )
    parser.add_argument(
        "-i", "--input", nargs='+', required=True,
        help="Input CSV files."
    )
    parser.add_argument(
        "-x", "--threshold", type=int, default=1,
        help="Threshold for Complete_BUSCOs."
    )
    parser.add_argument(
        "-o", "--output", default="n_seqs_above_x_buscos.tsv",
        help="Output file name (default: output.tsv)."
    )

    args = parser.parse_args()

    results = []
    for file in args.input:
        base_name = os.path.splitext(os.path.basename(file))[0]
        count = count_sequences_above_threshold(file, args.threshold)
        results.append([base_name, count])

    output_df = pd.DataFrame(results, columns=["Id", f"Num_Seqs_Above_{args.threshold}"])
    output_df.to_csv(args.output, sep="\t", index=False)
    print(f"Results written to {args.output}")

if __name__ == "__main__":
    main()
