import pandas as pd
import argparse
from pathlib import Path

# Example usage: python3 combine_results.py *.stats -o BRAT1_all_results.tsv

def combine_results(input_files, output_file):
    dfs = [pd.read_csv(file, sep="\t") for file in input_files]

    combined = pd.concat(dfs, ignore_index=True, sort=False)

    combined.to_csv(output_file, sep="\t", index=False)

    print(f"Combined {len(input_files)} files")
    print(f"Total rows: {len(combined):,}")
    print(f"Total columns: {len(combined.columns)}")
    print(f"Output: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine logistic regression/Fisher TSV results.")
    parser.add_argument("input_files", nargs="+", help="Input TSV files")
    parser.add_argument("-o", "--output", default="combined_results.tsv", help="Output TSV filename")

    args = parser.parse_args()

    combine_results(args.input_files, args.output)
