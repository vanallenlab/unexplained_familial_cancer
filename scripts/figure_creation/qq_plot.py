from qqman import qqman
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

def return_pvals(file, max_maf, annotation, column_of_interest="recalibrated_p"):
    df = pd.read_csv(file, sep='\t')
    df = df[(df['criteria'] == annotation) & (df['max_AF'] == max_maf)]
    if column_of_interest not in df:
        return []
    return df[column_of_interest].dropna().tolist()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data', required=True)
    parser.add_argument('-o', '--output', required=False, default='saige_results.png')
    args = parser.parse_args()

    df = pd.read_csv(args.data, sep='\t')

    max_mafs = [0.01, 0.001]
    groups = sorted(df['criteria'].unique())

    nrows = len(max_mafs)
    ncols = len(groups)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols, 3*nrows))
    axes = np.atleast_2d(axes)

    for i, af in enumerate(max_mafs):
        for j, grouping in enumerate(groups):

            p_vals = np.array(
                return_pvals(args.data, af, grouping), 
                dtype=float
            )

            if len(p_vals) == 0 or np.all(np.isnan(p_vals)):
                axes[i, j].set_title(f"{grouping} <= {af*100}% (no data)")
                axes[i, j].plot([], [])
                continue

            qqman.qqplot(p_vals, ax=axes[i, j])
            axes[i, j].set_title(f"{grouping} <= {af*100}%")

    plt.tight_layout()
    plt.savefig(args.output)
    plt.close()



if __name__ == "__main__":
	main()