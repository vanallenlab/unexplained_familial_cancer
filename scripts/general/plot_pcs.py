import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import os


def plot_pcs_part1(data, id_col, color_by, num_pcs, output_name):
    # Set a general style
    sns.set(style="whitegrid")

    # Load and merge data
    df = pd.read_csv(data, sep='\t', index_col=False)
    #df["#IID"] = df["#IID"].astype(str)

    #df2 = pd.read_csv("dfci-ufc.sample_meta.gatkhc_posthoc_outliers.tsv", sep='\t', index_col=False)
    #df = df.merge(df2, left_on="#IID", right_on="original_id")
    #df = df[df['intake_qc_pop'] == 'EUR']

    # Create a figure with 3 rows and 3 columns (9 plots for PC1–PC10 pairs)
    fig, axes = plt.subplots(3, 3, figsize=(18, 12))
    axes = axes.flatten()

    # Loop over adjacent PC pairs
    for i in range(num_pcs - 1):
        x_pc = f'PC{i+1}'
        y_pc = f'PC{i+2}'
        ax = axes[i]
        
        sns.scatterplot(
            data=df, x=x_pc, y=y_pc, hue=color_by, palette='tab10', s=20, ax=ax, alpha=0.8
        )
        ax.set_title(f'{x_pc} vs {y_pc}')
        ax.legend(loc='best', fontsize='small')
        ax.set_xlabel(x_pc)
        ax.set_ylabel(y_pc)

    # Remove any unused subplot (none in this case)
    # Adjust layout
    plt.tight_layout()
    plt.savefig(f"{output_name}.part1.png")

def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-d', '--data', required=True, help = 'File to use that has PC information and Sample information.')
    parser.add_argument('-i','--id-col', required=False, default = "#IID", help = "The name of the column to use for Sample ID.")
    parser.add_argument('-c', '--color-by',required=False, help = "Column to color by in the data.")
    parser.add_argument('-n', '--num-pcs', required=False, default = 10, type=int, help = "How many PCs to plot.")
    parser.add_argument('-o', '--output-name', required=False, default = "pc_plot", help = "Prefix to save files to.")
    args = parser.parse_args()

    plot_pcs_part1(args.data, args.id_col, args.color_by, args.num_pcs, args.output_name)

if __name__ == "__main__":
    main()