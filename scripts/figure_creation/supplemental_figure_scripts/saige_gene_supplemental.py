#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description="Create QQ plot from SAIGE gene stats")
    parser.add_argument("--cancer", required=True, help="Cancer type (e.g., breast)")
    parser.add_argument("--criteria", required=True, help="Criteria filter (e.g., T1;T2)")
    parser.add_argument("--max_AF", required=True, help="Max AF filter (001 or 0001)")
    parser.add_argument("--genes", required=False, default="", help="Comma separated genes to highlight")
    parser.add_argument("--color", required=False, default="#4d4d4d", help="Hex color for plot")
    parser.add_argument("--out", default="qq_plot.pdf", help="Output file")
    return parser.parse_args()


def convert_af(code):
    if code == "001":
        return 0.01
    elif code == "0001":
        return 0.001
    else:
        try:
            return float(code)
        except:
            raise ValueError("max_AF must be 001, 0001, or numeric")


def load_data(cancer):
    path = f"/Users/noah/Desktop/ufc_repository/results/SAIGE_GENE/{cancer}.saige.recalibrated.stats.tsv.gz"
    df = pd.read_csv(path, sep="\t", compression="gzip")
    return df


def make_qq(df, highlight_genes, color, out, cancer, criteria, af):
    
    pvals = df["recalibrated_p"].dropna()
    pvals = pvals[pvals > 0]

    observed = -np.log10(np.sort(pvals))
    expected = -np.log10(np.linspace(1/len(pvals), 1, len(pvals)))

    fig, ax = plt.subplots(figsize=(2,2))

    ax.scatter(expected, observed, s=6, color=color, alpha=0.7)

    maxv = max(expected.max(), observed.max())
    ax.plot([0, maxv], [0, maxv], linestyle="--", linewidth=1)
    ax.axhline(y=-np.log10(2.69e-6), linestyle="--", linewidth=1, color="red")
    ax.text(-0.25, -np.log10(2.69e-6) - 0.4, "Genome-wide significance", color="red", fontsize=5, va="bottom")
    # Highlight genes
    if highlight_genes:
        sub = df[df["#gene"].isin(highlight_genes)]
        for _, row in sub.iterrows():
            p = row["recalibrated_p"]
            if p <= 0:
                continue

            obs = -np.log10(p)

            rank = (pvals <= p).sum()
            exp = -np.log10(rank/len(pvals))

            ax.scatter(exp, obs, s=25,color="royalblue")

            ax.annotate(
                row["#gene"],
                xy=(exp, obs),
                xytext=(exp+0.2, obs+0.2),
                arrowprops=dict(arrowstyle="->", lw=0.8),
                fontsize=7,
                fontstyle='italic'
            )

    ax.set_xlabel("Expected -log10(P)",fontsize=5)
    ax.set_ylabel("Observed -log10(P)",fontsize=5)
    ax.tick_params(axis='both', labelsize=5)
    ax.set_title(f"{cancer.capitalize().replace("_patient_and_family_prs", " wtih 2+ FDRs").replace("_cell", " Cell")} {criteria.split(';')[0]}-{criteria.split(';')[-1].replace("T1-T1","T1")} (AF<{float(af) * 100}%)",fontsize=7)

    # Remove top/right spines
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(f"/Users/noah/Desktop/ufc_repository/results/supplementary_figures/saige_supplement/{cancer}.{criteria.split(';')[-1]}.{af}.pdf",dpi=300,pad_inches=0,bbox_inches="tight")
    plt.close()


def main():

    args = parse_args()

    af = convert_af(args.max_AF)
    highlight_genes = [g.strip() for g in args.genes.split(",") if g]

    df = load_data(args.cancer)

    df = df[df["criteria"] == args.criteria]
    df = df[df["max_AF"] == af]

    make_qq(df, highlight_genes, args.color, args.out, args.cancer, args.criteria, af)


if __name__ == "__main__":
    main()