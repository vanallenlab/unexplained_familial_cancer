#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import roc_curve, roc_auc_score


def main():
    p = argparse.ArgumentParser(description="Plot AUROC curve in Nature Genetics-style formatting.")
    p.add_argument("--tsv", required=True, help="Input TSV with columns: cancer_status, percent_homozygosed")
    p.add_argument("--out", required=True, help="Output PDF/PNG path (recommended: PDF)")
    p.add_argument("--figsize", nargs=2, type=float, default=[1.8, 1.8], metavar=("W", "H"),
                   help="Figure size in inches, e.g. --figsize 2.2 2.2")
    p.add_argument("--color", default="#98DF8A", help="ROC curve color (hex or named), default green")
    p.add_argument("--title", default=None, help="Optional title (Nat Gen usually prefers none)")
    args = p.parse_args()

    # -----------------------------
    # Load + clean
    # -----------------------------
    df = pd.read_csv(args.tsv, sep="\t", dtype=str)

    needed = {"cancer_status", "percent_homozygosed"}
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit(f"ERROR: missing columns in TSV: {sorted(missing)}")

    # handle "." and other junk safely
    df["cancer_status"] = pd.to_numeric(df["cancer_status"], errors="coerce")
    df["percent_homozygosed"] = pd.to_numeric(df["percent_homozygosed"], errors="coerce")
    df = df.dropna(subset=["cancer_status", "percent_homozygosed"]).copy()

    # enforce int labels 0/1
    df["cancer_status"] = df["cancer_status"].astype(int)

    if df["cancer_status"].nunique() < 2:
        raise SystemExit("ERROR: need both cases (1) and controls (0) to compute AUROC.")

    y = df["cancer_status"].values
    s = df["percent_homozygosed"].values

    # -----------------------------
    # AUROC + ROC curve
    # -----------------------------
    auc = roc_auc_score(y, s)
    fpr, tpr, _ = roc_curve(y, s)

    # -----------------------------
    # Nat Gen-ish formatting
    # -----------------------------
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 1,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=tuple(args.figsize))

    # diagonal
    ax.plot([0, 1], [0, 1], linewidth=0.8, linestyle="--", color="black", alpha=0.6)

    # roc curve
    ax.plot(fpr, tpr, linewidth=1.6, color=args.color)

    # axis formatting
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")

    # clean spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ticks
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_yticks([0, 0.5, 1.0])

    # AUROC label on plot
    ax.text(
        0.60, 0.08,
        f"AUROC = {auc:.2f}",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=7
    )

    if args.title:
        ax.set_title(args.title, fontsize=7)

    fig.tight_layout(pad=0.4)
    fig.savefig(args.out, bbox_inches="tight", facecolor="white",pad_inches=0)
    plt.close(fig)

    print(f"Saved: {args.out}")
    print(f"AUROC: {auc:.6f}")


#python3 create_figure.5d.py --tsv /Users/noah/Desktop/ufc_repository/results/analysis_1_roh/kidney_roh_results/kidney_roh.tsv  --out /Users/noah/Desktop/ufc_repository/results/analysis_1_roh/kidney_roh_results/kidney_roh.auroc.pdf



if __name__ == "__main__":
    main()
