#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# -------------------------------
# Global matplotlib settings
# -------------------------------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 0.4,
    "xtick.major.width": 0.4,
    "ytick.major.width": 0.4,
})

# -------------------------------
# 1. Load files
# -------------------------------
file_no_prs = "/Users/noah/Desktop/ufc_repository/results/epidemiological_results/patient_family_logistic_results.tsv"
file_with_prs = "/Users/noah/Desktop/ufc_repository/results/epidemiological_results/patient_family_logistic_with_prs.tsv"

df_no_prs = pd.read_csv(file_no_prs, sep="\t")
df_no_prs = df_no_prs[df_no_prs['n_intersection'] >= 5]
df_with_prs = pd.read_csv(file_with_prs, sep="\t")
df_with_prs = df_with_prs[df_with_prs['n_intersection'] >= 5]

# -------------------------------
# 2. Case-insensitive matching
# -------------------------------
for df in (df_no_prs, df_with_prs):
    df["patient_cancer"] = df["patient_cancer"].str.lower()
    df["family_cancer"] = df["family_cancer"].str.lower()

# -------------------------------
# 3. Normalize family cancer
# -------------------------------
def normalize_family(row):
    pc = row["patient_cancer"]
    fc = row["family_cancer"]

    if fc == "skin" and pc in {
        "melanoma",
        "basal_cell_carcinoma",
        "squamous_cell_carcinoma",
    }:
        return pc

    if fc == "blood_soft_tissue" and pc in {
        "hematologic",
        "non-hodgkin",
    }:
        return pc

    return fc

for df in (df_no_prs, df_with_prs):
    df["family_cancer"] = df.apply(normalize_family, axis=1)

# -------------------------------
# 4. Keep only matching cancers
# -------------------------------
df_no_prs = df_no_prs[
    df_no_prs["patient_cancer"] == df_no_prs["family_cancer"]
].copy()

df_with_prs = df_with_prs[
    df_with_prs["patient_cancer"] == df_with_prs["family_cancer"]
].copy()

# -------------------------------
# 5. Merge datasets
# -------------------------------
merged = df_no_prs.merge(
    df_with_prs,
    on=["patient_cancer", "family_cancer"],
    suffixes=("_no_prs", "_with_prs"),
)

# -------------------------------
# 6. Volcano plot function
# -------------------------------
def plot_volcano(
    merged,
    outfile,
    show_arrows=True,
    show_labels=True,
):
    fig, ax = plt.subplots(figsize=(3.5, 3.5))

    blue = "#1f77b4"
    green = "#2ca02c"

    for _, row in merged.iterrows():
        x0 = np.log2(row["odds_ratio_no_prs"])
        y0 = -np.log10(row["p_value_no_prs"])
        x1 = np.log2(row["odds_ratio_with_prs"])
        y1 = -np.log10(row["p_value_with_prs"])

        if show_arrows:
            ax.annotate(
                "",
                xy=(x1, y1),
                xytext=(x0, y0),
                arrowprops=dict(
                    arrowstyle="->",
                    linewidth=0.4,
                    color=blue,
                ),
                zorder=3,
            )

        ax.scatter(
            x0, y0,
            s=14,
            color=blue,
            linewidth=0.4,
            zorder=4,
        )

        ax.scatter(
            x1, y1,
            s=14,
            color=green,
            linewidth=0.4,
            zorder=5,
        )

        if show_labels:
            ax.text(
                x1,
                y1,
                row["patient_cancer"].capitalize(),
                fontsize=6,
                ha="left",
                va="bottom",
            )

    # Reference lines
    ax.axhline(-np.log10(0.05), linestyle="--", linewidth=0.4, color="black")
    ax.axhline(
        -np.log10(0.05 / len(merged)),
        linestyle="--",
        linewidth=0.4,
        color="black",
    )
    ax.axvline(0, linestyle="--", linewidth=0.4, color="black")

    # Labels
    ax.set_xlabel(r"$\log_2(\mathrm{OR})$ ", labelpad=2)
    ax.set_ylabel(r"$-\log_{10}P$", labelpad=2)
    

    # Legend
    legend_elements = [
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=blue, markersize=4,
               label="Without PRS"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=green, markersize=4,
               label="With PRS"),
    ]

    ax.legend(
        handles=legend_elements,
        frameon=False,
        loc="upper right",
        fontsize=6,
    )

    # Spine cleanup
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.title("Patient Cancer ~ Family Cancer + Sex + PRS (Cancer Specific)", fontsize=7)
    plt.tight_layout(pad=0.6)
    plt.savefig(outfile, dpi=300)
    plt.close()

# -------------------------------
# 7. Save outputs (exact paths)
# -------------------------------
plot_volcano(
    merged,
    outfile="/Users/noah/Desktop/ufc_repository/results/epidemiological_results/Figure_3E_labels.png",
    show_arrows=True,
    show_labels=True,
)

plot_volcano(
    merged,
    outfile="/Users/noah/Desktop/ufc_repository/results/epidemiological_results/Figure_3E_no_labels.png",
    show_arrows=False,
    show_labels=False,
)

plot_volcano(
    merged,
    outfile="/Users/noah/Desktop/ufc_repository/results/epidemiological_results/Figure_3E_no_labels.pdf",
    show_arrows=False,
    show_labels=False,
)
