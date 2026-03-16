#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

# Files
snp_file = "/Users/noah/Desktop/ufc_repository/results/supplementary_figures/step_5_count_variants/denovo_snp.counts.tsv"
indel_file = "/Users/noah/Desktop/ufc_repository/results/supplementary_figures/step_5_count_variants/denovo_indel.counts.tsv"

# Colors (same palette you've been using)
snv_color = "#D2B48C"
indel_color = "#8B5A2B"

# Load data
snp = pd.read_csv(snp_file, sep="\t")
indel = pd.read_csv(indel_file, sep="\t")

# Nature-style formatting
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out"
})

fig, ax = plt.subplots(figsize=(3.5, 2.5))

# Plot medians (Q2)
ax.plot(
    snp["TP_PROB"],
    snp["Q2"],
    color=snv_color,
    linewidth=1.5,
    label="De novo SNVs"
)

ax.plot(
    indel["TP_PROB"],
    indel["Q2"],
    color=indel_color,
    linewidth=1.5,
    label="De novo indels"
)

# QC cutoffs
ax.axvline(x=0.63, color=snv_color, linestyle="--", linewidth=1)
ax.axvline(x=0.76, color=indel_color, linestyle="--", linewidth=1)

ax.text(0.63, ax.get_ylim()[1]*0.95, "SNV QC cutoff",
        rotation=90, ha="right", va="top", fontsize=6, color=snv_color)

ax.text(0.76, ax.get_ylim()[1]*0.95, "Indel QC cutoff",
        rotation=90, ha="right", va="top", fontsize=6, color=indel_color)

# Labels
ax.set_xlabel("True positive probability threshold")
ax.set_ylabel("Median de novo variants per sample")

# Nature-style axes
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(loc="center right",frameon=False, fontsize=5, handlelength=0.5)

plt.tight_layout()

# Save
out_file = "/Users/noah/Desktop/ufc_repository/results/supplementary_figures/random_forest_figs/rf_denovo_counts.pdf"
plt.savefig(out_file)

print("Saved:", out_file)