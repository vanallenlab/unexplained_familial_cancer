#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

# Input files
snp_file = "/Users/noah/Desktop/ufc_repository/results/supplementary_figures/step_6_sensitivity/sensitivity_output_snp.tsv"
indel_file = "/Users/noah/Desktop/ufc_repository/results/supplementary_figures/step_6_sensitivity/sensitivity_output_indel.tsv"

# Load data
snp = pd.read_csv(snp_file, sep="\t")
indel = pd.read_csv(indel_file, sep="\t")

# Rename SNP label to SNV
snp["variant_type"] = "SNVs"
indel["variant_type"] = "Indels"

# Colors (tan + darker brown)
colors = {
    "SNVs": "#D2B48C",   # tan
    "Indels": "#8B5A2B"  # darker brown
}

# Nature-style plotting
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

ax.plot(
    snp["tp_prob"],
    snp["sensitivity"] * 100,
    color=colors["SNVs"],
    linewidth=1.5,
    label="SNVs"
)

ax.plot(
    indel["tp_prob"],
    indel["sensitivity"],
    color=colors["Indels"],
    linewidth=1.5,
    label="Indels"
)

# Labels
ax.set_xlabel("True positive probability threshold")
ax.set_ylabel("Sensitivity (%)")

# Style to match Nature journals
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(frameon=False)

ax.axvline(x=0.63, color="black", linestyle="--", linewidth=1)
ax.axvline(x=0.76, color="black", linestyle="--", linewidth=1)

ax.text(0.63, ax.get_ylim()[1]*0.5, "SNV QC cutoff",
        rotation=90, va="top", ha="right", fontsize=6, color=colors["SNVs"])

ax.text(0.76, ax.get_ylim()[1]*0.5, "Indel QC cutoff",
        rotation=90, va="top", ha="right", fontsize=6, color=colors["Indels"])

plt.tight_layout()

# Save
out_file = "/Users/noah/Desktop/ufc_repository/results/supplementary_figures/random_forest_figs/rf_sensitivity.pdf"
plt.savefig(out_file)

print(f"Saved: {out_file}")