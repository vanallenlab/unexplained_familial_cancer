#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------
# Load gene summary data (points)
# -----------------------
file = "/Users/noah/Desktop/ufc_repository/results/analysis_4f_prs_utility/Fig12.tsv"
df = pd.read_csv(file, sep="\t")

df["isolated_pct"] = df["isolated_pct_above"] * 100
df["familial_pct"] = df["familial_pct_above"] * 100
df["control_pct"] = df["control_pct_above"] * 100
df = df[(df['Gene'] != "NF1")]
df = df[(df['Gene'] != "RAD51D")]
df = df[(df['Gene'] != "CDH1")]

# -----------------------
# Load raw PGS data (for curves)
# -----------------------
pgs_file = "/Users/noah/Desktop/ufc_repository/results/analysis_4c_results/4c_prs_results/breast.PGS004688.analysis_4c.tsv"
pgs = pd.read_csv(pgs_file, sep="\t")

# -----------------------
# Define thresholds (x-axis)
# -----------------------
thresholds = np.linspace(pgs["PGS"].min(), pgs["PGS"].max(), 200)

# -----------------------
# Function to compute tail percentages
# -----------------------
def compute_curve(values, thresholds):
    return np.array([
        (values >= t).mean() * 100
        for t in thresholds
    ])

# -----------------------
# Split groups
# -----------------------
control_vals = pgs.loc[pgs["group"] == "Control", "PGS"].values
familial_vals = pgs.loc[pgs["group"] == "Familial Breast", "PGS"].values

# If isolated exists, include it safely
if "Isolated Breast" in pgs["group"].unique():
    isolated_vals = pgs.loc[pgs["group"] == "Isolated Breast", "PGS"].values
else:
    isolated_vals = np.array([])

# -----------------------
# Compute curves
# -----------------------
control_curve = compute_curve(control_vals, thresholds)
familial_curve = compute_curve(familial_vals, thresholds)

if len(isolated_vals) > 0:
    isolated_curve = compute_curve(isolated_vals, thresholds)

# -----------------------
# Plot
# -----------------------
fig, ax = plt.subplots(figsize=(3.2, 2.5))

# Curves
ax.plot(thresholds, control_curve, color="#B0B0B0", lw=1, alpha=0.9)
ax.plot(thresholds, familial_curve, color="#b24e44", lw=1.2, alpha=0.9)

if len(isolated_vals) > 0:
    ax.plot(thresholds, isolated_curve, color="#FF6F61", lw=1, alpha=0.9)

# Points (original)
ax.scatter(df["SD_for_gene"], df["control_pct"],
           s=6, color="#B0B0B0", zorder=3)

ax.scatter(df["SD_for_gene"], df["familial_pct"],
           s=6, color="#b24e44", zorder=3)

if "isolated_pct" in df.columns:
    ax.scatter(df["SD_for_gene"], df["isolated_pct"],
               s=6, color="#FF6F61", zorder=3)

# -----------------------
# Clean axes
# -----------------------
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_xlabel("Breast Cancer PRS Z-Score", fontsize=7)
ax.set_ylabel("% with PRS exceeding\npathogenic variant effect size", fontsize=7)

ax.tick_params(labelsize=7)
ax.set_ylim(bottom=0)
ax.set_xlim(0.5,3.55)
# =========================================================
# Secondary axis (OR scale)
# =========================================================

BASE_OR = 1.782

def sd_to_or(sd):
    return BASE_OR ** sd

def or_to_sd(or_val):
    or_val = np.asarray(or_val)
    or_val = np.where(or_val <= 0, np.nan, or_val)
    return np.log(or_val) / np.log(BASE_OR)

secax = ax.secondary_xaxis('bottom', functions=(sd_to_or, or_to_sd))
secax.set_xlabel("Equivalent Odds Ratio", fontsize=7)
secax.spines["bottom"].set_position(("outward", 30))

# Gene OR ticks
gene_or = np.array([
    1.82, 1.37, 7.62, 5.23, 2.47, 3.83, 1.20,
])

gene_or = gene_or[np.isfinite(gene_or) & (gene_or > 0)]
secax.set_xticks(gene_or)
secax.set_xticklabels([])

# -----------------------
# Save
# -----------------------
plt.tight_layout()

plt.savefig(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4f_prs_utility/Supplementary Figure 12a.pdf",
    dpi=600
)

plt.close()