#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------
# Load data
# -----------------------
file = "/Users/noah/Desktop/ufc_repository/results/analysis_4f_prs_utility/Fig12.tsv"
df = pd.read_csv(file, sep="\t")

df["isolated_pct"] = df["isolated_pct_above"] * 100
df["familial_pct"] = df["familial_pct_above"] * 100
df["control_pct"] = df["control_pct_above"] * 100

# -----------------------
# Plot
# -----------------------
fig, ax = plt.subplots(figsize=(3.2, 2.5))

ax.scatter(df["SD_for_gene"], df["isolated_pct"],
           label="Isolated", s=4, alpha=0.7, color="#FF6F61")

ax.scatter(df["SD_for_gene"], df["familial_pct"],
           label="Familial", s=4, alpha=0.9, color="#b24e44")

ax.scatter(df["SD_for_gene"], df["control_pct"],
           label="Control", s=4, alpha=0.8, color="#B0B0B0")

# -----------------------
# Clean axes
# -----------------------
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_xlabel("Breast Cancer PRS Z-Score", fontsize=7)
ax.set_ylabel("% with PRS exceeding\npathogenic variant effect size", fontsize=7)
ax.tick_params(labelsize=7)
ax.set_ylim(bottom=0)

# =========================================================
# FIXED SECONDARY AXIS (IMPORTANT CHANGE)
# =========================================================

BASE_OR = 1.782

# safety wrapper (prevents runtime warning)
def sd_to_or(sd):
    sd = np.asarray(sd)
    return BASE_OR ** sd

def or_to_sd(or_val):
    or_val = np.asarray(or_val)
    or_val = np.where(or_val <= 0, np.nan, or_val)  # prevents log(0)
    return np.log(or_val) / np.log(BASE_OR)

# CREATE TRUE secondary axis (NO "twiny")
secax = ax.secondary_xaxis('bottom', functions=(sd_to_or, or_to_sd))

secax.set_xlabel("Equivalent Odds Ratio", fontsize=7)

# push it below main axis so BOTH are visible
secax.spines["bottom"].set_position(("outward", 30))

# -----------------------
# Gene OR markers
# -----------------------
gene_or = np.array([
    1.82, 1.37, 7.62, 5.23, 2.50,
    2.47, 1.93, 3.83, 1.20, 1.72
])

gene_or = gene_or[np.isfinite(gene_or) & (gene_or > 0)]

secax.set_xticks(gene_or)
#secax.set_xticklabels([f"{x:.2f}" for x in gene_or], fontsize=6)
secax.set_xticklabels([])
# -----------------------
# Legend
# -----------------------
#ax.legend(loc="upper center",frameon=False, fontsize=7)

plt.tight_layout()

# -----------------------
# Save
# -----------------------
plt.savefig(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4f_prs_utility/Supplementary Figure 12a.pdf",
    dpi=600
)
plt.close()