#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# -----------------------
# Load data
# -----------------------
fig_file = "/Users/noah/Desktop/ufc_repository/results/analysis_4f_prs_utility/Fig12.tsv"
pgs_file = "/Users/noah/Desktop/ufc_repository/results/analysis_4c_results/4c_prs_results/breast.PGS004688.analysis_4c.tsv"

df = pd.read_csv(fig_file, sep="\t")
pgs = pd.read_csv(pgs_file, sep="\t")

# convert to %
df["isolated_pct"] = df["isolated_pct_above"] * 100
df["familial_pct"] = df["familial_pct_above"] * 100
df["control_pct"] = df["control_pct_above"] * 100

# BRCA points only
df_brca = df[df["Gene"].isin(["BRCA1", "BRCA2"])]

# -----------------------
# Extract PGS groups
# -----------------------
control = pgs.loc[pgs["group"] == "Control", "PGS"].values
familial = pgs.loc[pgs["group"] == "Familial Breast", "PGS"].values
isolated = pgs.loc[pgs["group"] == "Isolated Breast", "PGS"].values

# -----------------------
# Summary stats
# -----------------------
def summarize(name, vals):
    mu = np.mean(vals)
    sd = np.std(vals)
    return mu, sd

mu_c, sd_c = summarize("Control", control)
mu_f, sd_f = summarize("Familial Breast", familial)
mu_i, sd_i = summarize("Isolated Breast", isolated)

print("Control",mu_c,sd_c)
print("Familial",mu_f,sd_f)
print("Isolated",mu_i,sd_i)
# -----------------------
# X range (zoomed window)
# -----------------------
x = np.arange(2.8, 3.55 + 0.1, 0.1)

# tail probability: P(X >= x)
def tail(mu, sd):
    return (1 - norm.cdf(x, loc=mu, scale=sd)) * 100

y_c = tail(mu_c, sd_c)
y_f = tail(mu_f, sd_f)
y_i = tail(mu_i, sd_i)

# -----------------------
# Plot
# -----------------------
fig, ax = plt.subplots(figsize=(3.2, 2.5))

# curves
ax.plot(x, y_f, color="#b24e44", lw=1.6)
ax.plot(x, y_c, color="#B0B0B0", lw=1.6)
ax.plot(x, y_i, color="#FF6F61", lw=1.6)

# points
ax.scatter(df_brca["SD_for_gene"], df_brca["familial_pct"],
           label="Breast Dx, Breast FHx +", s=12, color="#b24e44", zorder=3)

ax.scatter(df_brca["SD_for_gene"], df_brca["isolated_pct"],
           label="Breast Dx, Breast FHx -", s=12, color="#FF6F61", zorder=3)

ax.scatter(df_brca["SD_for_gene"], df_brca["control_pct"],
           label="Control", s=12, color="#B0B0B0", zorder=3)

# -----------------------
# Axes
# -----------------------
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_xlabel("Breast Cancer PRS Z-Score", fontsize=7)

ax.set_xlim(2.8, 3.55)
ax.set_ylim(bottom=0)

ax.set_xticks([2.9, 3.2, 3.5])
ax.tick_params(labelsize=7)

# -----------------------
# Secondary axis
# -----------------------
BASE_OR = 1.782

def sd_to_or(sd):
    return BASE_OR ** sd

def or_to_sd(or_val):
    return np.log(or_val) / np.log(BASE_OR)

secax = ax.secondary_xaxis('bottom', functions=(sd_to_or, or_to_sd))
secax.set_xlabel("Equivalent Odds Ratio", fontsize=7)
secax.spines["bottom"].set_position(("outward", 30))

gene_or = np.array([7.62, 5.23])
secax.set_xticks(gene_or)
#secax.set_xticklabels([f"{v:.2f}" for v in gene_or], fontsize=7)
secax.set_xticklabels([])
# -----------------------
# Legend + save
# -----------------------
ax.legend(loc="upper right", frameon=False, fontsize=7)

plt.tight_layout()
plt.savefig(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4f_prs_utility/Supplementary Figure 12b.pdf",
    dpi=600
)
plt.close()