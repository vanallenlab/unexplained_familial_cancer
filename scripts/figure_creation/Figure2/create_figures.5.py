#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# -------------------------
# Load data
# -------------------------
df = pd.read_csv(
    "/Users/noah/Desktop/ufc_repository/results/analysis_5_gsea_results/cosmic_genes_jan1.tsv",
    sep="\t"
)

# -------------------------
# Clean data
# -------------------------
df = df.copy()
df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
df["AF"] = pd.to_numeric(df["AF"], errors="coerce")

df = df[
    df["beta"].notna() &
    df["p_value"].notna() &
    np.isfinite(df["beta"]) &
    np.isfinite(df["p_value"]) &
    (df["p_value"] > 0)
].copy()

# -------------------------
# Extract Tier
# -------------------------
df["Tier"] = df["Pathogenic_Threshold"].str.extract(r"(Tier\d+)")

# -------------------------
# Compute effect sizes
# -------------------------
df["OR"] = np.exp(df["beta"])
df["logOR"] = np.log2(df["OR"])
df["neglog10p"] = -np.log10(df["p_value"])

# -------------------------
# Styling
# -------------------------
sns.set_style("white")
sns.set_context("paper", font_scale=1.2)

fig, ax = plt.subplots(figsize=(8, 6))

# -------------------------
# Aesthetic mappings
# -------------------------
base_blue = "#0D49E8"
tier_alpha = {
    "Tier2": 1,
    "Tier3": 0.75,
    "Tier4": 0.45,
    "Tier0": 0.2,
}

af_sizes = {
    0.001: 130,  # rare = bigger
    0.01: 40,    # more common = smaller
}

# -------------------------
# Plot points
# -------------------------
for (tier, af), sub in df.groupby(["Tier", "AF"]):
    ax.scatter(
        sub["logOR"],
        sub["neglog10p"],
        color=base_blue,
        alpha=tier_alpha[tier],
        s=af_sizes[af],
        edgecolor="black",
        linewidth=1,
        zorder=3,
        label=f"{tier} / AF={af}"  # for reference, but we’ll make separate legends
    )

# -------------------------
# Reference lines
# -------------------------
sig_threshold = 0.05
ax.axhline(-np.log10(sig_threshold), ls="--", lw=1.1, color="#333333", zorder=2)
ax.axhline(-np.log10(0.00037), ls="--", lw=1.1, color="#333333", zorder=3)
ax.axvline(0, color="#555555", lw=1, zorder=2)

# -------------------------
# Axes limits
# -------------------------
xlim = max(abs(df["logOR"].quantile(0.01)), abs(df["logOR"].quantile(0.99)))
xlim = min(max(xlim, 2), 6)
ax.set_xlim(-xlim, xlim)
ax.set_ylim(0, df["neglog10p"].max() + 2)

# -------------------------
# Axes appearance
# -------------------------
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
for spine in ["left", "bottom"]:
    ax.spines[spine].set_linewidth(1.2)
    ax.spines[spine].set_color("black")

ax.set_xlabel(r"$log_{2}(OR)$ (cases vs. matched controls)", fontsize=12, fontweight="bold")
ax.set_ylabel(r"$-\log_{10} P$", fontsize=12)
# ax.set_title(
#     "Prevalence of Potentially PGVs in CPGs across 17 Cancers",
#     fontsize=13,
#     fontweight="bold",
#     pad=10
# )

# -------------------------
# Legends
# -------------------------
# 1) Tier (alpha)
tier_handles = [
    plt.Line2D(
        [0], [0],
        marker='o',
        color=base_blue,
        alpha=a,
        markersize=8,
        linestyle='None',
        markeredgecolor='black',
        markeredgewidth=0.4
    )
    for t, a in tier_alpha.items()
]
tier_labels = [t.replace("Tier0", "VUS") for t in tier_alpha.keys()]
tier_labels = [t for t in tier_labels if t != "VUS"] + ["VUS"]

legend1 = ax.legend(
    tier_handles, tier_labels,
    title="Tier Shade",
    frameon=False,
    fontsize=9,
    title_fontsize=10,
    loc="upper right"
)

# 2) AF (size)
af_handles = [
    plt.Line2D(
        [0], [0],
        marker='o',
        color='gray',
        markersize=np.sqrt(s),  # matplotlib marker size ~ sqrt(points)
        linestyle='None',
        markeredgecolor='black',
        markeredgewidth=0.4
    )
    for af, s in af_sizes.items()
]
af_labels = [f"AF ≤ {af * 100}%" for af in af_sizes.keys()]
legend2 = ax.legend(
    af_handles, af_labels,
    title="Allele Frequency",
    frameon=False,
    fontsize=9,
    title_fontsize=10,
    loc="upper left"
)

# Attach both legends
ax.add_artist(legend1)
ax.add_artist(legend2)

# -------------------------
# Tight layout
# -------------------------
plt.tight_layout(rect=[0, 0, 0.9, 1])

# -------------------------
# Save
# -------------------------
plt.savefig(
    "/Users/noah/Desktop/ufc_repository/results/analysis_5_gsea_results/Figure_2D.png",
    dpi=700,
    bbox_inches="tight"
)

print("Saved: volcano_plot_tier_af.png")
