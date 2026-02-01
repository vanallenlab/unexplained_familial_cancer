#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker

# -------------------------
# Load data
# -------------------------
df = pd.read_csv(
    "/Users/noah/Desktop/ufc_repository/results/analysis_5_gsea_results/cosmic_genes_jan21.tsv.gz",
    sep="\t"
)
top_hits = df.loc[df.groupby("cancer")["p_value"].idxmin()]

# -------------------------
# Clean data
# -------------------------
df = df.copy()
df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
df["AF"] = pd.to_numeric(df["AF"], errors="coerce")
df = df[~df['Pathogenic_Threshold'].str.contains("Tier5")]

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

fig, ax = plt.subplots(figsize=(4, 4.25))

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
    0.001: 65,  # rare = bigger
    0.01: 20,    # more common = smaller
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
# Highlight top association per cancer
# -------------------------
top_hits = df.loc[df.groupby("cancer")["p_value"].idxmin()]

ax.scatter(
    top_hits["logOR"],
    top_hits["neglog10p"],
    s=3,                 # small center dot
    c="red",
    edgecolors="none",
    zorder=4
)

# -------------------------
# Reference lines
# -------------------------
sig_threshold = 0.05

ax.axhline(-np.log10(0.05), ls="--", lw=1.1, color="#333333", zorder=3)
ax.axhline(-np.log10(0.00037), ls="--", lw=1.1, color="#333333", zorder=3)
ax.axvline(0, color="#555555", lw=1, zorder=2)

# -------------------------
# Axes limits
# -------------------------
xlim = max(abs(df["logOR"].quantile(0.01)), abs(df["logOR"].quantile(0.99)))
xlim = min(max(xlim, 2), 6)
ax.set_xlim(-xlim, xlim)
ax.set_ylim(0, df["neglog10p"].max() + 2)

ax.set_ylim(bottom=0)
ax.yaxis.set_major_locator(mticker.MultipleLocator(1.5))
# -------------------------
# Axes appearance
# -------------------------
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
for spine in ["left", "bottom"]:
    ax.spines[spine].set_linewidth(1.2)
    ax.spines[spine].set_color("black")

ax.set_xlabel(r"$log_{2}(OR)$ (cases vs. matched controls)", fontsize=7, fontweight="bold")
ax.set_ylabel(r"$-\log_{10} P$", fontsize=7)
ax.set_title(
    "Enrichment of rare non-synonymous variants in CPGs across 17 cancers",
    fontsize=7,
    fontweight="bold"
)

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

# legend1 = ax.legend(
#     tier_handles, tier_labels,
#     title="Tier Shade",
#     frameon=False,
#     fontsize=7,
#     title_fontsize=7,
#     loc="upper right"
# )

legend1 = ax.legend(
    tier_handles, tier_labels,
    title="Tier Shade",
    frameon=False,
    fontsize=7,
    title_fontsize=7,
    loc="upper left",                 # anchor point of the legend box
    bbox_to_anchor=(0.8, 0.98)       # (x, y) in axes fraction coords
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
        markeredgewidth=0.8
    )
    for af, s in af_sizes.items()
]
af_labels = [f"AF ≤ {af * 100}%" for af in af_sizes.keys()]
# legend2 = ax.legend(
#     af_handles, af_labels,
#     title="Allele Frequency",
#     frameon=False,
#     fontsize=7,
#     title_fontsize=7,
#     loc="upper left"
# )
legend2 = ax.legend(
    af_handles, af_labels,
    title="Allele Frequency",
    frameon=False,
    fontsize=7,
    title_fontsize=7,
    bbox_to_anchor=(0.27, 0.93)
)

# 3) Top association per cancer (red center dot)
top_hit_handle = plt.Line2D(
    [0], [0],
    marker='o',
    color='red',
    markersize=4,        # small, matches center dot
    linestyle='None',
    markeredgecolor='none'
)

legend3 = ax.legend(
    handles=[top_hit_handle],
    labels=["Most significant association\n(per cohort)"],
    frameon=False,
    fontsize=6,
    loc="lower right",
    title_fontsize=5
)

# Attach both legends
ax.add_artist(legend1)
ax.add_artist(legend2)
ax.add_artist(legend3)

# -------------------------
# Tight layout
# -------------------------
plt.xlim(-4,6)
plt.ylim(0,3.5)
plt.tight_layout(rect=[0, 0, 1, 1])

# -------------------------
# Save
# -------------------------
plt.savefig(
    "/Users/noah/Desktop/ufc_repository/results/analysis_5_gsea_results/Figure_2D.pdf",
    dpi=700,
    bbox_inches="tight",pad_inches=0
)

print("Saved: volcano_plot_tier_af.pdf")
