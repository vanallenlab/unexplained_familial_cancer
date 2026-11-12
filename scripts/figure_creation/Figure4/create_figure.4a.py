import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Manually editable odds ratios
or_data = {
     "CCDC27": [1.06, 0.89, 1.24, 1.98, 3],
    "BRCA1":  [8.48, 9.4, 11.71, 12.16],
    "BRCA2":  [6.05, 6.6, 5.68, 5.96],
    "PALB2":  [4.10, 3.74, 3.59, 5.13],
    "CHEK2":  [2.17, 1.15, 1.82, 1.4],
    "ATM":    [2.04, 2.91, 4.3, 4.96],
    "CDH1": [1.72,6.98, 34.55], 
    "MSH2":   [0.9, 1.15, 1.4, 1.65],
}

fdr_labels = ["≥0", "≥1", "≥2", "≥3", "≥3,\nno PGV"]
colors = ["#F8C9C7", "#F08078", "#E52B27", "#A00008", "#4DA4EA"]

genes = list(or_data.keys())

fig, ax = plt.subplots(figsize=(3.5, 2))

# Each gene occupies one fixed y-position
y_positions = {
    gene: len(genes) - 1 - i
    for i, gene in enumerate(genes)
}

for gene in genes:
    y = y_positions[gene]
    odds_ratios = or_data[gene]

    # Horizontal connecting lines: no vertical incline
    for i in range(len(odds_ratios) - 1):
        ax.plot(
            [odds_ratios[i], odds_ratios[i + 1]],
            [y, y],
            color=colors[i + 1],
            linewidth=2.5,
            alpha=0.45,
            zorder=1
        )

    # Four dots, all on the same horizontal row
    for i, odds_ratio in enumerate(odds_ratios):
        ax.scatter(
            odds_ratio,
            y,
            color=colors[i],
            s=30,
            edgecolor="none",
            zorder=2
        )

# Reference line for no association
ax.axvline(
    x=1,
    color="#999999",
    linestyle="--",
    linewidth=1.5,
    zorder=0
)

# Axes and labels
ax.set_yticks([y_positions[gene] for gene in genes])
ax.set_yticklabels(genes, fontsize=5, fontstyle="italic")

for label in ax.get_yticklabels():
    if label.get_text() == "CCDC27":
        label.set_fontweight("bold")
        
ax.set_xscale("log")
ax.set_xlim(0.8, 37)
ax.set_xticks([1, 3, 10, 35])
ax.set_xticklabels(["1", "3", "10", "35"], fontsize=7)

ax.set_xlabel('Breast Cancer Odds Ratios for "Tier 1" Variants in CPGs', fontsize=7, labelpad=10)
ax.set_ylabel("")
ax.set_ylim(-0.7, len(genes) - 0.3)

# Legend
legend_handles = [
    Line2D(
        [0], [0],
        marker="o",
        color="none",
        markerfacecolor=colors[i],
        markeredgecolor="none",
        markersize=9,
        label=fdr_labels[i]
    )
    for i in range(5)
]

ax.legend(
    handles=legend_handles,
    title="Relatives\nwith\ncancer",
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    frameon=True,
    fontsize=7,
    title_fontsize=7
)

# Clean styling
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(False)

plt.tight_layout()
plt.savefig("/Users/noah/Desktop/ufc_repository/results/Figure3/Figure_3a.pdf", dpi=300, bbox_inches="tight")
# plt.show()