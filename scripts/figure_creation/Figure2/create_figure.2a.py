#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import os

# ----------------------------
# Inputs / outputs
# ----------------------------
in_matrix = "/Users/noah/Desktop/ufc_repository/results/Figure_2/lof_heatmap_matrix.tsv"
out_pdf = "/Users/noah/Desktop/ufc_repository/results/Figure_2/lof_heatmap_selected_genes.pdf"

# ----------------------------
# Load matrix (gene x cancer, binary)
# ----------------------------
heat = pd.read_csv(in_matrix, sep="\t", index_col=0)

# ----------------------------
# Rename cancers
# ----------------------------
rename_map = {
    "Basal_Cell_Carcinoma": "BCC",
    "Squamous_Cell_Carcinoma": "SCC",
    "Non-Hodgkin": "NHL",
    "control": "Control"
}
heat = heat.rename(columns=rename_map)

# (optional) also clean underscores -> spaces for remaining cancers
heat.columns = [c.replace("_", " ") for c in heat.columns]

# ----------------------------
# Plot heatmap (same styling)
# ----------------------------
fig_width = max(5, len(heat.columns) * 0.22)
fig_height = max(1.8, len(heat) * 0.35)

plt.figure(figsize=(3.5, 1.75))

ax = sns.heatmap(
    heat,
    cmap=["#f0f0f0", "#0D49E880"],
    linewidths=0.25,
    linecolor="gray",
    cbar=False
)

# ----------------------------
# Manual annotation dictionaries
# ----------------------------
COSMIC_PAIRS = {
    "BRCA1": ["Breast", "Prostate"],
    "BAP1": ["Kidney", "Melanoma"],
    "MSH2": ["Colorectal", "Uterus", "Ovary"]
}

LITERATURE_PAIRS = {
    "ATM": ["Breast"],
    "TINF2": ["Thyroid"],
    "STK11": ["Breast"]
}

# ----------------------------
# Overlay COSMIC-associated pairs (red solid)
# ----------------------------
for gene, cancers in COSMIC_PAIRS.items():
    if gene not in heat.index:
        continue
    y = list(heat.index).index(gene)
    for cancer in cancers:
        if cancer not in heat.columns:
            continue
        x = list(heat.columns).index(cancer)
        ax.add_patch(Rectangle((x, y), 1, 1, fill=False, edgecolor="red", linewidth=1.2))

# ----------------------------
# Overlay literature-associated pairs (orange dashed)
# ----------------------------
for gene, cancers in LITERATURE_PAIRS.items():
    if gene not in heat.index:
        continue
    y = list(heat.index).index(gene)
    for cancer in cancers:
        if cancer not in heat.columns:
            continue
        x = list(heat.columns).index(cancer)
        ax.add_patch(Rectangle((x, y), 1, 1, fill=False, edgecolor="orange", linewidth=1.2, linestyle="--"))

# ----------------------------
# Legend
# ----------------------------
legend_elements = [
    Patch(facecolor="#0D49E880", edgecolor="black", label="≥1 individual"),
    Patch(facecolor="none", edgecolor="red", linewidth=1.2, label="COSMIC-associated pair"),
    Patch(facecolor="none", edgecolor="orange", linewidth=1.2, linestyle="--", label="Literature-associated pair")
]

# ax.legend(
#     handles=legend_elements,
#     loc="upper left",
#     bbox_to_anchor=(1.01, 1),
#     frameon=False,
#     fontsize=5
# )
ax.legend(
    handles=legend_elements,
    bbox_to_anchor=(0.25, -0.4),
    frameon=False,
    fontsize=5
)
# ----------------------------
# Axis formatting
# ----------------------------
ax.set_xticks([i + 0.5 for i in range(len(heat.columns))])
ax.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=5)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=5, fontstyle="italic")

plt.tight_layout(rect=[0,0,1.3,1.3])

# ----------------------------
# Save as PDF
# ----------------------------
#os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
plt.savefig(out_pdf, format="pdf", bbox_inches="tight",pad_inches=0)


print("Saved heatmap PDF to:", out_pdf)
