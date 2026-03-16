#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import os

# ----------------------------
# Inputs / outputs
# ----------------------------
in_matrix = "/Users/noah/Desktop/ufc_repository/results/Figure_2/fig2A_data.tsv"
out_pdf = "/Users/noah/Desktop/ufc_repository/results/Figure_2/Figure2A.pdf"
out_png = "/Users/noah/Desktop/ufc_repository/results/Figure_2/Figure2A.png"
# ----------------------------
# Load matrix (gene x cancer, binary)
# ----------------------------
heat = pd.read_csv(in_matrix, sep="\t", index_col=0)
heat = heat.loc[heat.index.isin(["BRCA1","ATM","BAP1","MSH2","TINF2"])]
heat["Cervix"] = 0
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
# Reorder x-axis
# ----------------------------
priority_order = ["Breast", "Prostate", "Colorect", "Colorectal", "Uterus", "Ovary", "Kidney", "Melanoma", "Thyroid"]
# Keep only columns present in heat
existing_cols = [c for c in priority_order if c in heat.columns]

# Remaining columns not in priority_order (excluding Control)
remaining_cols = sorted([c for c in heat.columns if c not in existing_cols + ["Control"]])

# Combine final order
final_order = existing_cols + remaining_cols
if "Control" in heat.columns:
    final_order.append("Control")

# Reorder heat columns
heat = heat[final_order]

# ----------------------------
# Plot heatmap
# ----------------------------
fig_width = max(4.5, len(heat.columns) * 0.18)  # slightly narrower
fig_height = max(1.8, len(heat) * 0.35)

plt.figure(figsize=(2.8, 1.75))

ax = sns.heatmap(
    heat,
    cmap=["#f0f0f0", "#0D49E880"],
    linewidths=0.25,
    linecolor="gray",
    cbar=False
)

# ----------------------------
# Remove y-label
# ----------------------------
ax.set_ylabel("")  # no y-axis label

# ----------------------------
# Axis formatting
# ----------------------------
ax.set_xticks([i + 0.5 for i in range(len(heat.columns))])
ax.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=5)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=5, fontstyle="italic")

# ----------------------------
# Manual annotation dictionaries
# ----------------------------
COSMIC_PAIRS = {
    "BRCA1": ["Breast"],
    "BAP1": ["Kidney", "Melanoma"],
    "MSH2": ["Colorectal", "Uterus", "Ovary"]
}

LITERATURE_PAIRS = {
    "BRCA1": ["Prostate"],
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
plt.savefig(out_png,bbox_inches="tight",pad_inches=0.02)


print("Saved heatmap PDF to:", out_pdf)
