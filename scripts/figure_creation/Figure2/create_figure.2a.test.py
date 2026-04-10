#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import matplotlib as mpl
mpl.rcParams["font.family"] = "Arial"

# ----------------------------
# Inputs / outputs
# ----------------------------
in_matrix = "/Users/noah/Desktop/ufc_repository/results/Figure_2/fig2A_data.tsv"
out_pdf = "/Users/noah/Desktop/ufc_repository/results/Figure_2/Figure2A.pdf"
out_png = "/Users/noah/Desktop/ufc_repository/results/Figure_2/Figure2A.png"

# ----------------------------
# Load and Clean Data
# ----------------------------
# Note: Ensure the file exists at the path above
heat = pd.read_csv(in_matrix, sep="\t", index_col=0)

# Filter for specific genes
genes_of_interest = ["BRCA1", "ATM", "BAP1", "MSH2", "TINF2"]
heat = heat.loc[heat.index.isin(genes_of_interest)]

# Add Cervix column if missing
if "Cervix" not in heat.columns:
    heat["Cervix"] = 0

rename_map = {
    "Basal_Cell_Carcinoma": "BCC",
    "Squamous_Cell_Carcinoma": "SCC",
    "Non-Hodgkin": "NHL",
    "control": "Control",
    "Neuroendocrine": "NETs"
}
heat = heat.rename(columns=rename_map)
heat.columns = [c.replace("_", " ") for c in heat.columns]

# ----------------------------
# Order columns
# ----------------------------
priority_order = [
    "Breast", "Prostate", "Colorectal", "Uterus", 
    "Ovary", "Kidney", "Melanoma", "Thyroid"
]

existing = [c for c in priority_order if c in heat.columns]
remaining = sorted([c for c in heat.columns if c not in existing and c != "Control"])

final_cols = existing + remaining
if "Control" in heat.columns:
    final_cols.append("Control")

heat = heat[final_cols]

# ----------------------------
# FIGURE SETUP
# ----------------------------
# Increased height slightly to accommodate legend comfortably
fig, ax = plt.subplots(figsize=(3.1, 1.7))

sns.heatmap(
    heat,
    ax=ax,
    cmap=["#f0f0f0", "#0D49E880"],
    linewidths=0.5,
    linecolor="gray",
    cbar=False
)

# ----------------------------
# TICKS
# ----------------------------
ax.set_xticks([i + 0.5 for i in range(len(heat.columns))])
ax.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=6)

ax.set_yticks([i + 0.5 for i in range(len(heat.index))])
ax.set_yticklabels(heat.index, fontsize=6, fontstyle="italic")
ax.set_ylabel("")

# ----------------------------
# ANNOTATIONS (Rectangles)
# ----------------------------
COSMIC_PAIRS = {
    "BRCA1": ["Breast"],
    "BAP1": ["Kidney", "Melanoma"],
    "MSH2": ["Colorectal", "Uterus", "Ovary"]
}

LITERATURE_PAIRS = {
    "BRCA1": ["Prostate"],
    "ATM": ["Breast"],
    "TINF2": ["Thyroid"]
}

def add_highlights(pairs, color, linestyle):
    for g, categories in pairs.items():
        if g in heat.index:
            y_idx = list(heat.index).index(g)
            for c in categories:
                if c in heat.columns:
                    x_idx = list(heat.columns).index(c)
                    ax.add_patch(Rectangle((x_idx, y_idx), 1, 1, 
                                           fill=False, edgecolor=color, 
                                           lw=1.2, linestyle=linestyle, clip_on=False))

add_highlights(COSMIC_PAIRS, "red", "-")
add_highlights(LITERATURE_PAIRS, "orange", "--")

# ----------------------------
# THE LEGEND
# ----------------------------
legend_elements = [
    Patch(facecolor="#0D49E880", edgecolor="black", label="≥1 LoF SV carrier"),
    Patch(facecolor="none", edgecolor="red", linewidth=1.2, label="COSMIC-associated pair"),
    Patch(facecolor="none", edgecolor="orange", linewidth=1.2, linestyle="--", label="Literature-associated pair")
]

# loc="upper center" and bbox_to_anchor=(0.5, -0.4) places it below the X-axis labels
ax.legend(
    handles=legend_elements,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.5), 
    frameon=False,
    fontsize=5,
    ncol=1  # Stack them vertically for clarity
)

# ----------------------------
# FINAL LAYOUT & SAVE
# ----------------------------
# Use tight_layout but leave room at the bottom for the legend
plt.tight_layout()
# Adjusting the bottom specifically to prevent the legend from being cut off
plt.subplots_adjust(bottom=0.35)

plt.savefig(out_pdf, format="pdf", bbox_inches="tight")
plt.savefig(out_png, dpi=300, bbox_inches="tight")

print(f"Success! Saved to: {out_pdf}")