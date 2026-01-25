#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.patches import Patch

# -----------------------
# Color parsing
# -----------------------
def parse_color_file(filepath):
    color_map = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            k, v = line.split(":", 1)
            color_map[k.lower()] = v.strip().strip('"')
    return color_map

cancer_color = parse_color_file(
    "/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml"
)

# -----------------------
# Inputs
# -----------------------
or_table = pd.read_csv(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.OR.tsv",
    sep="\t"
)
pval_table = pd.read_csv(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.pvalue.tsv",
    sep="\t"
)

pgs_dir = Path("/Users/noah/Desktop/ufc_repository/results/analysis_4a_prs/")

cancer_configs = {
    "breast": ["PGS004688"],
    "cervix": ["PGS000784"],
    "colorectal": ["PGS000785"],
    "kidney": ["PGS000787"],
    "melanoma": ["PGS003382"],
    "non-hodgkin": ["PGS000791"],
    "prostate": ["PGS004694"],
    "thyroid": ["PGS000797"],
    "hematologic": ["PGS000788"],
    "lung": ["PGS000789"],
    "ovary": ["PGS004249"],
    "uterus": ["PGS004244"],
    "brain": ["PGS003384"],
    "bladder": ["PGS004687"],
}

# -----------------------
# Style (tight, Arial)
# -----------------------
plt.rcParams.update({
    "font.family": "Arial",
    "axes.titlesize": 5,
    "axes.labelsize": 5,
    "xtick.labelsize": 4.5,
    "ytick.labelsize": 4.5,
    "axes.linewidth": 0.4,
})

# -----------------------
# Figure
# -----------------------
fig, axes = plt.subplots(
    2, 7,
    figsize=(7.5, 1.5),
    sharex=True,
    sharey=True
)

axes = axes.flatten()

# -----------------------
# Plot loop
# -----------------------
for i, (cancer, pgs_ids) in enumerate(cancer_configs.items()):
    pgs_id = pgs_ids[0]
    ax = axes[i]

    # OR and p-value
    or_val = float(
        or_table.loc[or_table.PGS_ID == pgs_id, cancer].values[0].split()[0]
    )
    p_val = float(
        pval_table.loc[pval_table.PGS_ID == pgs_id, cancer].values[0]
    )

    pgs_file = list(pgs_dir.glob(f"{cancer}.{pgs_id}*.anon_pgs"))[0]
    df = pd.read_csv(pgs_file, sep="\t")

    palette = {
        1: cancer_color[cancer],  # case
        0: "#B0B0B0"               # control
    }

    # IMPORTANT: normalize densities per group
    sns.kdeplot(
        data=df,
        x="PGS",
        hue="case",
        ax=ax,
        fill=True,
        alpha=0.7,
        bw_adjust=0.8,
        linewidth=0,
        palette=palette,
        legend=False,
        common_norm=False   # 🔑 fixes density interpretation
    )

    # Annotation
    ax.text(
        0.98, 0.92,
        f"Or={or_val:.2f}\np={p_val:.1e}",
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=4.5,
        bbox=dict(
            facecolor="white",
            alpha=0.6,
            pad=0.15,
            linewidth=0
        )
    )

    ax.set_title(cancer.replace("_", " ").capitalize(), pad=1)

    # Axis labels (minimal)
    if i % 7 == 0:
        ax.set_ylabel("Density", labelpad=1)
    else:
        ax.set_ylabel("")

    if i >= 7:
        ax.set_xlabel("PGS", labelpad=1)
    else:
        ax.set_xlabel("")

    ax.tick_params(length=2, width=0.8)

# -----------------------
# Remove unused axes
# -----------------------
for j in range(len(cancer_configs), len(axes)):
    fig.delaxes(axes[j])

# -----------------------
# Legend (control only)
# -----------------------
legend_elements = [
    Patch(facecolor="#B0B0B0", label="Control")
]

fig.legend(
    handles=legend_elements,
    loc="center right",
    bbox_to_anchor=(1.005, 0.55),
    frameon=False,
    fontsize=5,
    title_fontsize=5
)

# -----------------------
# Layout tuning
# -----------------------
plt.subplots_adjust(
    left=0.035,
    right=0.935,
    top=0.95,
    bottom=0.15,
    wspace=0.1,
    hspace=0.3
)

# -----------------------
# Save
# -----------------------
plt.savefig(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4a_prs/Figure_4A.png",
    dpi=600
)
plt.savefig(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4a_prs/Figure_4A.pdf",
    dpi=600
)
plt.close()
