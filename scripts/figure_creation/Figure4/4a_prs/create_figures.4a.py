#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.patches import Patch
import matplotlib.colors as mcolors

# -----------------------
# Color parsing
# -----------------------
def parse_color_file(filepath: str) -> dict:
    """Parse color YAML and return lowercase name → hex mapping."""
    color_map = {}
    with open(filepath, "r") as infile:
        for line in infile:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" in line:
                key, value = line.split(":", 1)
                color_map[key.lower()] = value.strip().strip('"').strip("'")
    return color_map

cancer_color = parse_color_file("/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml")

# -----------------------
# Darken color function
# -----------------------
def darken_color(color, factor=0.7):
    """Darken a matplotlib color by a factor (0 < factor < 1)."""
    rgb = mcolors.to_rgb(color)
    return tuple(channel * factor for channel in rgb)

# -----------------------
# Inputs
# -----------------------
or_table = pd.read_csv(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.OR.tsv",
    sep="\t", index_col=False
)
pval_table = pd.read_csv(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.pvalue.tsv",
    sep="\t", index_col=False
)
pgs_dir = Path("/Users/noah/Desktop/ufc_repository/results/analysis_4a_prs/")

# -----------------------
# Cancer PGS IDs
# -----------------------
cancer_configs = {
    "breast": ["PGS004688"],
    "cervix": ["PGS000784"],
    "colorectal": ["PGS000785"],
    "kidney": ["PGS000787"],
    "melanoma": ["PGS003382"],
    "non-hodgkin": ["PGS000791"],
    "prostate": ["PGS004694"],
    "thyroid": ["PGS000797"],
    "hematologic":["PGS000788"],
    "lung":["PGS000789"],
    "ovary":["PGS004249"],
    "uterus": ["PGS004244"],
    "brain": ["PGS003384"],
    "bladder":["PGS004687"]
}

# -----------------------
# Style
# -----------------------
plt.rcParams.update({
    'axes.labelsize': 5,
    'axes.labelweight': 'bold',
    'xtick.labelsize': 4.5,
    'ytick.labelsize': 4.5,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'axes.titlesize': 5,
    'axes.titleweight': 'bold',
})

# -----------------------
# Figure (wider + shorter)
# -----------------------
fig, axes = plt.subplots(
    2, 7, figsize=(7.5, 1.5), sharex=True, sharey=True
)
axes = axes.T.flatten()  # fill left-to-right

# OR/p-value offset from top-right of subplot (adjustable)
PVAL_X = 0.03  # relative to axes
OR_X = 0.95
OR_PVAL_Y = 0.95

# -----------------------
# Plot loop
# -----------------------
plot_idx = 0
for cancer, pgs_ids in cancer_configs.items():
    for pgs_id in pgs_ids:
        if not pgs_id:
            plot_idx += 1
            continue

        # OR / p-value
        or_val = float(or_table.loc[or_table['PGS_ID'] == pgs_id, cancer].values[0].split()[0])
        p_val = float(pval_table.loc[pval_table['PGS_ID'] == pgs_id, cancer])

        # Load PGS file
        pgs_file = list(pgs_dir.glob(f"{cancer}.{pgs_id}*.anon_pgs"))
        if not pgs_file:
            print(f"Missing: {cancer}.{pgs_id}*.anon_pgs")
            plot_idx += 1
            continue
        df = pd.read_csv(pgs_file[0], sep="\t")

        if "case" not in df.columns or "PGS" not in df.columns:
            plot_idx += 1
            continue

        ax = axes[plot_idx]

        # Colors
        case_fill = cancer_color[cancer.lower()]
        case_edge = darken_color(case_fill, 0.7)
        control_fill = "#B0B0B0"
        control_edge = "#4F4F4F"

        # Plot controls first, then cases on top
        sns.kdeplot(
            data=df[df['case'] == 0],
            x="PGS",
            ax=ax,
            fill=True,
            alpha=0.4,
            linewidth=1,
            color=control_fill,
            edgecolor=control_edge
        )
        sns.kdeplot(
            data=df[df['case'] == 1],
            x="PGS",
            ax=ax,
            fill=True,
            alpha=0.4,
            linewidth=1,
            color=case_fill,
            edgecolor=case_edge
        )

        # OR / p-value left annotation
        ax.text(
            PVAL_X, OR_PVAL_Y,
            f"P={p_val:.1e}\nOR={or_val:.1f}",
            transform=ax.transAxes,
            ha='left', va='top',
            fontsize=5
        )
        # OR / p-value left annotation
        # ax.text(
        #     OR_X, OR_PVAL_Y,
        #     f"OR={or_val:.1f}",
        #     transform=ax.transAxes,
        #     ha='left', va='top',
        #     fontsize=5,
        #     bbox=dict(facecolor="white", alpha=0.6, pad=0.05, linewidth=0)
        # )

        ax.set_title(f"{cancer.upper()} PRS",pad=0.5)
        ax.tick_params(width=0.8, length=2)
        ax.set_xlabel("")  # remove local x-label
        ax.set_ylabel("")  # remove local y-label
        plot_idx += 1

# -----------------------
# Hide unused axes
# -----------------------
for i in range(plot_idx, len(axes)):
    fig.delaxes(axes[i])

# -----------------------
# Axis labels (one for all)
# -----------------------
fig.text(0.5, 0, "PGS (Z-score)", ha='center', fontsize=6, fontweight='bold')
fig.text(0.02, 0.55, "Density", va='center', rotation='vertical', fontsize=6, fontweight='bold')

# -----------------------
# Legend
# -----------------------
legend_elements = [Patch(facecolor=control_fill, edgecolor=control_edge, label="Control")]
# fig.legend(
#     handles=legend_elements,
#     loc='lower left',
#     bbox_to_anchor=(1.01, 0.55),
#     frameon=False,
#     fontsize=5,
#     title="Group",
#     title_fontsize=5
# )

# -----------------------
# Layout + save
# -----------------------
plt.subplots_adjust(left=0.06, right=0.96, top=0.95, bottom=0.15, wspace=0.1, hspace=0.3)
plt.savefig("/Users/noah/Desktop/ufc_repository/results/analysis_4a_prs/Figure_4A.png", bbox_inches="tight",dpi=600,pad_inches=0)
plt.savefig("/Users/noah/Desktop/ufc_repository/results/analysis_4a_prs/Figure_4A.pdf", bbox_inches="tight",dpi=600, pad_inches=0)
plt.close()
