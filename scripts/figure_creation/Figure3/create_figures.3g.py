
#!/usr/bin/env python3

import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# -----------------------
# Paths
# -----------------------
BASE_DIR = Path(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all"
)

COLOR_FILE = Path(
    "/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml"
)

OUT = BASE_DIR / "cancer_or_comparison"

# -----------------------
# Cohorts
# Order: UFC, Replication, CDRv9
# -----------------------
cohorts = [
    "Discovery",
    "Replication",
    "Population",
]

# -----------------------
# Placeholder OR values
# Order: UFC, Replication, CDRv9
# Replace with actual ORs
# -----------------------
or_values = {
    "Prostate":      [2.18, 1.50, 1.25],
    "Breast":        [1.82, 1.55, 1.20],
    "Melanoma":      [1.77, 1.45, 1.12],
    "Thyroid":       [1.61, 1.15, 1.05],
    "Leukemia":      [1.50, 1.25, 1.10],
    "Colorectal":    [1.47, 1.15, 1.08],
    "Non-Hodgkin":   [1.43, 1.20, 1.08],
    "Kidney":        [1.39, 1.20, 1.06],
    "Endometrial":   [1.32, 1.10, 1.03],
    "Bladder":       [1.23, 1.10, 1.05],
    "Lung":          [1.17, 1.05, 1.02],
    "Ovarian":       [1.07, 1.30, 1.10],
}

cancer_types = list(or_values.keys())

# -----------------------
# Load color scheme
# -----------------------
with open(COLOR_FILE, "r") as f:
    color_scheme = yaml.safe_load(f)

# Assumes these are top-level YAML keys
cohort_colors = {
    "Discovery": color_scheme["Discovery"],
    "Replication": color_scheme["Replication"],
    "Population": color_scheme["Population"],
}

# -----------------------
# Plot styling
# -----------------------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
})

fig, ax = plt.subplots(figsize=(2, 3))

y = np.arange(len(cancer_types))

# Three bars per cancer type
bar_height = 0.23
offsets = [-bar_height, 0, bar_height]

# -----------------------
# Grouped horizontal bar chart
# -----------------------
for i, cohort in enumerate(cohorts):
    ors = [or_values[cancer][i] for cancer in cancer_types]

    ax.barh(
        y + offsets[i],
        ors,
        height=bar_height * 0.9,
        color=cohort_colors[cohort],
        label=cohort,
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
    )

# Reference line: OR = 1
ax.axvline(
    1,
    color="black",
    linestyle="--",
    linewidth=0.9,
    alpha=0.7,
    zorder=2,
)

# -----------------------
# Axes formatting
# -----------------------
ax.set_yticks(y)
ax.set_yticklabels(cancer_types)
ax.invert_yaxis()

ax.set_xlabel("Odds Ratio (OR) per SD")
# ax.set_ylabel("Cancer Type")

# Start x-axis at zero
max_or = max(
    max(values) for values in or_values.values()
)
ax.set_xlim(0.9, max_or)

# -----------------------
# Legend
# -----------------------
# ax.legend(
#     frameon=False,
#     loc="lower right",
#     ncol=1,
# )
leg = ax.legend(
    frameon=False,
    loc="lower right",
    bbox_to_anchor=(1.1, 0.0),
    bbox_transform=ax.transAxes,
    ncol=1,
)
# -----------------------
# Grid and spines
# -----------------------
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
#ax.spines["left"].set_visible(False)

ax.grid(
    axis="x",
    linestyle=":",
    linewidth=0.6,
    alpha=0.4,
)

ax.set_axisbelow(True)

# -----------------------
# Layout
# -----------------------
plt.subplots_adjust(
    left=0.1,   # Reduce left margin; adjust for cancer labels
    right=0.78,  # Reserve space for the legend inside the figure
    top=0.98,
    bottom=0.16,
)

# -----------------------
# Save
# -----------------------
plt.savefig(f"/Users/noah/Desktop/ufc_repository/results/Figure4/4g_comparison.png", dpi=300, bbox_inches="tight")
plt.savefig(f"/Users/noah/Desktop/ufc_repository/results/Figure4/4g_comparison.pdf", bbox_inches="tight")

#plt.show()
