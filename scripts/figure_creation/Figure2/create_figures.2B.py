import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.lines import Line2D
from matplotlib import font_manager as fm

# -----------------------------
# Load data
# -----------------------------
file_path = Path("/Users/noah/Desktop/ufc_repository/results/Figure_2/sv_lof_fisher_summary.tsv")
df = pd.read_csv(file_path, sep="\t")

# -----------------------------
# Filter cohorts
# -----------------------------
df = df[
    ~df["cohort"].str.contains("isolated", case=False, na=False) &
    ~df["cohort"].str.contains("patient_and_family", case=False, na=False)
].copy()

# -----------------------------
# Clean display names (for axis only)
# -----------------------------
cohort_display = []
for c in df["cohort"].unique():
    name = c.lower()
    if "non-hodgkin" in name:
        cohort_display.append("NHL")
    elif "neuroendocrine" in name:
        cohort_display.append("NETs")
    elif "basal_cell" in name:
        cohort_display.append("BCC")
    elif "squamous_cell" in name:
        cohort_display.append("SCC")
    else:
        name = name.replace("_", " ")
        cohort_display.append(name.title())

cohort_mapping = dict(zip(df["cohort"].unique(), cohort_display))
df["cohort_display"] = df["cohort"].map(cohort_mapping)

# -----------------------------
# Order cohorts by rare OR
# -----------------------------
order = (
    df[df["sv_type"] == "singleton"]
    .sort_values("odds_ratio", ascending=False)["cohort_display"]
)
df["cohort_display"] = pd.Categorical(df["cohort_display"], categories=order, ordered=True)
df = df.sort_values(["cohort_display", "sv_type"])

# -----------------------------
# Colors and offsets
# -----------------------------
colors = {
    "rare": "#1f77b4",         # blue
    "ultra_rare": "#ff7f0e",   # orange
    "singleton": "#2ca02c"     # green
}
offsets = {"rare": -0.25, "ultra_rare": 0, "singleton": 0.25}
offsets = {"singleton": -0.25, "ultra_rare": 0, "rare": 0.25}
# -----------------------------
# Plot
# -----------------------------
fig, ax = plt.subplots(figsize=(3.5, 2))

cohorts = df["cohort_display"].cat.categories
x_base = np.arange(len(cohorts))
bar_width = 0.2  # width of bars

for sv_type in ["rare", "ultra_rare", "singleton"]:
    sub = df[df["sv_type"] == sv_type]
    x_positions = [
        x_base[list(cohorts).index(c)] + offsets[sv_type]
        for c in sub["cohort_display"]
    ]
    bars = ax.bar(
        x_positions,
        sub["odds_ratio"],
        width=bar_width,
        color=colors[sv_type],
        edgecolor="black",
        linewidth=0.4
    )

    # Add * for p < 0.05
    for x, or_val, pval in zip(x_positions, sub["odds_ratio"], sub["p_value"]):
        if pval < 0.05:
            ax.text(x, or_val + 0.3, "*", ha="center", va="bottom", fontsize=6)

# Null line at OR=1
ax.axhline(1, linestyle="--", linewidth=0.6, color="black")

# -----------------------------
# Axis and labels
# -----------------------------
ax.set_xticks(x_base)
ax.set_xticklabels(cohorts, rotation=45, ha="right", fontsize=5)
ax.set_ylabel("Odds Ratio", fontsize=5)
ax.set_xlabel("")
ax.tick_params(axis='both', labelsize=5)
ax.margins(x=0)
ax.set_xlim(-0.6, len(cohorts) - 0.4)
ax.set_ylim(0, 10)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# -----------------------------
# Legend (bottom left)
# -----------------------------
legend_elements = [
    Line2D([0], [0], marker='s', color=colors["singleton"], label='Singleton',
           markersize=3, linestyle='None'),
    Line2D([0], [0], marker='s', color=colors["ultra_rare"], label='AF < 0.1%',
           markersize=3, linestyle='None'),
    Line2D([0], [0], marker='s', color=colors["rare"], label='AF < 1%',
           markersize=3, linestyle='None')
]

ax.legend(
    handles=legend_elements,
    frameon=False,
    loc="upper right",
    prop=fm.FontProperties(family="Arial", size=5)
)

ax.text(15.5, 7, '* p ≤ 0.05', fontsize=5, ha='center', va='bottom', fontfamily="Arial")
ax.text(8, 9, 'Odds Ratios of LoF SVs by Cancer', fontsize=5, ha='center', va='bottom', fontfamily="Arial")
plt.tight_layout(pad=0)

# -----------------------------
# Save
# -----------------------------
output_path = file_path.parent / "Figure2B.pdf"
plt.savefig(output_path)
plt.close()

print(f"Saved clustered bar forest plot with significance stars to {output_path}")