#!/usr/bin/env python3
import sys, os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import ttest_ind

# -------------------------------
# Styling constants
# -------------------------------
FIG_W, FIG_H = 1.7, 1.3
BASE_FONTSIZE = 6
TICK_FONTSIZE = 5
STAR_FONTSIZE = 5
LINEWIDTH = 0.4

plt.rcParams.update({
    "font.size": BASE_FONTSIZE,
    "axes.labelsize": BASE_FONTSIZE,
    "xtick.labelsize": TICK_FONTSIZE,
    "ytick.labelsize": TICK_FONTSIZE,
    "axes.linewidth": LINEWIDTH,
    "font.family": "Arial",
})

# -------------------------------
# Helper functions
# -------------------------------
def darken_color(color, factor=0.7):
    c = mcolors.to_rgb(color)
    return tuple(channel * factor for channel in c)

def parse_color_file(filepath):
    color_map = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            key, value = line.split(":", 1)
            color_map[key.strip().lower()] = value.strip().strip('"')
    return color_map

def paired_label(n, other_n):
    if n <= 20 or other_n <= 20:
        return "≤20" if n <= 20 else "≥20"
    return f"={n}"

def significance_text(p):
    if p > 0.05:
        return "ns"
    elif p > 1e-2:
        return "*"
    elif p > 1e-5:
        return "**"
    return "***"

def cohen_d(x1, x2):
    n1, n2 = len(x1), len(x2)
    pooled_sd = np.sqrt(
        ((n1 - 1) * x1.var() + (n2 - 1) * x2.var()) / (n1 + n2 - 2)
    )
    return (x1.mean() - x2.mean()) / pooled_sd

# -------------------------------
# Inputs / outputs
# -------------------------------
input_file = sys.argv[1]
base = os.path.splitext(os.path.basename(input_file))[0]
output_dir = "/Users/noah/Desktop/ufc_repository/results/analysis_4c_results"

png_file = f"{output_dir}/pngs/{base}.results.png"
pdf_file = f"{output_dir}/pngs/{base}.results.pdf"
stats_file = f"{output_dir}/stats/{base}.results.stats"

cancer_type, PGS_ID, *_ = base.split(".")
cancer_type_cap = cancer_type.capitalize()

# -------------------------------
# Load data
# -------------------------------
df = pd.read_csv(input_file, sep="\t")
groups = [
    "Control",
    f"Isolated {cancer_type_cap}",
    f"Familial {cancer_type_cap}",
]
df = df[df["group"].isin(groups)].copy()
counts = df["group"].value_counts().to_dict()

# -------------------------------
# Stats
# -------------------------------
pairs = [
    (groups[0], groups[1]),
    (groups[0], groups[2]),
    (groups[1], groups[2]),
]

p_values, or_values = {}, {}
for g1, g2 in pairs:
    x1 = df.loc[df.group == g1, "PGS"]
    x2 = df.loc[df.group == g2, "PGS"]
    stat, p2 = ttest_ind(x1, x2, equal_var=False)
    p1 = p2 / 2 if x1.mean() < x2.mean() else 1 - p2 / 2
    p_values[(g1, g2)] = p1
    d = cohen_d(x1, x2)
    or_values[(g1, g2)] = np.exp(d * np.pi / np.sqrt(3))

# -------------------------------
# Plot
# -------------------------------
colors = parse_color_file("/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml")
plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")

# Boxplot
sns.boxplot(
    data=df,
    x="group",
    y="PGS",
    order=groups,
    showfliers=False,
    linewidth=LINEWIDTH,
    palette={
        groups[0]: colors["control"],
        groups[1]: colors[cancer_type],
        groups[2]: darken_color(colors[cancer_type], 0.6),
    },
)

# Overlay points
n_iso = counts.get(groups[1], 0)
n_fam = counts.get(groups[2], 0)
for g in groups:
    if g == "Control" or (n_iso > 20 and n_fam > 20):
        sns.stripplot(
            data=df[df.group == g],
            x="group",
            y="PGS",
            order=groups,
            color="black",
            size=1,
            alpha=0.5,
        )

# -------------------------------
# Significance brackets (FIXED)
# -------------------------------
y_min, y_max = df.PGS.min(), df.PGS.max()
y_range = y_max - y_min

bracket_step = 0.12 * y_range     # vertical spacing between brackets
bracket_height = 0.02 * y_range  # height of bracket ticks
star_gap = -0.25       # distance from bracket to stars
star_gap_ns = 0

for i, (g1, g2) in enumerate(pairs):
    x1, x2 = groups.index(g1), groups.index(g2)
    y = y_max + (i + 1) * bracket_step

    # draw bracket
    plt.plot(
        [x1, x1, x2, x2],
        [y - bracket_height, y, y, y - bracket_height],
        lw=LINEWIDTH,
        c="black",
    )

    # get significance string
    sig_label = significance_text(p_values[(g1, g2)])

    # set gap
    if sig_label == "ns":
        gap = star_gap_ns
    else:
        gap = star_gap

    # draw stars
    plt.text(
        (x1 + x2) / 2,
        y + gap,
        sig_label,
        ha="center",
        va="bottom",
        fontsize=STAR_FONTSIZE,
    )



# -------------------------------
# Axes / labels
# -------------------------------
plt.ylabel(f"{cancer_type.replace('_'," ").title()} PRS", fontsize=BASE_FONTSIZE, labelpad=-1.3)

xtick_labels = [
    f"{groups[0]}\n(n={counts.get(groups[0],0)})",
    f"{groups[1].replace('Isolated ','Isolated\n').replace('_',' ').title()}\n(n{paired_label(n_iso, n_fam)})",
    f"{groups[2].replace('Familial ','Familial\n').replace('_',' ').title()}\n(n{paired_label(n_fam, n_iso)})",
]

ax = plt.gca()
plt.xticks(range(3), xtick_labels, rotation=20, ha="center")
ax.tick_params(axis='x', pad=-1)   # <-- bring labels closer to graph


# plt.xticks(range(3), xtick_labels, rotation=20)
# plt.tick_params(axis="y", length=2)

# Remove top/right spines
plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)

# Tight layout with minimal padding
plt.subplots_adjust(left=0.17, right=0.98, top=0.95, bottom=0.28)
plt.ylim(y_min - 0.05 * y_range, y_max + (len(pairs) + 1) * bracket_step)

# Save
plt.savefig(pdf_file, facecolor="white", edgecolor="white")
plt.savefig(png_file, facecolor="white", edgecolor="white", dpi=300)
plt.close()

# -------------------------------
# Write stats
# -------------------------------
summary = df.groupby("group")["PGS"].agg(["count", "mean", "std"])
with open(stats_file, "w") as f:
    f.write("# Summary statistics\n")
    f.write(summary.to_string() + "\n\n")
    f.write("# Pairwise comparisons\n")
    for (g1, g2), p in p_values.items():
        f.write(f"{g1} vs {g2}:\tp={p:.3e}\tOR≈{or_values[(g1,g2)]:.3f}\n")
