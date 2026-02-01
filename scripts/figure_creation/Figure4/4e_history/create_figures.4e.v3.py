#!/usr/bin/env python3
import sys, os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import ttest_ind

# -------------------------------
# Styling constants (MATCHED)
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
    return tuple(ch * factor for ch in c)

def parse_color_file(filepath):
    cmap = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            k, v = line.split(":", 1)
            cmap[k.strip().lower()] = v.strip().strip('"')
    return cmap

def paired_label(n, other):
    if n <= 20:
        return "≤20"
    if other <= 20:
        return "≥20"
    return f"={n}"

def significance_text(p):
    if p > 0.05:
        return "ns"
    elif p > 1e-2:
        return "*"
    elif p > 1e-5:
        return "**"
    elif p > 1e-10:
        return "***"
    return "****"

def cohen_d(x1, x2):
    n1, n2 = len(x1), len(x2)
    pooled = np.sqrt(
        ((n1 - 1) * x1.var() + (n2 - 1) * x2.var()) / (n1 + n2 - 2)
    )
    return (x1.mean() - x2.mean()) / pooled

# -------------------------------
# Inputs / outputs
# -------------------------------
input_file = sys.argv[1]
base = os.path.splitext(os.path.basename(input_file))[0]

out_root = "/Users/noah/Desktop/ufc_repository/results/analysis_4e_results"
png_file = f"{out_root}/pngs/{base}.results.png"
pdf_file = f"{out_root}/pngs/{base}.results.pdf"
stats_file = f"{out_root}/stats/{base}.results.stats"

cancer_type, PGS_ID, *_ = base.split(".")
cancer_cap = cancer_type.replace("_", " ").title()

colors = parse_color_file(
    "/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml"
)

# -------------------------------
# Load data
# -------------------------------
df = pd.read_csv(input_file, sep="\t")
df["group"] = df["group"].str.replace("control", "Control")

groups = [
    "Control",
    f"Not-Inherited {cancer_type.capitalize()}",
    f"Familial {cancer_type.capitalize()}",
]

df = df[df["group"].isin(groups)].copy()
counts = df.group.value_counts().to_dict()

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

    _, p2 = ttest_ind(x1, x2, equal_var=False)
    p1 = p2 / 2 if x1.mean() < x2.mean() else 1 - p2 / 2
    p_values[(g1, g2)] = p1

    d = cohen_d(x1, x2)
    or_values[(g1, g2)] = np.exp(d * np.pi / np.sqrt(3))

# -------------------------------
# Plot
# -------------------------------
plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")

ax = sns.boxplot(
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

# -------------------------------
# Hatch discordant box
# -------------------------------
disc_idx = groups.index(groups[1])
box = ax.patches[disc_idx]
box.set_facecolor(colors["control"])
box.set_edgecolor(colors[cancer_type])
box.set_hatch("///")
box.set_linewidth(LINEWIDTH)

# -------------------------------
# Overlay points (paired ≥20 logic)
# -------------------------------
n_disc = counts.get(groups[1], 0)
n_fam = counts.get(groups[2], 0)
paired_ok = n_disc > 20 and n_fam > 20

for g in groups:
    if g == "Control" or paired_ok:
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
# Significance brackets (STARS)
# -------------------------------
y_min, y_max = df.PGS.min(), df.PGS.max()
y_range = y_max - y_min

bracket_step = 0.12 * y_range
bracket_height = 0.02 * y_range
star_gap = -0.25       # distance from bracket to stars
star_gap_ns = 0

for i, (g1, g2) in enumerate(pairs):
    x1, x2 = groups.index(g1), groups.index(g2)
    y = y_max + (i + 1) * bracket_step

    ax.plot(
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
        
    
    ax.text(
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
ax.set_ylabel(f"{cancer_type.title()} PRS", labelpad=-1.3)
ax.set_xlabel("")

xticks = [
    f"Control\n(n={counts.get(groups[0],0)})",
    f"{cancer_cap} Family\nDiscordant Case\n(n{paired_label(n_disc, n_fam)})",
    f"Familial \n {cancer_cap}\n(n{paired_label(n_fam, n_disc)})",
]

ax.set_xticks(range(3))
ax.set_xticklabels(xticks, rotation=20)
ax.tick_params(axis="x", pad=-1)

# Remove top/right spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Layout
plt.subplots_adjust(left=0.18, right=0.98, top=0.95, bottom=0.28)
ax.set_ylim(
    y_min - 0.05 * y_range,
    y_max + (len(pairs) + 1) * bracket_step,
)

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
