#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# -----------------------------
# CONFIG
# -----------------------------
TSV = "/Users/noah/Desktop/ufc_repository/results/Figure_2/sv_carrier_rates.case_vs_control.riaz_genes.pvalues.tsv"
OUT_PDF = "/Users/noah/Desktop/ufc_repository/results/Figure_2/sv_carrier_rates.case_vs_control.riaz_genes.NG.pdf"
OUT_PNG = "/Users/noah/Desktop/ufc_repository/results/Figure_2/sv_carrier_rates.case_vs_control.riaz_genes.NG.png"

TITLE = "SV carrier rate in riaz genes (QC-pass only)"
YLAB = "% SV carriers"
SHOW_TITLE = False  # Nature Genetics usually prefers figure legends over big titles

# publication defaults
DPI = 450
FIG_W, FIG_H = 3,2   # inches; adjust height if many cancers
FONT_FAMILY = "Arial"

# -----------------------------
# Load
# -----------------------------
df = pd.read_csv(TSV, sep="\t")

required = ["cancer", "pct_case", "pct_ctrl"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit(f"Missing required columns in TSV: {missing}")

# Order cancers by pct_case
df = df.sort_values("pct_case", ascending=False).reset_index(drop=True)

# -----------------------------
# Styling (Nature Genetics-ish)
# -----------------------------
plt.rcParams.update({
    "font.family": FONT_FAMILY,
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,  # embed fonts properly
    "ps.fonttype": 42
})

# -----------------------------
# Plot
# -----------------------------
n = df.shape[0]
x = np.arange(n)
bar_w = 0.38

# dynamically adjust width/height if many cancers
# fig_h = max(FIG_H, 0.18 * n) if n > 12 else FIG_H
# fig_w = max(FIG_W, 0.35 * n) if n > 18 else FIG_W

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

# Colors: clean muted palette
case_color = "#3B6FB6"   # muted blue
ctrl_color = "#B0B0B0"   # light gray

bars_case = ax.bar(
    x - bar_w/2,
    df["pct_case"].values,
    width=bar_w,
    label="Cases",
    edgecolor="black",
    linewidth=0.4,
    color=case_color
)

bars_ctrl = ax.bar(
    x + bar_w/2,
    df["pct_ctrl"].values,
    width=bar_w,
    label="Matched controls",
    edgecolor="black",
    linewidth=0.4,
    color=ctrl_color
)

# Axes formatting
ax.set_ylabel(YLAB)
ax.set_xticks(x)
ax.set_xticklabels(df["cancer"].values, rotation=45, ha="right")

# Remove top/right spines (clean journal look)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# subtle y grid
ax.yaxis.grid(True, linestyle="-", linewidth=0.4, alpha=0.25)
ax.set_axisbelow(True)

# Y limit with headroom for "ns"
ymax = float(max(df["pct_case"].max(), df["pct_ctrl"].max()))
if ymax < 1:
    ax.set_ylim(0, 1.2)
else:
    ax.set_ylim(0, ymax * 1.28)

# Add "ns" above each cancer (centered between bars)
# Put it slightly above the max of the two bars for that cancer
for i, row in df.iterrows():
    y = max(row["pct_case"], row["pct_ctrl"])
    ax.text(
        i,
        y + (ax.get_ylim()[1] * 0.02),
        "ns",
        ha="center",
        va="bottom",
        fontsize=7
    )

# Legend
ax.legend(frameon=False, loc="upper right")

# Optional title (usually off)
if SHOW_TITLE:
    ax.set_title(TITLE)

# Tight layout (minimize wasted margins)
fig.tight_layout(pad=0.6)

# Save
Path(OUT_PDF).parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PDF, dpi=DPI, bbox_inches="tight", facecolor="white",pad_inches=0)
fig.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight", facecolor="white",pad_inches=0)
plt.close(fig)

print("Saved:", OUT_PDF)
print("Saved:", OUT_PNG)
