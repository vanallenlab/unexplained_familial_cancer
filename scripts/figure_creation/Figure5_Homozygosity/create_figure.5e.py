#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Patch

# ------------------------------------------------------------
# Hardcoded inputs / outputs
# ------------------------------------------------------------
IN_TSV = "/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/thyroid_roh_results/thyroid.chr3_113500001_113950000.roh_signal.tsv.gz"
OUT_PDF = "/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/thyroid_roh_results/thyroid.chr3_113500001_113950000.roh_signal.plot.pdf"
OUT_PNG = "/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/thyroid_roh_results/thyroid.chr3_113500001_113950000.roh_signal.plot.png"

# Locus (hardcoded based on your haplotype string)
CHR = "chr3"
START = 113500001
END   = 113950000

# Plot settings
DOWNSAMPLE_BP = 10          # plot every 100 bp (fast + looks identical)
FIGSIZE = (3,2)         # inches (wide and short = good for genome tracks)
DPI = 450

# Colors
CASE_COLOR = "#9EDAE5"       # nice green
CTRL_COLOR = "#B0B0B0"       # clean gray

FONT_FAMILY = "Arial"

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------
df = pd.read_csv(IN_TSV, sep="\t", compression="gzip")

required = {"chr", "pos", "case_roh_fraction", "control_roh_fraction"}
missing = required - set(df.columns)
if missing:
    raise SystemExit(f"Missing required columns: {sorted(missing)}")

# Force types (prevents weird string bugs)
df["pos"] = df["pos"].astype(int)
df["case_roh_fraction"] = df["case_roh_fraction"].astype(float)
df["control_roh_fraction"] = df["control_roh_fraction"].astype(float)

# Restrict to requested locus (just in case the file contains extra positions)
df = df[(df["chr"] == CHR) & (df["pos"] >= START) & (df["pos"] <= END)].copy()
df = df.sort_values("pos")

if df.empty:
    raise SystemExit("No rows found in the requested locus range — check CHR/START/END.")

# ------------------------------------------------------------
# Downsample
#   Adjacent bases have identical values most of the time.
#   Plotting every 100 bp is visually the same but much faster.
# ------------------------------------------------------------
df = df.iloc[::DOWNSAMPLE_BP, :].copy()

# ------------------------------------------------------------
# Nat Gen-ish style
# ------------------------------------------------------------
plt.rcParams.update({
    "font.family": FONT_FAMILY,
    "font.size": 7,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=FIGSIZE)

ax.axvspan(113600001,113850000, color='lightblue', alpha=0.4, zorder=0)

x_start = 113600001
x_end   = 113850000
y = 0.025

# double-headed arrow
ax.annotate(
    "",
    xy=(x_start, y),
    xytext=(x_end, y),
    arrowprops=dict(
        arrowstyle="<->",
        linewidth=1.5,
        color="black"
    ),
    zorder=3
)

# centered label above the arrow
ax.text(
    (x_start + x_end) / 2,
    y * 1.05,   # small vertical offset
    "250 kb",
    ha="center",
    va="bottom",
    fontsize=7
)
# Lines (simple + clean)
line_case, = ax.plot(
    df["pos"],
    df["case_roh_fraction"],
    linewidth=1.2,
    color=CASE_COLOR,
    label="Thyroid\nCases"
)

line_control, = ax.plot(
    df["pos"],
    df["control_roh_fraction"],
    linewidth=1.2,
    color=CTRL_COLOR,
    label="Controls"
)

# Axis formatting
ax.set_xlim(START, END)
ax.set_ylabel("% homozygosed (RoH ≥100 kb)")


# X ticks: keep them readable
xticks = [START,(START + END) // 2, END]
ax.set_xticks(xticks)
ax.set_xticklabels([f"{x:,}" for x in xticks])
ax.set_xlabel(f"chr3")

# Clean look
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Subtle horizontal grid only
ax.yaxis.grid(True, linewidth=0.4, alpha=0.25)
ax.set_axisbelow(True)

# Legend (no box)
# Create a light blue patch for the legend
sig_patch = Patch(facecolor='lightblue', alpha=0.4, label='Significant\nregion')

#ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.8, 1.00))
ax.legend(handles=[line_case, line_control, sig_patch], frameon=False, loc="upper left", bbox_to_anchor=(0.8, 1.0))

# Tight margins (kills the annoying whitespace)
fig.tight_layout(pad=0.6)

# Plot genes
# SPICE1
plt.plot([113442718,113515187], [0.2, 0.2], color='black', lw=3,solid_capstyle='round')

# SIDT1
plt.plot([113532555,113629575], [0.2, 0.2], color='black', lw=3,solid_capstyle='round')

# ATP6V1A
plt.plot([113747033,113812056], [0.25, 0.25], color='black', lw=3,solid_capstyle='round')

# USF3 
plt.plot([113648385,113696646], [0.2, 0.2], color='red', lw=3,solid_capstyle='round')

# GRAMD1C
plt.plot([113828182,113947174], [0.2, 0.2], color='black', lw=3,solid_capstyle='round')

# ZDHHC23
plt.plot([113947901,113965401], [0.23, 0.25], color='black', lw=3,solid_capstyle='round')

#NAA50 3:113716458-113746300
plt.plot([113716458,113746300], [0.15, 0.15], color='black', lw=3,solid_capstyle='round')

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
Path(OUT_PDF).parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PDF, dpi=DPI, bbox_inches="tight", facecolor="white",pad_inches=0)
fig.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("Saved:", OUT_PDF)
print("Saved:", OUT_PNG)
print(f"Plotted {len(df):,} points (every {DOWNSAMPLE_BP} bp)")