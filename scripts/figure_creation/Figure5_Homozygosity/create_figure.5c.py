#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Patch

# ------------------------------------------------------------
# Hardcoded inputs / outputs
# ------------------------------------------------------------
IN_TSV = "/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/kidney_roh_results/kidney.chr22_32100001_32550000.roh_signal.tsv.gz"
IN_TSV = "/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/kidney_roh_results/kidney.chr22_32000001_32650000.roh_signal.tsv.gz"
OUT_PDF = "/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/kidney_roh_results/kidney.chr22_32100001_32550000.roh_signal.plot.pdf"
OUT_PNG = "/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/kidney_roh_results/kidney.chr22_32100001_32550000.roh_signal.plot.png"

# Locus (hardcoded based on your haplotype string)
CHR = "chr22"
START = 32000001
END   = 32660000

# Plot settings
DOWNSAMPLE_BP = 10          # plot every 100 bp (fast + looks identical)
FIGSIZE = (4.5,3.5)         # inches (wide and short = good for genome tracks)
DPI = 450

# Colors
CASE_COLOR = "#98DF8A"       # nice green
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

ax.axvspan(32200001, 32450000, color='lightblue', alpha=0.4, zorder=0)

x_start = 32200001
x_end   = 32450000
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

# # Add a right axis
# ax_r = ax.twinx()

# # keep it visually aligned with the left axis
# #ax_r.set_ylim(ax.get_ylim())

# # right-axis labels only
# ax_r.set_yticks([0.1, 0.2, 0.3])
# #ax_r.set_yticklabels(["0", "2", "4"])

# ax_r.set_ylabel(r"$-\log_{10}(P)$")
### finish adding right axis


# Lines (simple + clean)
line_case, = ax.plot(
    df["pos"],
    df["case_roh_fraction"],
    linewidth=1.2,
    color=CASE_COLOR,
    label="Kidney\nCases"
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
xticks = [START, START + (32650000 - START) // 4,(START + 32650000) // 2, START + 3 * (32650000 - START) // 4, 32650000]
ax.set_xticks(xticks)
ax.set_xticklabels([f"{x:,}" for x in xticks])
ax.set_xlabel(f"chr22")

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
import numpy as np
# BPIFC
pval = 0.1 + (-np.log10(0.000488006777292993))/4 * 0.2
plt.plot([32413845, 32464484], [pval, pval], color='blue', lw=3,solid_capstyle='round')
ax.text(32412582, pval - 0.02, "BPIFC", color="black",fontstyle="italic", fontsize=7)

#RTCB
pval = 0.1 + (-np.log10(0.675322635873349))/4 * 0.2
plt.plot([32387582,32412248], [pval, pval], color='black', lw=3,solid_capstyle='round')
ax.text(32377582, pval - 0.02, "RTCB", color="black",fontstyle="italic", fontsize=7)

#FBOX7 32474676-32498829
plt.plot([32474676,32498829], [0.075, 0.075], color='darkgray', lw=3,solid_capstyle='round')
ax.text(32459676, 0.075 - 0.02, "FBOX7", color="black",fontstyle="italic", fontsize=7)

#SYN3 32507820-33058381
# pval = 0.1 + (-np.log10(0.675322635873349))/4 * 0.2
# plt.plot([332507820,33058381], [pval, pval], color='black', lw=3,solid_capstyle='round')

#RFPL3 22:32354885-32361161
pval = 0.1 + (-np.log10(0.790594840301569))/4 * 0.2
plt.plot([32354885,32361161], [pval, pval], color='black', lw=3,solid_capstyle='round')
ax.text(32321885, pval - 0.02, "RFPL3", color="black",fontstyle="italic", fontsize=7)

#SLC5A4 22:32218464-32255347
pval = 0.1 + (-np.log10(0.299226453279006))/4 * 0.2
plt.plot([32218464,32255347], [pval, pval], color='black', lw=3,solid_capstyle='round')
ax.text(32208464, pval - 0.02, "SLC5A4", color="black",fontstyle="italic", fontsize=7)

#RFPL2 22:32190435-32205073
pval = 0.1 + (-np.log10(0.221030545692232))/4 * 0.2
plt.plot([32190435,32205073], [pval, pval], color='black', lw=3,solid_capstyle='round')
ax.text(32155435, pval - 0.02, "RFPL2", color="black",fontstyle="italic", fontsize=7)

# C22orf42 22:32149006-32159322
pval = 0.1 + (-np.log10(0.839865239796788))/4 * 0.2
plt.plot([32149006,32159322], [pval, pval], color='black', lw=3,solid_capstyle='round')
ax.text(32120000, pval - 0.02, "C22orf42", color="black",fontstyle="italic", fontsize=7)

# SLC5A1 22:32043261-32113029
pval = 0.1 + (-np.log10(0.850403772168157))/4 * 0.2
plt.plot([32043261,32113029], [pval, pval], color='black', lw=3,solid_capstyle='round')
ax.text(32043261, pval - 0.02, "SLC5A1", color="black",fontstyle="italic", fontsize=7)

# SYN3 32507820-33058381
pval = 0.1 + (-np.log10(0.811517764065727))/4 * 0.2
plt.plot([32507820,32547820], [pval, pval], color='black', lw=3,solid_capstyle='round')
ax.text(32537820, pval - 0.02, "SYN3", color="black",fontstyle="italic", fontsize=7)

ax.text(32010001, 0.42, "P=1.02e-6", color="black", fontsize=7)
ax.text(32010001, 0.4, "OR=9.63", color="black", fontsize=7)

ax.annotate(
    "",
    xy=(32635000, pval),     # arrow tip
    xytext=(32547820, pval),
    arrowprops=dict(
        arrowstyle="-|>",
        linewidth=3,
        color="black",
        shrinkA=0,
        shrinkB=0
    ),
    zorder=3
)

# Add a right axis
plt.plot([32640000,32640000], [0.1, 0.3], color='black', lw=1,solid_capstyle='round')
#ticks
plt.plot([32640000,32645000], [0.1, 0.1], color='black', lw=1,solid_capstyle='round')
ax.text(32645000 + 2000, 0.1, "0", va="center", ha="left", fontsize=7)

plt.plot([32640000,32645000], [0.2, 0.2], color='black', lw=1,solid_capstyle='round')
ax.text(32645000 + 2000, 0.2, "2", va="center", ha="left", fontsize=7)

plt.plot([32640000,32645000], [0.3, 0.3], color='black', lw=1,solid_capstyle='round')
ax.text(32645000 + 2000, 0.3, "4", va="center", ha="left", fontsize=7)
ax.text(32660000, 0.2, r"$-\log_{10} (P)$ Kidney Cancer in the UK Biobank",
        rotation=90, va="center", ha="left", fontsize=7)
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
