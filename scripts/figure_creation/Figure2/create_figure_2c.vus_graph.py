#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------
# Load data
# -----------------------
df = pd.read_csv("/Users/noah/Desktop/ufc_repository/results/analysis_5_gsea_results/tier_variant_case_control_summary.tsv", sep="\t")

# Ensure correct dtypes
df["Tier"] = df["Tier"].astype(str)
df["AF"] = df["AF"].astype(float)

# Only tiers that actually exist
TIERS = ["Tier0", "Tier1", "Tier2", "Tier3", "Tier4"]
df = df[df["Tier"].isin(TIERS)]

# -----------------------
# Split by AF and reindex
# -----------------------
df_001 = (
    df[df["AF"] == 0.001]
    .set_index("Tier")
    .reindex(TIERS)
    .fillna(0)
)

df_01 = (
    df[df["AF"] == 0.01]
    .set_index("Tier")
    .reindex(TIERS)
    .fillna(0)
)

# -----------------------
# Compute percentages
# -----------------------
cases_001 = df_001["Cases_with_variant"] / df_001["Total_cases"] * 100
cases_01  = df_01["Cases_with_variant"]  / df_01["Total_cases"]  * 100

ctrls_001 = df_001["Controls_with_variant"] / df_001["Total_controls"] * 100
ctrls_01  = df_01["Controls_with_variant"]  / df_01["Total_controls"]  * 100

# Increment from AF ≤ 0.001 → AF ≤ 0.01
cases_inc = cases_01 - cases_001
ctrls_inc = ctrls_01 - ctrls_001

# -----------------------
# Plot
# -----------------------
fig, ax = plt.subplots(figsize=(10, 6))

bar_width = 0.25
x = np.arange(len(TIERS))

x_cases = x - bar_width / 2
x_ctrls = x + bar_width / 2

# --- Cases ---
ax.bar(
    x_cases,
    cases_001,
    width=bar_width,
    color="firebrick",
    label="Cases (AF ≤ 0.001)"
)
ax.bar(
    x_cases,
    cases_inc,
    bottom=cases_001,
    width=bar_width,
    color="salmon",
    alpha=0.7,
    label="Cases (0.001 < AF ≤ 0.01)"
)

# --- Controls ---
ax.bar(
    x_ctrls,
    ctrls_001,
    width=bar_width,
    color="gray",
    label="Controls (AF ≤ 0.001)"
)
ax.bar(
    x_ctrls,
    ctrls_inc,
    bottom=ctrls_001,
    width=bar_width,
    color="darkgray",
    alpha=0.7,
    label="Controls (0.001 < AF ≤ 0.01)"
)

# -----------------------
# Annotate bar tops
# -----------------------
def annotate_bar_tops(ax, x_positions, heights, side="right"):
    for x_pos, h in zip(x_positions, heights):
        if h <= 0:
            continue

        if side == "right":
            ax.plot(
                [x_pos + bar_width / 2, x_pos + bar_width / 2 + 0.05],
                [h, h],
                lw=1,
                color="black"
            )
            text_x = x_pos + bar_width / 2 + 0.06
            ha = "left"
        else:
            ax.plot(
                [x_pos - bar_width / 2, x_pos - bar_width / 2 - 0.05],
                [h, h],
                lw=1,
                color="black"
            )
            text_x = x_pos - bar_width / 2 - 0.06
            ha = "right"

        ax.text(
            text_x,
            h,
            f"{int(round(h))}%",
            va="center",
            ha=ha,
            fontsize=11
        )

# Final totals (AF ≤ 0.01)
annotate_bar_tops(ax, x_cases, cases_01, side="left")
annotate_bar_tops(ax, x_cases, cases_001, side="left")
annotate_bar_tops(ax, x_ctrls, ctrls_01, side="right")
annotate_bar_tops(ax, x_ctrls, ctrls_001, side="right")

# -----------------------
# Formatting
# -----------------------
ax.set_xticks(x)
ax.set_xticklabels(TIERS,fontsize=14)

#ax.set_xlabel("Variant Tier")
#ax.set_ylabel("Percent of individuals with ≥1 variant", fontsize=14)

ax.set_ylim(0, 100)
ax.set_yticks(np.arange(0, 101, 20))
ax.set_yticklabels([f"{y}%" for y in range(0, 101, 20)],fontsize=14)

#ax.set_title("Percentage of Individuals Carrying VUS's in 91 CPGs")

ax.legend(frameon=False, ncol=2)

##### Remove top and right spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Keep left and bottom spines
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

# Optional: make axes slightly thicker
ax.spines["left"].set_linewidth(1.2)
ax.spines["bottom"].set_linewidth(1.2)

# Ensure ticks are only on left/bottom
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")
#####

plt.tight_layout()
plt.savefig("tier_side_by_side_plot.png", dpi=300)
plt.show()
