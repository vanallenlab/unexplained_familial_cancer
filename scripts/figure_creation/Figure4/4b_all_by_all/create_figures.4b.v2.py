#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# -----------------------
# Global style
# -----------------------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 5,
    "axes.labelsize": 7,
    "xtick.labelsize": 5,
    "ytick.labelsize": 5
})

# -----------------------
# Load data
# -----------------------
pval_df = pd.read_csv(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.pvalue.tsv",
    sep="\t"
)

or_df = pd.read_csv(
    "/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.OR.tsv",
    sep="\t"
)

# -----------------------
# Filter PGS
# -----------------------
pgs_keep = [
    "PGS000785","PGS000789","PGS000791","PGS000784",
    "PGS000787","PGS003382","PGS004694","PGS000797","PGS004687",
    "PGS003384","PGS000788","PGS004249","PGS004244","PGS004688"
]

def filter_pgs(df):
    df = df[~df["PGS_ID"].str.contains("-euro", case=False, na=False)]
    return df[df["PGS_ID"].isin(pgs_keep)]

pval_df = filter_pgs(pval_df)
or_df = filter_pgs(or_df)

# -----------------------
# Melt to long format
# -----------------------
pval_long = pval_df.melt(
    id_vars=["Intended_Cancer_Type", "PGS_ID"],
    var_name="Cancer",
    value_name="p_value"
)

or_long = or_df.melt(
    id_vars=["Intended_Cancer_Type", "PGS_ID"],
    var_name="Cancer",
    value_name="OR_raw"
)

# -----------------------
# Parse OR
# -----------------------
def parse_or(x):
    m = re.match(r"([\d\.Ee+-]+)", str(x))
    return float(m.group(1)) if m else np.nan

or_long["OR"] = or_long["OR_raw"].apply(parse_or)

# -----------------------
# Merge
# -----------------------
df = pd.merge(
    pval_long,
    or_long[["Intended_Cancer_Type", "PGS_ID", "Cancer", "OR"]],
    on=["Intended_Cancer_Type", "PGS_ID", "Cancer"]
)

df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
df = df.dropna(subset=["p_value", "OR"])

# -----------------------
# Clean labels
# -----------------------
df["Cancer"] = df["Cancer"].str.replace("_", " ").str.title()
df["Intended_Cancer_Type"] = df["Intended_Cancer_Type"].str.replace("_", " ").str.title()

# Compact but informative x labels
df["x_label"] = df["Intended_Cancer_Type"] + " PRS"

# -----------------------
# Transform values
# -----------------------
df["neglog10p"] = (-np.log10(df["p_value"])).clip(0, 8)
df["OR"] = df["OR"].clip(0.5, 4)

# -----------------------
# Axis ordering
# -----------------------
x_order = sorted(df["x_label"].unique())
y_order = sorted(df["Cancer"].unique())

df["x_label"] = pd.Categorical(df["x_label"], categories=x_order, ordered=True)
df["Cancer"] = pd.Categorical(df["Cancer"], categories=y_order, ordered=True)

df["x_pos"] = df["x_label"].cat.codes
df["y_pos"] = df["Cancer"].cat.codes

# -----------------------
# Significance
# -----------------------
sig = (df["p_value"] <= 0.05) & (df["OR"] > 1)

# -----------------------
# Plot
# -----------------------
fig, ax = plt.subplots(figsize=(4, 2))

sc = ax.scatter(
    df["x_pos"],
    df["y_pos"],
    s=df["OR"] ** 3 * 4,
    c=df["neglog10p"],
    cmap="Reds",
    alpha=0.9,
    linewidth=0
)

# Significant overlay
ax.scatter(
    df.loc[sig, "x_pos"],
    df.loc[sig, "y_pos"],
    s=df.loc[sig, "OR"] ** 3 * 5,
    facecolors="none",
    edgecolors="blue",
    linewidth=0.5
)

#### Added OR scaling value ####
or_vals = [0.5, 1, 2]

size_handles = [
    plt.scatter(
        [], [], 
        s=or_val ** 3 * 4,
        color="white",
        edgecolor="black",
        linewidth=0.6,
        label=f"OR = {or_val}"
    )
    for or_val in or_vals
]

# Create a legend handle for significant points
sig_handle = plt.scatter(
    [], [], 
    s=5,                # choose a size similar to your points
    c="white",            # fill color (matches non-significant points if needed)
    edgecolor="blue",    # green edge indicates significance
    linewidth=1,        # thickness of green edge
    label="P ≤ 0.05"
)

all_handles = size_handles + [sig_handle]

plt.legend(
    handles=all_handles,
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(0.83, 0.24),
    bbox_transform=plt.gcf().transFigure,
    title_fontsize=5,
    fontsize=5,
    labelspacing=0.8,     # ⬅️ increases vertical spacing between entries
    handletextpad=0.8,    # space between circle and text
    borderaxespad=0.0
)



# -----------------------
# Axes formatting
# -----------------------
ax.set_xticks(range(len(x_order)))
ax.set_xticklabels(x_order, rotation=45, ha="right")

ax.set_yticks(range(len(y_order)))
ax.set_yticklabels(y_order)

#ax.set_xlabel("PRS Model", labelpad=6)
ax.set_ylabel("Patient Cancer", labelpad=4)

# -----------------------
# Colorbar (tight)
# -----------------------
cbar = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
cbar.set_label(r"$-\log_{10}P$", fontsize=6)
cbar.ax.tick_params(labelsize=5, width=0.5, length=2)
cbar.outline.set_linewidth(0.5)

# -----------------------
# Manual layout (critical for 4×2)
# -----------------------
plt.subplots_adjust(
    left=0.2,
    right=0.9,
    bottom=0.27,
    top=0.98
)

# -----------------------
# Save
# -----------------------
out = "/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/pgs_cancer_circle_matrix"
plt.savefig(out + ".png", dpi=300)
plt.savefig(out + ".pdf")
plt.close()
