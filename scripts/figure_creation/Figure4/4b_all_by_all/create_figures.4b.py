import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# --- Load files ---
pval_df = pd.read_csv("/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.pvalue.tsv", sep="\t")
or_df = pd.read_csv("/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.OR.tsv", sep="\t")

# --- Keep these PGS IDs ---
pgs_keep = [
    "PGS000785","PGS000789","PGS000791","PGS000784",
    "PGS000787","PGS003382","PGS004694","PGS000797","PGS004687",
    "PGS003384","PGS000788","PGS004249","PGS004244","PGS004688"
]

def filter_pgs(df):
    df = df[~df["PGS_ID"].str.contains("-euro", case=False)]
    df = df[df["PGS_ID"].isin(pgs_keep)]
    return df

pval_df = filter_pgs(pval_df)
or_df = filter_pgs(or_df)

# --- Melt both DataFrames ---
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

# --- Extract numeric OR from string like '1.077 (0.954, 1.215)' ---
def parse_or(value):
    try:
        # Extract the first numeric part
        match = re.match(r"([\d\.Ee+-]+)", str(value))
        return float(match.group(1)) if match else np.nan
    except Exception:
        return np.nan

or_long["OR"] = or_long["OR_raw"].apply(parse_or)

# --- Merge ---
merged = pd.merge(pval_long, or_long, on=["Intended_Cancer_Type", "PGS_ID", "Cancer"])
# Option 1: reassign
merged['Intended_Cancer_Type'] = merged['Intended_Cancer_Type'].str.replace("_", " ")
merged['Cancer'] = merged['Cancer'].str.replace("_", " ").str.title()

# --- Create ordered axis labels ---
merged["x_label"] = (
    merged["Intended_Cancer_Type"].str.title()
    + " ("
    + merged["PGS_ID"]
    + ")"
)

# Sort labels alphabetically
x_order = sorted(merged["x_label"].unique())
y_order = sorted(merged["Cancer"].unique())

# Convert to ordered categoricals
merged["x_label"] = pd.Categorical(
    merged["x_label"], categories=x_order, ordered=True
)
merged["Cancer"] = pd.Categorical(
    merged["Cancer"], categories=y_order, ordered=True
)

# --- Transform values ---
merged["p_value"] = merged["p_value"].astype(float)
merged["neglog10p"] = -np.log10(merged["p_value"])
merged["neglog10p"] = merged["neglog10p"].clip(0, 10)
merged["OR"] = merged["OR"].clip(0, 5)

# --- Map categories to numeric positions ---
x_map = {label: i for i, label in enumerate(x_order)}
y_map = {label: i for i, label in enumerate(y_order)}

merged["x_pos"] = merged["x_label"].map(x_map)

# --- Plot ---
plt.figure(figsize=(14, 7))
# plt.scatter(
#     x=merged["Intended_Cancer_Type"].str.title() + " (" + merged["PGS_ID"] + ")",
#     y=merged["Cancer"],
#     s=merged["OR"] ** 1.5 * 80,          # circle size by OR
#     c=merged["neglog10p"],               # color by -log10(p)
#     cmap="Reds",
#     alpha=0.9,
#     linewidth=0.5
# )
plt.scatter(
    x=merged["x_pos"],
    y=merged["Cancer"],
    s=merged["OR"] ** 1.5 * 80,
    c=merged["neglog10p"],
    cmap="Reds",
    alpha=0.9,
    linewidth=0.5
)
plt.colorbar(label=r"$-\log_{10}P$")
plt.yticks(fontsize=14)
#plt.xticks(fontsize=14,rotation=45, ha='right')


### Plot different colors

sig = sig = (merged["p_value"] <= 0.05) & (merged["OR"] > 1)
nonsig = merged["p_value"] > 0.05

# Base layer: all points, NO edges
plt.scatter(
    x=merged["x_pos"],
    y=merged["Cancer"],
    s=merged["OR"] ** 1.5 * 80,
    c=merged["neglog10p"],
    cmap="Reds",
    alpha=0.9,
    linewidth=0,
    edgecolors="none",
    zorder=1
)

# Overlay: significant points only, green edge
plt.scatter(
    x=merged.loc[sig, "x_pos"],
    y=merged.loc[sig, "Cancer"],
    s=merged.loc[sig, "OR"] ** 1.5 * 80,
    c=merged.loc[sig, "neglog10p"],
    cmap="Reds",
    alpha=0.9,
    linewidth=1.8,
    edgecolors="green",
    zorder=2
)
###

#### Added OR scaling value ####
or_vals = [0.5, 1, 2]

size_handles = [
    plt.scatter(
        [], [], 
        s=or_val ** 1.5 * 80,
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
    s=100,                # choose a size similar to your points
    c="white",            # fill color (matches non-significant points if needed)
    edgecolor="green",    # green edge indicates significance
    linewidth=1.2,        # thickness of green edge
    label="P ≤ 0.05 & OR > 1"
)

all_handles = size_handles + [sig_handle]

plt.legend(
    handles=all_handles,
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(0.80, 0.35),  # adjust x,y manually
    bbox_transform=plt.gcf().transFigure,
    title_fontsize=10
)



# Add it to the legend (can be separate from size legend)
# plt.legend(
#     handles=[sig_handle],
#     title="Significance",
#     frameon=False,
#     loc="upper left",
#     bbox_to_anchor=(0.82, 0.15),  # adjust below colorbar
#     bbox_transform=plt.gcf().transFigure,
#     title_fontsize=10
# )

# plt.legend(
#     handles=size_handles,
#     title="Odds Ratio",
#     frameon=False,
#     loc="lower right",
#     bbox_to_anchor=(1.35, 1)
# )
#### Added OR scaling value ####

# --- Apply ordered ticks ---
plt.xticks(
    ticks=range(len(x_order)),
    labels=x_order,
    rotation=45,
    ha="right",
    fontsize=14
)
plt.xlabel("PRS Inteneded Cancer",fontsize=20)
plt.title("Polygenic Scores Across Cancers", fontsize=16, fontweight='bold')
plt.grid(False)
plt.tight_layout()
plt.savefig("/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/pgs_cancer_circle_matrix.png", dpi=300)
plt.show()
