import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# --- Load files ---
pval_df = pd.read_csv("/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.pvalue.tsv", sep="\t")
or_df = pd.read_csv("/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.OR.tsv", sep="\t")

# --- Keep these PGS IDs ---
pgs_keep = [
    "PGS003386","PGS000789","PGS000791","PGS000784",
    "PGS004690","PGS003382","PGS004694","PGS000797","PGS004687",
    "PGS003384","PGS000788","PGS000793","PGS004244","PGS004688"
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
plt.colorbar(label="-log10(p-value)")
plt.yticks(fontsize=14)
#plt.xticks(fontsize=14,rotation=45, ha='right')
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
plt.savefig("pgs_cancer_circle_matrix.png", dpi=300)
plt.show()
