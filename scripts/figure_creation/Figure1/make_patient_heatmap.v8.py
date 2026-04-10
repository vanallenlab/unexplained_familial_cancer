#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib as mpl

# -----------------------
# Global font settings
# -----------------------
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 6,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
})

# -----------------------
# Load matrix file
# -----------------------
file_path = "/Users/noah/Desktop/ufc_repository/results/epidemiological_results/family_history_matrix.tsv"
df_all = pd.read_csv(file_path, sep="\t", dtype=str, index_col=0)

# Fill missing cells with 0
df_all = df_all.fillna("0")

# -----------------------
# Separate totals from matrix
# -----------------------
# The last row = column totals for labels
col_totals = df_all.iloc[-1, :-1].astype(int)

# The last column = row totals for labels
row_totals = df_all.iloc[:-1, -1].astype(int)

# The heatmap matrix = all rows except last, all columns except last
df = df_all.iloc[:-1, :-1].astype(int)

# Ensure all indices/columns are strings
df.index = df.index.astype(str)
df.columns = df.columns.astype(str)

# -----------------------
# Pretty label function
# -----------------------
def clean_label(x):
    x = str(x)
    return (
        x.replace("_", " ")
         .replace("non-hodgkin", "NHL")
         .replace("basal cell carcinoma", "BCC")
         .replace("squamous cell carcinoma", "SCC")
         .replace("neuroendocrine", "NETs")
         .title()
    )

y_labels = [f"{clean_label(i)}\n(n={row_totals[i]})" for i in df.index]
x_labels = [f"{clean_label(j)}\n(n={col_totals[j]})" for j in df.columns]

# -----------------------
# Plot heatmap
# -----------------------
fig, ax = plt.subplots(figsize=(8,6))
sns.heatmap(df, cmap="Reds", linewidths=0.4, linecolor="lightgray", annot=True, fmt="d", ax=ax)

ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=6)
ax.set_yticklabels(y_labels, rotation=0, fontsize=6)

ax.set_xlabel("Family Cancer History", fontsize=7, fontweight="bold")
ax.set_ylabel("Patient Cancer", fontsize=7, fontweight="bold")

plt.tight_layout()
plt.savefig("/Users/noah/Downloads/output.pdf")
plt.show()