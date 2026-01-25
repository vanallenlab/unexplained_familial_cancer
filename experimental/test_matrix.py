#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------
# 1. Load fake data
# -----------------------
df = pd.read_csv("fake_paf_data.tsv", sep="\t")

# -----------------------
# 2. Construct binary risk factors
# -----------------------
risk_df = pd.DataFrame({
    "paf_prs_ge_1": (df["paf_prs"] >= 1).astype(int),
    "paf_prs_ge_2": (df["paf_prs"] >= 2).astype(int),
    "paf_VUS": df["paf_VUS"].astype(int),
    "paf_SV": df["paf_SV"].astype(int),
    "paf_GeneX": df["paf_GeneX"].astype(int),
})

factors = risk_df.columns

# -----------------------
# 3. Compute Jaccard
# -----------------------
def jaccard(a, b):
    union = ((a == 1) | (b == 1)).sum()
    if union == 0:
        return 0.0
    intersection = ((a == 1) & (b == 1)).sum()
    return intersection / union

jaccard_mat = pd.DataFrame(index=factors, columns=factors, dtype=float)

for f1 in factors:
    for f2 in factors:
        j = jaccard(risk_df[f1], risk_df[f2])
        jaccard_mat.loc[f1, f2] = j

# -----------------------
# 4. Apply masking rules
# -----------------------
# Mask diagonal
np.fill_diagonal(jaccard_mat.values, np.nan)

# Mask PRS vs PRS comparisons
prs_cols = ["paf_prs_ge_1", "paf_prs_ge_2"]
for c1 in prs_cols:
    for c2 in prs_cols:
        jaccard_mat.loc[c1, c2] = np.nan

# Keep only bottom triangle
mask_upper = np.triu(np.ones_like(jaccard_mat, dtype=bool))
jaccard_mat = jaccard_mat.mask(mask_upper)

# -----------------------
# 5. Plot heatmap
# -----------------------
plt.figure(figsize=(7, 6))
sns.set(style="white")

ax = sns.heatmap(
    jaccard_mat,
    cmap="Greens",
    vmin=0,
    vmax=1,
    square=True,
    linewidths=0.5,
    linecolor="white",
    cbar_kws={"label": "Jaccard similarity"},
    annot=True,
    fmt=".2f"
)

ax.set_title("Co-occurrence of Genetic Risk Factors (Filtered)", fontsize=14)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

# Clean spines
for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()
plt.show()
