#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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
df_all = df_all.fillna("0")

# -----------------------
# Separate totals from matrix
# -----------------------
col_totals = df_all.iloc[-1, :-1].astype(int)
row_totals = df_all.iloc[:-1, -1].astype(int)

col_totals_female = [437,352,219,218,195,96,119,76,63,65,52,60,52,29]
col_totals_male = [336,260,239,170,176,105,78,70,63,47,53,43,43,16]

df_counts = df_all.iloc[:-1, :-1].astype(int)

jaccard = pd.DataFrame(index=df_counts.index, columns=df_counts.columns, dtype=float)
# for i in df_counts.index:
#     for j in df_counts.columns:
#         inter = df_counts.loc[i, j]
#         union = row_totals[i] + col_totals[j] - inter
#         jaccard.loc[i, j] = (inter / union if union > 0 else 0)

# Define sex-specific cancer groups
female_cancers = {"breast", "uterus", "ovary", "cervix"}
male_cancers = {"prostate"}

# Convert lists to dicts keyed by column name (IMPORTANT)
col_totals_female_dict = dict(zip(df_counts.columns, col_totals_female))
col_totals_male_dict = dict(zip(df_counts.columns, col_totals_male))

for i in df_counts.index:
    for j in df_counts.columns:

        inter = df_counts.loc[i, j]

        # Choose correct column total
        if i in female_cancers:
            col_total = col_totals_female_dict[j]
        elif i in male_cancers:
            col_total = col_totals_male_dict[j]
        else:
            col_total = col_totals[j]  # fallback (your original)

        union = row_totals[i] + col_total - inter
        jaccard.loc[i, j] = inter / union if union > 0 else 0

# Scale colors by max Jaccard across all cells
vmax = jaccard.max().max()

# -----------------------
# Label cleanup
# -----------------------
def clean_label(x):
    x = str(x).capitalize()
    replacements = {
        "Basal_cell": "BCC",
        "Squamous_cell": "SCC",
        "Non-hodgkin": "NHL",
        "Neuroendocrine": "NETs",
        "Blood_soft_tissue": "B & ST"
    }
    for k,v in replacements.items():
        x = x.replace(k, v)
    x = x.replace("_", " ")
    return x#.title()

y_labels = [f"{clean_label(i)}\n(n={row_totals[i]})" for i in df_counts.index]
x_labels = [f"{clean_label(j)}\n(n={col_totals[j]})" for j in df_counts.columns]

# -----------------------
# Plot heatmap
# -----------------------
fig, ax = plt.subplots(figsize=(5.1, 4))
# sns.heatmap(
#     jaccard,
#     cmap="Reds",
#     linewidths=0.4,
#     linecolor="lightgray",
#     annot=df_counts,
#     fmt="d",
#     vmax=vmax,
#     ax=ax
# )

import numpy as np
# Create annotation labels
annot_labels = np.where((df_counts <= 20) & (df_counts > 0), ".", df_counts.astype(int).astype(str))

# sns.heatmap(
#     jaccard,
#     cmap="Reds",
#     linewidths=0.4,
#     linecolor="lightgray",
#     annot=annot_labels,   # use modified labels
#     fmt="",               # IMPORTANT: since we're passing strings
#     vmax=vmax,
#     ax=ax
# )

sns.heatmap(
    jaccard,
    cmap="Reds",
    linewidths=0.4,
    linecolor="lightgray",
    annot=annot_labels,
    fmt="",
    vmax=vmax,
    ax=ax,
    cbar_kws={
        "label": "Jaccard Index",  # add title
        "pad": 0.02,               # default ~0.05, smaller = closer
        "aspect": 40     # thinner bar
    }
)
# -----------------------
# Highlight matching cells
# -----------------------
CYAN = "#007C7C"

# -----------------------
# Fixed cancer lists (alphabetical)
# -----------------------
PATIENT_CANCERS = [
    "Basal_Cell_Carcinoma","Breast","Prostate","Squamous_Cell_Carcinoma","Melanoma",
    "Lung","Colorectal","Thyroid","Hematologic","Non-Hodgkin","Uterus","Bladder","Sarcoma",
    "Ovary","Kidney","Cervix","Neuroendocrine","Brain","Control"
]

FAMILY_CANCERS = [
    "Skin","Breast","Prostate","Lung","Colorectal","Blood_soft_tissue","Brain","Bladder",
    "Ovary","Thyroid","Bone","Cervix","Kidney","Uterus"
]

# -----------------------
# Patient → family match map
# -----------------------
MATCH_MAP = {
    "Basal_Cell_Carcinoma": ["Skin"],
    "Squamous_Cell_Carcinoma": ["Skin"],
    "Melanoma": ["Skin"],
    "Non-Hodgkin": ["Blood_soft_tissue"],
    "Hematologic": ["Blood_soft_tissue"],
    "Sarcoma":["Bone"],
    "Breast":["Breast"],
    "Lung":["Lung"],
    "Prostate":["Prostate"],
    "Thyroid":["Thyroid"],
    "Ovary":["Ovary"],
    "Kidney":["Kidney"],
    "Brain":["Brain"],
    "Uterus":["Uterus"],
    "Colorectal":["Colorectal"],
    "Cervix":["Cervix"]

}
from matplotlib.patches import Rectangle
for i, pc in enumerate(PATIENT_CANCERS):
    if pc not in MATCH_MAP:
        continue
    for fc in MATCH_MAP[pc]:
        if fc not in FAMILY_CANCERS:
            continue
        j = FAMILY_CANCERS.index(fc)
        ax.add_patch(
            Rectangle(
                (j, i), 1, 1,
                fill=False,
                edgecolor=CYAN,
                linewidth=2.2
            )
        )


legend_patch = Rectangle(
    (0, 0), 1, 1,
    fill=False,
    edgecolor=CYAN,
    linewidth=1.5,
)
ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=5)
ax.set_yticklabels(y_labels, rotation=0, fontsize=5)
ax.set_xlabel("Cancer History of First-Degree Relatives", fontsize=7, fontweight="bold")
ax.set_ylabel("18 Cancer Cohorts in This Study", fontsize=7, fontweight="bold")

plt.tight_layout()
plt.savefig("/Users/noah/Desktop/ufc_repository/results/Figure1/patient_family_matrix.pdf",pad_inches=0, bbox_inches="tight")
#plt.show()