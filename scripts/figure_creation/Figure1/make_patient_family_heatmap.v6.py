#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

# -----------------------
# Global font settings (Arial only)
# -----------------------
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 6,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
})

# -----------------------
# Fixed cancer lists (alphabetical)
# -----------------------
PATIENT_CANCERS = [
    "Basal_Cell_Carcinoma","Breast","Prostate","Squamous_Cell_Carcinoma","Melanoma",
    "Lung","Colorectal","Thyroid","Hematologic","Non-Hodgkin","Uterus","Bladder","Sarcoma",
    "Cervix","Ovary","Kidney","Brain","Neuroendocrine","Control"
]

FAMILY_CANCERS = [
    "Skin","Breast","Prostate","Lung","Colorectal","Blood_Soft_Tissue","Brain","Bladder",
    "Ovary","Thyroid","Bone","Cervix","Kidney","Uterus"
]

# Sex-restricted cancers
SEX_RESTRICTED = {
    "Breast": "Female",
    "Cervix": "Female",
    "Ovary": "Female",
    "Uterus": "Female",
    "Prostate": "Male",
}

# -----------------------
# Load data
# -----------------------
df = pd.read_csv(
    "/Users/noah/Desktop/ufc_repository/results/epidemiological_results/merged_sex_dx_family.tsv.gz",
    sep="\t",
    dtype=str
)

# Remove excluded individuals
# remove_ids = pd.read_csv(
#     "samples_to_exclude.jan29.list",
#     header=None,
#     dtype=str
# )[0].str.strip()

# df = df[~df["Sample"].isin(remove_ids)]

# Normalize strings
for col in ["original_dx", "maternal_family_dx", "paternal_family_dx"]:
    df[col] = df[col].fillna("").str.lower()

df["inferred_sex"] = df["inferred_sex"].str.capitalize()

# -----------------------
# Helper functions
# -----------------------
def split_dx(dx):
    return [x.strip() for x in dx.split(";") if x.strip()]

def has_cancer(dx_string, cancer):
    return cancer.lower() in split_dx(dx_string)

def family_has_any(row):
    return bool(split_dx(row["maternal_family_dx"]) or
                split_dx(row["paternal_family_dx"]))

def family_has(row, cancer):
    return (
        has_cancer(row["maternal_family_dx"], cancer) or
        has_cancer(row["paternal_family_dx"], cancer)
    )

# -----------------------
# Remove controls with ANY family cancer
# -----------------------
is_control = df["original_dx"].apply(lambda x: has_cancer(x, "control"))
has_fam_any = df.apply(family_has_any, axis=1)
df = df[~(is_control & has_fam_any)]

# -----------------------
# Build matrices
# -----------------------
count_matrix = pd.DataFrame(
    0, index=PATIENT_CANCERS, columns=FAMILY_CANCERS, dtype=int
)

jaccard_matrix = pd.DataFrame(
    0.0, index=PATIENT_CANCERS, columns=FAMILY_CANCERS
)

# Patient masks with sex restriction
patient_masks = {}
for pc in PATIENT_CANCERS:
    base_mask = df["original_dx"].apply(lambda x: has_cancer(x, pc))

    if pc in SEX_RESTRICTED:
        sex = SEX_RESTRICTED[pc]
        base_mask &= df["inferred_sex"] == sex

    patient_masks[pc] = base_mask

# Family masks (no sex restriction)
family_masks = {
    fc: df.apply(lambda r: family_has(r, fc), axis=1)
    for fc in FAMILY_CANCERS
}

patient_sizes = {pc: m.sum() for pc, m in patient_masks.items()}
family_sizes = {fc: m.sum() for fc, m in family_masks.items()}

# -----------------------
# Jaccard calculation
# -----------------------
for pc in PATIENT_CANCERS:
    for fc in FAMILY_CANCERS:
        intersection = (patient_masks[pc] & family_masks[fc]).sum()
        union = (patient_masks[pc] | family_masks[fc]).sum()

        count_matrix.loc[pc, fc] = intersection
        jaccard_matrix.loc[pc, fc] = intersection / union if union > 0 else 0

# -----------------------
# Axis labels with n
# -----------------------
y_labels = [
    f"{pc.replace('_',' ').replace('Non-Hodgkin', 'NHL').replace('Basal Cell Carcinoma','BCC').replace('Squamous Cell Carcinoma','SCC')}\n"
    f"(n={patient_sizes[pc]})"
    for pc in PATIENT_CANCERS
]

x_labels = [
    f"{fc.replace('Blood_Soft_Tissue','B & ST').replace('_',' ')}\n"
    f"(n={family_sizes[fc]})"
    for fc in FAMILY_CANCERS
]

# -----------------------
# Plot
# -----------------------
fig, ax = plt.subplots(figsize=(5.1, 4))

sns.heatmap(
    jaccard_matrix,
    cmap="Reds",
    linewidths=0.4,
    linecolor="lightgray",
    cbar_kws={"label": "Jaccard Index"},
    ax=ax
)

# Colorbar formatting
cbar = ax.collections[0].colorbar
cbar.set_label("Jaccard Index", fontsize=5, fontweight="bold")
cbar.ax.tick_params(labelsize=6, width=1, length=4)
cbar.outline.set_linewidth(1)


# Count / dot overlay
for i, pc in enumerate(PATIENT_CANCERS):
    for j, fc in enumerate(FAMILY_CANCERS):
        val = count_matrix.loc[pc, fc]

        if val == 0:
            txt = "0"
            color = "black"
        elif val <= 20:
            txt = "·"
            color = "black"
        else:
            txt = str(val)
            color = "white" if val > 190 else "black"

        ax.text(
            j + 0.5, i + 0.5,
            txt,
            ha="center", va="center",
            fontsize=5,
            color=color
        )

ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=5)
ax.set_yticklabels(y_labels, rotation=0, fontsize=5)
ax.tick_params(axis="x", pad=0)   # default is ~3–4


ax.set_xlabel("Family Cancer History", fontsize=7, fontweight="bold", labelpad=-1)
ax.set_ylabel("Patient Cancer", fontsize=7, fontweight="bold", labelpad=1)

# Tighten layout (especially right side past colorbar)
plt.subplots_adjust(left=0.15, right=1.01, bottom=0.12, top=0.99)

plt.savefig("/Users/noah/Desktop/ufc_repository/results/epidemiological_results/family_patient_cancer_matrix.png", dpi=300, facecolor="white")
plt.savefig("/Users/noah/Desktop/ufc_repository/results/epidemiological_results/family_patient_cancer_matrix.pdf", dpi=300, facecolor="white")
plt.close()
