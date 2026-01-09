#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------
# Fixed cancer lists (alphabetical)
# -----------------------
PATIENT_CANCERS = [
    "Basal_Cell_Carcinoma","Breast","Prostate","Squamous_Cell_Carcinoma","Melanoma",
    "Lung","Colorectal","Thyroid","Hematologic","Non-Hodgkin","Uterus","Bladder","Sarcoma","Cervix","Ovary","Kidney",
    "Brain","Neuroendocrine","Control"
]

FAMILY_CANCERS = [
    "Skin","Breast","Prostate","Lung","Colorectal","Blood_Soft_Tissue","Brain","Bladder","Ovary","Thyroid","Bone","Cervix",
    "Kidney","Uterus"
]

# -----------------------
# Load data
# -----------------------
df = pd.read_csv(
    "dfci-ufc.aou.phenos.v2.tsv",
    sep="\t",
    dtype=str
)

# Remove excluded individuals
remove_ids = pd.read_csv(
    "samples_to_exclude.jan29.list",
    header=None,
    dtype=str
)[0].str.strip()

df = df[~df["Sample"].astype(str).isin(remove_ids)]

# Normalize strings
for col in ["original_dx", "maternal_family_dx", "paternal_family_dx"]:
    df[col] = df[col].fillna("").str.lower()

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
#   rows = PATIENT
#   cols = FAMILY
# -----------------------
count_matrix = pd.DataFrame(
    0, index=PATIENT_CANCERS, columns=FAMILY_CANCERS, dtype=int
)

jaccard_matrix = pd.DataFrame(
    0.0, index=PATIENT_CANCERS, columns=FAMILY_CANCERS
)

patient_masks = {
    pc: df["original_dx"].apply(lambda x: has_cancer(x, pc))
    for pc in PATIENT_CANCERS
}

family_masks = {
    fc: df.apply(lambda r: family_has(r, fc), axis=1)
    for fc in FAMILY_CANCERS
}

patient_sizes = {pc: m.sum() for pc, m in patient_masks.items()}
family_sizes = {fc: m.sum() for fc, m in family_masks.items()}

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
    f"{pc.replace('_',' ').replace('Non-Hodgkin', 'Non-Hodgkin Lymphoma')}\n"
    f"(n={patient_sizes[pc]})"
    for pc in PATIENT_CANCERS
]

x_labels = [
    f"{fc.replace('Blood_Soft_Tissue','Blood & Soft Tissue').replace('_',' ')}\n"
    f"(n={family_sizes[fc]})"
    for fc in FAMILY_CANCERS
]

# -----------------------
# Plot
# -----------------------
plt.figure(figsize=(18, 12))

ax = sns.heatmap(
    jaccard_matrix,
    cmap="Reds",
    linewidths=0.5,
    linecolor="lightgray",
    cbar_kws={"label": "Jaccard Index"},
    annot=False
)

# Get colorbar
cbar = ax.collections[0].colorbar

cbar.set_label("Jaccard Index", fontsize=16, fontweight="bold")
cbar.ax.tick_params(labelsize=13, width=2, length=8)
cbar.outline.set_linewidth(2)


# Count / dot overlay
for i, pc in enumerate(PATIENT_CANCERS):      # y-axis
    for j, fc in enumerate(FAMILY_CANCERS):   # x-axis
        val = count_matrix.loc[pc, fc]

        if val == 0:
            txt = "0"
        elif 1 <= val <= 20:
            txt = "·"
        else:
            txt = str(val)
            color = "white" if val > 190 else "black"

        ax.text(
            j + 0.5, i + 0.5,
            txt,
            ha="center", va="center",
            fontsize=12,
            color=color
        )

ax.set_xticklabels(x_labels, rotation=45, ha="right",fontsize=12)
ax.set_yticklabels(y_labels, rotation=0,fontsize=12)

#plt.title("Patient Cancer vs Family Cancer (Jaccard Index)")
#plt.xlabel("Family Cancer History",fontsize=16,fontweight="bold")
ax.set_xlabel("Family Cancer History", fontsize=16, fontweight="bold", labelpad=12)
plt.ylabel("Patient Cancer",fontsize=16,fontweight="bold")
plt.tight_layout()
plt.savefig("family_patient_cancer_matrix.png", dpi=300)
plt.close()

