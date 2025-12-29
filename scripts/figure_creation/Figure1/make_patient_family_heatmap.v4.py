#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------
# Fixed cancer lists (alphabetical)
# -----------------------
PATIENT_CANCERS = sorted([
    "Basal_Cell_Carcinoma","Squamous_Cell_Carcinoma","Melanoma",
    "Breast","Male_Breast","Prostate","Lung","Colorectal","Kidney",
    "Hematologic","Non-Hodgkin","Thyroid","Ovary","Pancreas","Bladder",
    "Brain","Cervix","Uterus","Neuroendocrine","Sarcoma","Control"
])

FAMILY_CANCERS = sorted([
    "Breast","Male_Breast","Prostate","Lung","Cervix","Colorectal",
    "Kidney","Thyroid","Ovary","Pancreas","Bladder","Brain","Bone",
    "Skin","Uterus","Blood_Soft_Tissue"
])

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
    "samples_with_pvs.nov7.list",
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
# -----------------------
count_matrix = pd.DataFrame(
    0, index=FAMILY_CANCERS, columns=PATIENT_CANCERS, dtype=int
)
jaccard_matrix = pd.DataFrame(
    0.0, index=FAMILY_CANCERS, columns=PATIENT_CANCERS
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

for fc in FAMILY_CANCERS:
    for pc in PATIENT_CANCERS:
        intersection = (patient_masks[pc] & family_masks[fc]).sum()
        union = (patient_masks[pc] | family_masks[fc]).sum()

        count_matrix.loc[fc, pc] = intersection
        jaccard_matrix.loc[fc, pc] = intersection / union if union > 0 else 0

# -----------------------
# Axis labels with n
# -----------------------
x_labels = [f"{pc}\n(n={patient_sizes[pc]})" for pc in PATIENT_CANCERS]
y_labels = [f"{fc}\n(n={family_sizes[fc]})" for fc in FAMILY_CANCERS]

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

# Dot/count overlay
for i, fc in enumerate(FAMILY_CANCERS):
    for j, pc in enumerate(PATIENT_CANCERS):
        val = count_matrix.loc[fc, pc]
        if val == 0:
            txt = "0"
        elif 1 <= val <= 20:
            txt = "·"
        else:
            txt = str(val)

        ax.text(j + 0.5, i + 0.5, txt,
                ha="center", va="center", fontsize=9)

ax.set_xticklabels(x_labels, rotation=45, ha="right")
ax.set_yticklabels(y_labels, rotation=0)

plt.title("Family Cancer vs Patient Cancer (Jaccard Index)")
plt.xlabel("Patient Cancer")
plt.ylabel("Family Cancer History")
plt.tight_layout()
plt.savefig("family_patient_cancer_matrix.png",dpi=300)
