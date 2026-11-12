import pandas as pd

METADATA_FILE = "CDRv7.ufc_original.metadata"

CANCERS = {
    "Breast": ("breast_cancer", "fdrs_w_breast_cancer"),
    "Prostate": ("prostate_cancer", "fdrs_w_prostate_cancer"),
    "Melanoma": ("melanoma_cancer", "fdrs_w_skin_cancer"),
    "Lung": ("lung_cancer", "fdrs_w_lung_cancer"),
    "Colorectal": ("colon_cancer", "fdrs_w_colon_cancer"),
    "Thyroid": ("thyroid_cancer", "fdrs_w_thyroid_cancer"),
    "Hematologic": ("leukemia_cancer", "fdrs_w_blood_and_soft_tissue_cancer"),
    "NHL": ("non_hodgkin_lymphoma_cancer", "fdrs_w_blood_and_soft_tissue_cancer"),
    "Endometrial": ("endometrial_cancer", "fdrs_w_endometrial_cancer"),
    "Bladder": ("bladder_cancer", "fdrs_w_bladder_cancer"),
    "Ovarian": ("ovarian_cancer", "fdrs_w_ovarian_cancer"),
    "Kidney": ("kidney_cancer", "fdrs_w_kidney_cancer")
}

SEX_RESTRICTIONS = {
    "Breast": "Female",
    "Prostate": "Male",
    "Endometrial": "Female",
    "Ovarian": "Female",
}

df = pd.read_csv(
    METADATA_FILE,
    sep="\t",
    dtype={"person_id": str},
)

# Missing cancer flags and family-history counts are treated as zero.
for cancer_col, fdr_col in CANCERS.values():
    df[cancer_col] = df[cancer_col].fillna(0)
    df[fdr_col] = df[fdr_col].fillna(0)

# Define controls using the pancancer phenotype.
controls = df[df["pancancer_cancer"].isna()]

cancer_names = list(CANCERS.keys())

matrix = pd.DataFrame(
    0,
    index=cancer_names + ["Controls"],
    columns=cancer_names,
    dtype=int,
)

# Count family-history overlap among individuals with each cancer.
for patient_cancer, (patient_col, _) in CANCERS.items():
    patient_df = df

    if patient_cancer in SEX_RESTRICTIONS:
        patient_df = df[
            df["sex_at_birth"] == SEX_RESTRICTIONS[patient_cancer]
        ]

    affected = patient_df[patient_df[patient_col] == 1]

    for family_cancer, (_, fdr_col) in CANCERS.items():
        matrix.loc[patient_cancer, family_cancer] = int(
            (affected[fdr_col] >= 1).sum()
        )

# Count family-history overlap among pancancer-negative controls.
for family_cancer, (_, fdr_col) in CANCERS.items():
    matrix.loc["Controls", family_cancer] = int(
        (controls[fdr_col] >= 1).sum()
    )

# Add row and column sums.
matrix["Total"] = matrix.sum(axis=1)
matrix.loc["Total"] = matrix.sum(axis=0)

matrix.to_csv("/Users/noah/Desktop/ufc_repository/results/Figure1/Figure_1D.tsv",sep='\t',index=False)

# print(matrix.to_string())
# print(f"\nSaved overlap matrix to: {OUTPUT_FILE}")
