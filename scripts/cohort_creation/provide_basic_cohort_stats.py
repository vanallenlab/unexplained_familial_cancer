import pandas as pd

cancer_cols = [
    "breast_cancer", "leukemia_cancer", "prostate_cancer", "kidney_cancer",
    "bladder_cancer", "colon_cancer", "thyroid_cancer", "ovarian_cancer",
    "endometrial_cancer", "melanoma_cancer", "nhl_cancer", "lung_cancer"
]

metadata = pd.read_csv("CDRv7.ufc_original.metadata", sep="\t", dtype={"person_id": str})

# Count individuals with at least one listed cancer
cancer_status = metadata[cancer_cols].fillna(0).eq(1)
has_cancer = cancer_status.any(axis=1)
multiple_cancers = cancer_status.sum(axis=1) > 1

print(f"Individuals with at least one listed cancer: {metadata.loc[has_cancer, 'person_id'].nunique():,}")
print(f"Individuals with more than one listed cancer: {metadata.loc[multiple_cancers, 'person_id'].nunique():,}")

# Count individuals with pancancer_cancer == 0 or NA
pancancer_zero_or_na = metadata["pancancer_cancer"].isna() | metadata["pancancer_cancer"].eq(0)
print(f"Individuals with pancancer_cancer == 0 or NA: {metadata.loc[pancancer_zero_or_na, 'person_id'].nunique():,}")

# Retain individuals with at least one listed cancer OR no pancancer diagnosis
keep = has_cancer | pancancer_zero_or_na
study_metadata = metadata.loc[keep].copy()

print(f"Total individuals in study cohort: {study_metadata['person_id'].nunique():,}")

# Define cases and controls within the study cohort
study_case_ids = set(
    study_metadata.loc[study_metadata[cancer_cols].fillna(0).eq(1).any(axis=1), "person_id"].dropna()
)
study_control_ids = set(
    study_metadata.loc[study_metadata["pancancer_cancer"].isna() | study_metadata["pancancer_cancer"].eq(0), "person_id"].dropna()
)

print(f"Study cases (at least one listed cancer): {len(study_case_ids):,}")
print(f"Study controls (pancancer_cancer == 0 or NA): {len(study_control_ids):,}")

study_metadata.to_csv("CDRv7.ufc_12cancers.metadata", sep="\t", index=False)

# Read SV QC-passing sample list
sv_samples = pd.read_csv(
    "all_ufc_samples.march4.sv_qc_pass.2695.list",
    header=None,
    dtype=str
)[0].str.strip()

sv_sample_set = set(sv_samples.dropna())

# Calculate SV availability in the study cohort, split into cases and controls
study_ids = set(study_metadata["person_id"].dropna())
study_with_sv = study_ids & sv_sample_set

cases_with_sv = study_case_ids & sv_sample_set
controls_with_sv = study_control_ids & sv_sample_set

n_study = len(study_ids)
n_with_sv = len(study_with_sv)
percent_with_sv = (n_with_sv / n_study * 100) if n_study else 0

n_cases = len(study_case_ids)
n_controls = len(study_control_ids)

print(f"\nStudy individuals with SV QC-passing calls: {n_with_sv:,} / {n_study:,} ({percent_with_sv:.2f}%)")
print(
    f"  Cases with SV QC-passing calls: {len(cases_with_sv):,} / {n_cases:,} "
    f"({len(cases_with_sv) / n_cases * 100 if n_cases else 0:.2f}%)"
)
print(
    f"  Controls with SV QC-passing calls: {len(controls_with_sv):,} / {n_controls:,} "
    f"({len(controls_with_sv) / n_controls * 100 if n_controls else 0:.2f}%)"
)