import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from itertools import product

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------
df = pd.read_csv("dfci-ufc.aou.phenos.v2.tsv", sep="\t")
df = df[df['original_dx'] != "control"]
# Rename original_id → Sample
df = df.rename(columns={"original_id": "Sample"})

# Drop rows with missing sex
#df = df[df["reported_sex"].isin(["male", "female"])].copy()

# ------------------------------------------------------------
# 1b. Remove samples with known pathogenic variants
# ------------------------------------------------------------
pvs_samples = set(
    pd.read_csv(
        "samples_with_pvs.nov7.list",
        header=None
    )[0].astype(str)
)

# Filter out those samples
df = df[~df["Sample"].astype(str).isin(pvs_samples)].copy()

# ------------------------------------------------------------
# 2. Cancer definitions
# ------------------------------------------------------------
PATIENT_CANCERS = {
    "Breast", "Prostate", "Colorectal", "Cervix", "Thyroid",
    "Sarcoma", "Neuroendocrine", "Lung", "Bladder", "Kidney",
    "Uterus", "Ovary", "Hematologic", "Non-Hodgkin",
    "Basal_Cell_Carcinoma", "Squamous_Cell_Carcinoma",
    "Melanoma", "Brain"
}

FAMILY_CANCERS = {
    "Bladder", "Blood_Soft_Tissue", "Bone", "Brain", "Cervix",
    "Colorectal", "Kidney", "Lung", "Ovary", "Prostate",
    "Skin", "Thyroid", "Uterus", "Breast"
}

FEMALE_ONLY = {"Breast", "Uterus", "Ovary", "Cervix"}
MALE_ONLY = {"Prostate"}

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------
def split_dx(val):
    if pd.isna(val):
        return set()
    return {
        x.strip()
        for x in str(val).split(";")
        if x.strip() != "Male_Breast"
    }

def has_cancer(dx_set, cancer):
    return cancer in dx_set

# Parse diagnosis columns
df["patient_dx"] = df["original_dx"].apply(split_dx)
df["maternal_dx"] = df["maternal_family_dx"].apply(split_dx)
df["paternal_dx"] = df["paternal_family_dx"].apply(split_dx)

df["family_dx"] = df.apply(
    lambda r: r["maternal_dx"] | r["paternal_dx"],
    axis=1
)

# ------------------------------------------------------------
# 4. Fisher exact test scan
# ------------------------------------------------------------
results = []

for patient_cancer, family_cancer in product(PATIENT_CANCERS, FAMILY_CANCERS):

    sub = df.copy()

    # Sex restrictions
    if patient_cancer in FEMALE_ONLY:
        sub = sub[sub["reported_sex"] == "female"]
    if patient_cancer in MALE_ONLY:
        sub = sub[sub["reported_sex"] == "male"]

    if sub.empty:
        continue

    # Boolean vectors
    patient_has = sub["patient_dx"].apply(
        lambda x: has_cancer(x, patient_cancer)
    )
    family_has = sub["family_dx"].apply(
        lambda x: has_cancer(x, family_cancer)
    )

    # Contingency table
    a = int((patient_has & family_has).sum())
    b = int((patient_has & ~family_has).sum())
    c = int((~patient_has & family_has).sum())
    d = int((~patient_has & ~family_has).sum())

    # Require minimal signal
    if a < 3:
        continue

    table = [[a, b], [c, d]]

    try:
        odds_ratio, p_value = fisher_exact(table, alternative="greater")
    except Exception:
        continue

    results.append({
        "patient_cancer": patient_cancer,
        "family_cancer": family_cancer,
        "n_intersection": a,
        "odds_ratio": odds_ratio,
        "p_value": p_value,
        "n_patients_tested": int(len(sub))
    })

# ------------------------------------------------------------
# 5. Output
# ------------------------------------------------------------
out = pd.DataFrame(results)

out = out.sort_values(
    by=["p_value", "odds_ratio"],
    ascending=[True, False]
)

out.to_csv(
    "patient_family_fisher_results.tsv",
    sep="\t",
    index=False
)

print(f"Wrote {len(out)} rows to patient_family_fisher_results.tsv")

