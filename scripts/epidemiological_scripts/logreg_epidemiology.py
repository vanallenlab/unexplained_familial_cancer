import pandas as pd
import numpy as np
import statsmodels.api as sm
from itertools import product

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------
df = pd.read_csv("dfci-ufc.aou.phenos.v2.tsv", sep="\t")
df = df[df['original_dx'] != "control"]
pheno_df = pd.read_csv("dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz",sep='\t',usecols = ['original_id','inferred_sex','cohort','intake_qc_pop'])
pheno_df = pheno_df[pheno_df['cohort'] == "aou"]
df['Sample'] = df['Sample'].astype(str)
pheno_df['Sample'] = pheno_df['original_id'].astype(str)
df = df.merge(pheno_df, on="Sample",how="left")
df = df[df["inferred_sex"].isin(["male", "female"])].copy()
df["sex_binary"] = (df["inferred_sex"] == "female").astype(int)
df = df[df['intake_qc_pop'] == 'EUR']
# ------------------------------------------------------------
# 1b. Remove samples with known pathogenic variants
# ------------------------------------------------------------
pvs_samples = set(
    pd.read_csv(
        "samples_to_exclude.jan29.list",
        header=None
    )[0].astype(str)
)

print(df.columns)
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

df["patient_dx"] = df["original_dx"].apply(split_dx)
df["maternal_dx"] = df["maternal_family_dx"].apply(split_dx)
df["paternal_dx"] = df["paternal_family_dx"].apply(split_dx)

df["family_dx"] = df.apply(
    lambda r: r["maternal_dx"] | r["paternal_dx"],
    axis=1
)

# ------------------------------------------------------------
# 4. Logistic regression scan
# ------------------------------------------------------------
results = []

for patient_cancer, family_cancer in product(PATIENT_CANCERS, FAMILY_CANCERS):

    sub = df.copy()

    # Sex restrictions
    if patient_cancer in FEMALE_ONLY:
        sub = sub[sub["inferred_sex"] == "female"]
    if patient_cancer in MALE_ONLY:
        sub = sub[sub["inferred_sex"] == "male"]

    if sub.empty:
        continue

    # Binary outcome and predictor
    sub["patient_has"] = sub["patient_dx"].apply(
        lambda x: has_cancer(x, patient_cancer)
    ).astype(int)

    sub["family_has"] = sub["family_dx"].apply(
        lambda x: has_cancer(x, family_cancer)
    ).astype(int)

    # Require minimal signal (cases with family history)
    if sub.loc[sub["patient_has"] == 1, "family_has"].sum() < 2:
        continue

    # Design matrix
    covars = ["family_has"]
    if patient_cancer not in ("Breast","Prostate","Cervix","Uterus","Ovary"):
        covars = covars + ['sex_binary']
    X = sm.add_constant(sub[["family_has"]])
    y = sub["patient_has"]

    try:
        model = sm.Logit(y, X).fit(disp=False)
    except Exception:
        continue

    beta = model.params["family_has"]
    se = model.bse["family_has"]
    pval = model.pvalues["family_has"]

    odds_ratio = np.exp(beta)
    ci_low, ci_high = np.exp(
        model.conf_int().loc["family_has"]
    )

    results.append({
        "patient_cancer": patient_cancer,
        "family_cancer": family_cancer,
        "n_intersection": int(
            ((sub["patient_has"] == 1) & (sub["family_has"] == 1)).sum()
        ),
        "odds_ratio": odds_ratio,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "p_value": pval,
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
    "patient_family_logistic_results.tsv",
    sep="\t",
    index=False
)

print(f"Wrote {len(out)} rows to patient_family_logistic_results.tsv")

