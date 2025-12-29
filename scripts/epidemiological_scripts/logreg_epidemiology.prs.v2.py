import pandas as pd
import numpy as np
import statsmodels.api as sm
from itertools import product
import os

# ============================================================
# Configuration
# ============================================================
RAW_PGS_DIR = "raw_pgs"
MAP_FILE = os.path.join(RAW_PGS_DIR, "cancer_pgs.map")

# ============================================================
# 1. Load cancer → PGS map
# ============================================================
pgs_map = pd.read_csv(
    MAP_FILE,
    sep=r"\s+",
    header=None,
    names=["cancer", "pgs_id"]
)

pgs_map["cancer"] = pgs_map["cancer"].str.lower()

# Melanoma equivalents
#MELANOMA_EQUIV = {
#    "melanoma": "skin",
#    "melanoma": "basal_cell_carcinoma",
#    "melanoma": "squamous_cell_carcinoma"
#}

#pgs_map["cancer"] = pgs_map["cancer"].replace(MELANOMA_EQUIV)
PGS_LOOKUP = dict(zip(pgs_map["cancer"], pgs_map["pgs_id"]))

# ============================================================
# 2. Load phenotype table
# ============================================================
df = pd.read_csv("dfci-ufc.aou.phenos.v2.tsv", sep="\t")
df = df[df["original_dx"] != "control"].copy()

df["original_id"] = df["Sample"].astype(str)

# ============================================================
# 3. Remove pathogenic variant samples
# ============================================================
pvs = set(
    pd.read_csv("samples_with_pvs.nov7.list", header=None)[0].astype(str)
)
df = df[~df["original_id"].isin(pvs)].copy()

# ============================================================
# 4. Sex filtering + encoding
# ============================================================
df = df[df["reported_sex"].isin(["male", "female"])].copy()
df["sex_binary"] = (df["reported_sex"] == "female").astype(int)

# ============================================================
# 5. Normalize diagnosis strings (NO SETS)
# ============================================================
def normalize_dx(val):
    if pd.isna(val) or isinstance(val, bool):
        return ""
    return str(val).lower()

df["patient_dx"] = df["original_dx"].apply(normalize_dx)
df["maternal_dx"] = df["maternal_family_dx"].apply(normalize_dx)
df["paternal_dx"] = df["paternal_family_dx"].apply(normalize_dx)

# ============================================================
# 6. Cancer definitions
# ============================================================
PATIENT_CANCERS = {
    "breast", "prostate", "colorectal", "cervix", "thyroid",
    "sarcoma", "neuroendocrine", "lung", "bladder", "kidney",
    "uterus", "ovary", "hematologic", "non-hodgkin",
    "basal_cell_carcinoma", "squamous_cell_carcinoma",
    "melanoma", "brain"
}

FAMILY_CANCERS = {
    "bladder", "blood_soft_tissue", "bone", "brain", "cervix",
    "colorectal", "kidney", "lung", "ovary", "prostate",
    "skin", "thyroid", "uterus", "breast"
}

FEMALE_ONLY = {"breast", "uterus", "ovary", "cervix"}
MALE_ONLY = {"prostate"}

# ============================================================
# 7. Safe cancer membership helper
# ============================================================
def has_dx(dx_string, cancer):
    if not dx_string:
        return False
    return cancer in {
        x.strip()
        for x in dx_string.split(";")
        if x.strip()
    }
# ============================================================
# 8. PRS loader with case-only z-scoring
# ============================================================
def load_prs(cancer, case_ids):
    if cancer not in PGS_LOOKUP:
        return None

    pgs_id = PGS_LOOKUP[cancer]
    path = os.path.join(RAW_PGS_DIR, f"{pgs_id}.raw.pgs")

    if not os.path.exists(path):
        return None

    pgs = pd.read_csv(path, sep="\t")
    
    if "sample" in pgs.columns:
        pgs = pgs.rename(columns={"sample": "original_id"})

    pgs["original_id"] = pgs["original_id"].astype(str)

    if "PGS" not in pgs.columns:
        return None

    non_cases = pgs[~pgs["original_id"].isin(case_ids)]

    mu = non_cases["PGS"].mean()
    sd = non_cases["PGS"].std(ddof=0)

    if sd == 0 or np.isnan(sd):
        return None

    pgs["prs_z"] = (pgs["PGS"] - mu) / sd

    return pgs[["original_id", "prs_z"]].rename(
        columns={"prs_z": f"{cancer}_prs"}
    )

# ============================================================
# 9. Main logistic regression scan
# ============================================================
results = []

for patient_cancer, family_cancer in product(PATIENT_CANCERS, FAMILY_CANCERS):
    sub = df.copy()

    if patient_cancer in FEMALE_ONLY:
        sub = sub[sub["reported_sex"] == "female"]
    if patient_cancer in MALE_ONLY:
        sub = sub[sub["reported_sex"] == "male"]
    if sub.empty:
        continue

    # ---- Outcomes & predictors (SAFE) ----
    sub["patient_has"] = sub["patient_dx"].apply(
        lambda x: has_dx(x, patient_cancer)
    ).astype(int)

    sub["family_has"] = sub.apply(
        lambda r: has_dx(r["maternal_dx"], family_cancer)
                  or has_dx(r["paternal_dx"], family_cancer),
        axis=1
    ).astype(int)

    if sub.loc[sub["patient_has"] == 1, "family_has"].sum() < 3:
        continue

    case_ids = set(
        sub.loc[sub["patient_has"] == 1, "original_id"]
    )

    covars = ["family_has", "sex_binary"]
    # ---- Patient PRS ----
    patient_prs = load_prs(patient_cancer, case_ids)
    if patient_prs is not None:
        sub = sub.merge(patient_prs, on="original_id", how="left")
        covars.append(f"{patient_cancer}_prs")

    # ---- Family PRS ----
    SKIN_PATIENT_SUBTYPES = {"basal_cell_carcinoma","melanoma","squamous_cell_carcinoma"}

    family_prs = load_prs(family_cancer, case_ids)
    if family_prs is not None and patient_cancer.lower() != family_cancer.lower() and not (family_cancer.lower() == "skin" and patient_cancer.lower() in SKIN_PATIENT_SUBTYPES):
        sub = sub.merge(family_prs, on="original_id", how="left")
        covars.append(f"{family_cancer}_prs")
    if patient_cancer in ("prostate","breast","ovary","uterus","cervix"):
        covars = [c for c in covars if c != "sex_binary"]
    try:
        X = sm.add_constant(sub[covars])
    except KeyError as e:
        missing = [c for c in covars if c not in sub.columns]
        print(
            f"[COVAR ERROR] patient_cancer={patient_cancer} "
            f"family_cancer={family_cancer}"
            f"missing={missing}"
        )
        continue
    y = sub["patient_has"]

    try:
        model = sm.Logit(y, X).fit(disp=False)
    except Exception as e:
        print(
          f"[MODEL ERROR] patient_cancer={patient_cancer} "
          f"family_cancer={family_cancer} "
          f"error={repr(e)}"
        )
        continue

    beta = model.params["family_has"]

    results.append({
        "patient_cancer": patient_cancer,
        "family_cancer": family_cancer,
        "odds_ratio": np.exp(beta),
        "p_value": model.pvalues["family_has"],
        "n_tested": len(sub),
        "n_intersection": int(
            ((sub["patient_has"] == 1) & (sub["family_has"] == 1)).sum()
        ),
        "covariates": ",".join(covars)
    })

# ============================================================
# 10. Output
# ============================================================
out = pd.DataFrame(results).sort_values(
    ["p_value", "odds_ratio"],
    ascending=[True, False]
)

out.to_csv(
    "patient_family_logistic_with_prs.tsv",
    sep="\t",
    index=False
)

print(f"Wrote {len(out)} rows to patient_family_logistic_with_prs.tsv")

