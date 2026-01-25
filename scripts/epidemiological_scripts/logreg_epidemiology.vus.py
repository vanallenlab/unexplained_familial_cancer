#!/usr/bin/env python3

import pandas as pd
import numpy as np
import statsmodels.api as sm
from itertools import product
import os

# ============================================================
# Configuration
# ============================================================
RAW_VUS_DIR = "raw_vus"
GENE_MATRIX_FILE = os.path.join(RAW_VUS_DIR, "gene_by_patient_binary_matrix.tsv.gz")
COSMIC_FILE = "raw_vus/cosmic_ufc.v3.tsv"

# ============================================================
# 1. Load phenotype table
# ============================================================
df = pd.read_csv("dfci-ufc.aou.phenos.v2.tsv", sep="\t")
df = df[df["original_dx"] != "control"].copy()

pheno_df = pd.read_csv(
    "dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz",
    sep="\t",
    usecols=["original_id", "inferred_sex", "cohort",'intake_qc_pop']
)
pheno_df = pheno_df[pheno_df["cohort"] == "aou"]

df["Sample"] = df["Sample"].astype(str)
pheno_df["Sample"] = pheno_df["original_id"].astype(str)

df = df.merge(pheno_df, on="Sample", how="left")
df["original_id"] = df["Sample"].astype(str)
df = df[df['intake_qc_pop'] == 'EUR']
# ============================================================
# 2. Remove pathogenic variant samples
# ============================================================
pvs = set(
    pd.read_csv("samples_to_exclude.jan29.list", header=None)[0].astype(str)
)
df = df[~df["original_id"].isin(pvs)].copy()

# ============================================================
# 3. Sex filtering + encoding
# ============================================================
df = df[df["inferred_sex"].isin(["male", "female"])].copy()
df["sex_binary"] = (df["inferred_sex"] == "female").astype(int)

# ============================================================
# 4. Normalize diagnosis strings
# ============================================================
def normalize_dx(val):
    if pd.isna(val) or isinstance(val, bool):
        return ""
    return str(val).lower()

df["patient_dx"] = df["original_dx"].apply(normalize_dx)
df["maternal_dx"] = df["maternal_family_dx"].apply(normalize_dx)
df["paternal_dx"] = df["paternal_family_dx"].apply(normalize_dx)

# ============================================================
# 5. Cancer definitions
# ============================================================
PATIENT_CANCERS = {
    "breast", "prostate", "cervix","colorectal", "thyroid",
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
# 6. Diagnosis membership helper
# ============================================================
def has_dx(dx_string, cancer):
    if not dx_string:
        return False
    return cancer in {
        x.strip() for x in dx_string.split(";") if x.strip()
    }

# ============================================================
# 7. Load gene-by-patient binary matrix
# ============================================================
gene_mat = pd.read_csv(GENE_MATRIX_FILE, sep="\t", compression="gzip", index_col=0)
gene_mat = gene_mat.reset_index().rename(columns={"index": "Sample"})
# Ensure 'Sample' is string
gene_mat["Sample"] = gene_mat["Sample"].astype(str)
# ============================================================
# 8. Load COSMIC gene–cancer mapping
# ============================================================
cosmic = pd.read_csv(
    COSMIC_FILE,
    sep="\t",
    header=None,
    names=["gene", "cancers"]
)

cosmic["cancers"] = cosmic["cancers"].str.split(",")
# Build lookup: cancer -> set of genes
cancer_to_genes = {}
for _, row in cosmic.iterrows():
    for cancer in row["cancers"]:
        cancer_to_genes.setdefault(cancer, set()).add(row["gene"])

# ============================================================
# 9. Main logistic regression scan
# ============================================================
results = []

for patient_cancer, family_cancer in product(PATIENT_CANCERS, FAMILY_CANCERS):
    sub = df.copy()

    #if patient_cancer != "breast" or family_cancer != "breast":
    #    continue
    # -----------------------
    # Sex restrictions
    # -----------------------
    if patient_cancer in FEMALE_ONLY:
        sub = sub[sub["inferred_sex"] == "female"]
    if patient_cancer in MALE_ONLY:
        sub = sub[sub["inferred_sex"] == "male"]
    if sub.empty:
        continue

    # -----------------------
    # Outcomes
    # -----------------------
    sub["patient_has"] = sub["patient_dx"].apply(
        lambda x: has_dx(x, patient_cancer)
    ).astype(int)

    sub["family_has"] = sub.apply(
        lambda r: has_dx(r["maternal_dx"], family_cancer)
                  or has_dx(r["paternal_dx"], family_cancer),
        axis=1
    ).astype(int)

    # Require minimal signal
    if sub.loc[sub["patient_has"] == 1, "family_has"].sum() < 2:
        continue

    covars = ["family_has", "sex_binary"]

    # -----------------------
    # Gene-based covariates
    # -----------------------
    relevant_genes = cancer_to_genes.get(patient_cancer, set())
    #print(patient_cancer, relevant_genes)
    #exit(0)
    if relevant_genes:
        gene_cols = [g for g in relevant_genes if g in gene_mat.columns]
        if gene_cols:
            #print(gene_cols)
            sub = sub.merge(
                gene_mat[["Sample"] + gene_cols],
                on="Sample",
                how="left"
            )
            # Keep only genes with ≥1 carrier among cases
            case_mask = sub["patient_has"] == 1
            active_genes = [
                g for g in gene_cols
                if sub.loc[case_mask, g].sum() > 0
            ]

            covars.extend(active_genes)
            #print(covars)
            #exit(0)
    # Remove sex if redundant
    if patient_cancer in FEMALE_ONLY | MALE_ONLY:
        covars = [c for c in covars if c != "sex_binary"]

    # -----------------------
    # Model
    # -----------------------
    try:
        X = sm.add_constant(sub[covars])
        y = sub["patient_has"]
        model = sm.Logit(y, X).fit(disp=False)
    except Exception as e:
        print(
            f"[MODEL ERROR] patient={patient_cancer} family={family_cancer} "
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
        "gene_covariates": ",".join(
            [c for c in covars if c not in {"family_has", "sex_binary"}]
        )
    })
    #print(results)
# ============================================================
# 10. Output
# ============================================================
out = pd.DataFrame(results).sort_values(
    ["p_value", "odds_ratio"],
    ascending=[True, False]
)

out.to_csv(
    "patient_family_logistic_with_vus.tsv",
    sep="\t",
    index=False
)

print(f"Wrote {len(out)} rows to patient_family_logistic_with_vus.tsv")

