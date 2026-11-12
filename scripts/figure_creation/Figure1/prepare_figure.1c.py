import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


METADATA_FILES = {
    "CDRv7": "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/CDRv7.metadata",
    "UFC": "gs://dataproc-staging-wb-cordial-diamond-9893/UFC_REFERENCE_FILES/CDRv7.ufc_original.metadata",
}

OUTPUT_DIR = Path("familial_cancer_comparison")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CANCERS = {
    "Breast": {
        "cancer_col": "breast_cancer",
        "fdr_col": "fdrs_w_breast_cancer",
        "sex": "Female",
    },
    "Colon": {
        "cancer_col": "colon_cancer",
        "fdr_col": "fdrs_w_colon_cancer",
        "sex": None,
    },
    "Lung": {
        "cancer_col": "lung_cancer",
        "fdr_col": "fdrs_w_lung_cancer",
        "sex": None,
    },
    "Prostate": {
        "cancer_col": "prostate_cancer",
        "fdr_col": "fdrs_w_prostate_cancer",
        "sex": "Male",
    },
}


def read_metadata(gs_path):
    """Download a metadata file from GCS and read it as a TSV."""
    with tempfile.NamedTemporaryFile(suffix=".tsv") as tmp:
        subprocess.run(
            ["gsutil", "cp", gs_path, tmp.name],
            check=True,
        )
        return pd.read_csv(
            tmp.name,
            sep="\t",
            dtype={"person_id": str},
        )


def get_cancer_counts(df, cohort_name):
    """Count cancer phenotypes and concordant family history by cancer."""
    results = []

    for cancer, config in CANCERS.items():
        cancer_col = config["cancer_col"]
        fdr_col = config["fdr_col"]

        cancer_df = df.copy()

        if config["sex"] is not None:
            cancer_df = cancer_df[
                cancer_df["sex_at_birth"] == config["sex"]
            ].copy()

        # Treat missing phenotype and family-history values as zero.
        cancer_df[cancer_col] = cancer_df[cancer_col].fillna(0)
        cancer_df[fdr_col] = cancer_df[fdr_col].fillna(0)

        affected = cancer_df[cancer_df[cancer_col] == 1]

        n_cancer = len(affected)
        n_familial = int((affected[fdr_col] >= 1).sum())
        n_without_familial_history = n_cancer - n_familial

        results.append({
            "cancer": cancer,
            "cohort": cohort_name,
            "sex_restriction": config["sex"] or "None",
            "n_cancer_affected": n_cancer,
            "n_with_concordant_fdr": n_familial,
            "n_without_concordant_fdr": n_without_familial_history,
            "familial_cancer_rate": (
                n_familial / n_cancer if n_cancer > 0 else float("nan")
            ),
        })

    return results


def compare_cohorts(ufc_df, non_ufc_df):
    """Compare concordant family-history rates among cancer-affected people."""
    results = []

    for cancer, config in CANCERS.items():
        cancer_col = config["cancer_col"]
        fdr_col = config["fdr_col"]

        cohort_counts = {}

        for cohort_name, df in [
            ("UFC", ufc_df),
            ("Non-UFC CDRv7", non_ufc_df),
        ]:
            cancer_df = df.copy()

            if config["sex"] is not None:
                cancer_df = cancer_df[
                    cancer_df["sex_at_birth"] == config["sex"]
                ].copy()

            cancer_df[cancer_col] = cancer_df[cancer_col].fillna(0)
            cancer_df[fdr_col] = cancer_df[fdr_col].fillna(0)

            affected = cancer_df[cancer_df[cancer_col] == 1]

            n_cancer = len(affected)
            n_familial = int((affected[fdr_col] >= 1).sum())
            n_no_familial = n_cancer - n_familial

            cohort_counts[cohort_name] = {
                "n_cancer": n_cancer,
                "n_familial": n_familial,
                "n_no_familial": n_no_familial,
                "rate": (
                    n_familial / n_cancer
                    if n_cancer > 0
                    else float("nan")
                ),
            }

        ufc = cohort_counts["UFC"]
        non_ufc = cohort_counts["Non-UFC CDRv7"]

        # Rows: UFC, non-UFC CDRv7
        # Columns: >=1 concordant FDR, no concordant FDR
        table = [
            [ufc["n_familial"], ufc["n_no_familial"]],
            [non_ufc["n_familial"], non_ufc["n_no_familial"]],
        ]

        odds_ratio, p_value = fisher_exact(
            table,
            alternative="two-sided",
        )

        results.append({
            "cancer": cancer,
            "sex_restriction": config["sex"] or "None",
            "UFC_n_cancer_affected": ufc["n_cancer"],
            "UFC_n_with_concordant_fdr": ufc["n_familial"],
            "UFC_n_without_concordant_fdr": ufc["n_no_familial"],
            "UFC_familial_cancer_rate": ufc["rate"],
            "Non_UFC_n_cancer_affected": non_ufc["n_cancer"],
            "Non_UFC_n_with_concordant_fdr": non_ufc["n_familial"],
            "Non_UFC_n_without_concordant_fdr": non_ufc["n_no_familial"],
            "Non_UFC_familial_cancer_rate": non_ufc["rate"],
            "Fisher_OR_UFC_vs_non_UFC": odds_ratio,
            "Fisher_p_value": p_value,
        })

    return results


# Load both metadata files.
cdrv7_df = pd.read_csv("CDRv7.metadata",sep='\t',index_col=False, dtype = {"person_id": str}) 
ufc_df = pd.read_csv("CDRv7.ufc_original.metadata",sep='\t',index_col=False, dtype = {"person_id": str})

# The UFC cohort is contained within CDRv7. Remove UFC participants
# from CDRv7 to create an independent comparison population.
ufc_ids = set(ufc_df["person_id"].dropna())

cdrv7_df = cdrv7_df[
    ~cdrv7_df["person_id"].isin(ufc_ids)
].copy()

# Summarize cancer counts and family-history rates in all three groups.
summary_results = []
summary_results.extend(get_cancer_counts(cdrv7_df, "CDRv7 overall"))
summary_results.extend(get_cancer_counts(ufc_df, "UFC"))

summary_df = pd.DataFrame(summary_results)

# Compare UFC against the non-overlapping remainder of CDRv7.
comparison_df = pd.DataFrame(
    compare_cohorts(ufc_df, cdrv7_df)
)

# Save files for downstream graphing.
summary_path = OUTPUT_DIR / "familial_cancer_counts.csv"
comparison_path = OUTPUT_DIR / "familial_cancer_comparison.csv"

summary_df.to_csv(summary_path, index=False)
comparison_df.to_csv(comparison_path, index=False)

# Print results.
print("\nFAMILIAL CANCER COUNTS AND RATES")
print(summary_df.to_string(index=False))

print("\nUFC VS NON-UFC CDRv7: FISHER'S EXACT TEST")
print(comparison_df.to_string(index=False))

print(f"\nSaved: {summary_path}")
print(f"Saved: {comparison_path}")
