import pandas as pd
from scipy.stats import fisher_exact
from pathlib import Path

# -----------------------------
# Inputs
# -----------------------------
sv_files = {
    "rare": "dfci-ufc.aou.sv_lof.rare.tsv",
    "ultra_rare": "dfci-ufc.aou.sv_lof.ultra_rare.tsv",
    "singleton": "dfci-ufc.aou.sv_lof.singleton.tsv"
}

metadata_dir = Path("metadatas/")
metadata_files = list(metadata_dir.glob("*.metadata"))

qc_file = "all_samples_in_ufc_sv.qc_pass.list"
qc_pass = pd.read_csv(qc_file, header=None)[0].astype(str)

# -----------------------------
# Function to load SV table and mark LoF
# -----------------------------
def load_sv_lof(file_path):
    df = pd.read_csv(file_path, sep="\t")
    df["original_id"] = df["original_id"].astype(str)
    # Mark whether individual has ≥1 SV_LOF
    df["has_sv_lof"] = df["SV_LOF"].notna()
    return df[["original_id", "has_sv_lof"]].drop_duplicates(subset=["original_id"])

# -----------------------------
# Function to compute Fisher stats
# -----------------------------
import numpy as np

def fisher_stats(df, cohort_name, sv_type):
    n_cases = df["is_case"].sum()
    n_controls = (~df["is_case"]).sum()
    a = df.loc[df["is_case"], "has_sv_lof"].sum()
    c = df.loc[~df["is_case"], "has_sv_lof"].sum()
    b = n_cases - a
    d = n_controls - c

    table = [[a, b], [c, d]]
    odds_ratio, p_value = fisher_exact(table)

    # Add 0.5 continuity correction if any cell is zero
    if 0 in [a, b, c, d]:
        a += 0.5
        b += 0.5
        c += 0.5
        d += 0.5

    se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    ci_lower = np.exp(np.log(odds_ratio) - 1.96 * se)
    ci_upper = np.exp(np.log(odds_ratio) + 1.96 * se)

    return {
        "cohort": cohort_name,
        "sv_type": sv_type,
        "n_cases": n_cases,
        "n_controls": n_controls,
        "cases_with_sv": a,
        "controls_with_sv": c,
        "odds_ratio": odds_ratio,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "p_value": p_value
    }
"""
def fisher_stats(df, cohort_name, sv_type):
    n_cases = (df["is_case"]).sum()
    n_controls = (~df["is_case"]).sum()
    cases_with_sv = df.loc[df["is_case"], "has_sv_lof"].sum()
    controls_with_sv = df.loc[~df["is_case"], "has_sv_lof"].sum()

    # Build table and run Fisher
    table = [
        [cases_with_sv, n_cases - cases_with_sv],
        [controls_with_sv, n_controls - controls_with_sv]
    ]
    odds_ratio, p_value = fisher_exact(table)

    # Return a dict for summary
    return {
        "cohort": cohort_name,
        "sv_type": sv_type,
        "n_cases": n_cases,
        "n_controls": n_controls,
        "odds_ratio": odds_ratio,
        "p_value": p_value
    }
"""
# -----------------------------
# Load metadata and loop
# -----------------------------
summary_stats = []

for meta_file in metadata_files:
    cohort_name = meta_file.stem
    print(f"Processing cohort {cohort_name}...")

    # Load cohort metadata
    pheno = pd.read_csv(meta_file, sep="\t")
    pheno["original_id"] = pheno["Sample"].astype(str)
    pheno = pheno[pheno["original_id"].isin(qc_pass)]
    pheno["is_case"] = pheno["original_dx"] != "control"

    # Loop over all 3 SV frequencies
    for sv_type, sv_file in sv_files.items():
        sv_lof = load_sv_lof(sv_file)

        # Merge with cohort metadata
        merged = pheno.merge(
            sv_lof,
            on="original_id",
            how="left"
        )
        # If null, mark as False
        merged["has_sv_lof"] = merged["has_sv_lof"].fillna(False)

        # Compute Fisher stats
        stats = fisher_stats(merged, cohort_name, sv_type)
        summary_stats.append(stats)

# -----------------------------
# Create summary table
# -----------------------------
summary_df = pd.DataFrame(summary_stats)
summary_df.to_csv("sv_lof_fisher_summary.tsv", sep="\t", index=False)
print("Saved sv_lof_fisher_summary.tsv")
