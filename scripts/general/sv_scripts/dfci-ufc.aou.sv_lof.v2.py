import pandas as pd
from scipy.stats import fisher_exact
from pathlib import Path

# -----------------------------
# Parameters for runs
# -----------------------------
runs = [
    {"name": "rare", "AF": 0.01, "AC": None},
    {"name": "ultra_rare", "AF": 0.001, "AC": None},
    {"name": "singleton", "AF": 0.001, "AC": 1},
]

# Input files
sv_file = "ufc.rare_svs.tsv.gz"
riaz_genes_file = "riaz_genes.list"
pheno_file = "dfci-ufc.aou.phenos.v2.tsv.gz"
exclude_file = "samples_to_exclude.jan29.list"
qc_file = "all_samples_in_ufc_sv.qc_pass.list"

# Load gene list
with open(riaz_genes_file) as f:
    riaz_genes = set(line.strip() for line in f if line.strip())

# Load phenotype
pheno = pd.read_csv(pheno_file, sep="\t", compression="gzip")
pheno["Sample"] = pheno["Sample"].astype(str)

# Exclude unwanted samples
exclude = pd.read_csv(exclude_file, header=None)[0].astype(str)
pheno = pheno[~pheno["Sample"].isin(exclude)]

# Keep QC-passed samples only
qc_pass = pd.read_csv(qc_file, header=None)[0].astype(str)
pheno = pheno[pheno["Sample"].isin(qc_pass)]

# -----------------------------
# Function to process SV table
# -----------------------------
def process_sv(df, AF, AC=None):
    df = df[df["AF"] == AF].copy()
    if AC is not None:
        df = df[df["AC"] == AC].copy()

    gene_cols = [
        "PREDICTED_INTRAGENIC_EXON_DUP",
        "PREDICTED_LOF",
        "PREDICTED_PARTIAL_EXON_DUP"
    ]
    df[gene_cols] = df[gene_cols].fillna("")
    
    # Combine gene columns
    df["SV_LOF"] = df[gene_cols].agg(",".join, axis=1)

    # Tokenize genes safely
    df["SV_LOF"] = df["SV_LOF"].apply(lambda x: list({g.strip() for g in x.split(",") if g.strip()}))

    # Keep rows with ≥1 exact match to riaz_genes
    df = df[df["SV_LOF"].apply(lambda genes: any(g in riaz_genes for g in genes))]
    df["SV_LOF"] = df["SV_LOF"].apply(lambda genes: [g for g in genes if g in riaz_genes])

    # Explode genes
    df = df.explode("SV_LOF")

    # Explode samples
    df["SAMPLES"] = df["SAMPLES"].fillna("").apply(lambda x: [s.strip() for s in x.split(",") if s.strip()])
    df = df.explode("SAMPLES")

    # Final table
    final = df[["ID", "AC", "SV_LOF", "SAMPLES"]].copy()
    final = final.rename(columns={"ID": "SV_ID", "SAMPLES": "original_id"})
    final = final.drop_duplicates(subset=["SV_LOF", "original_id"])
    
    return final

# -----------------------------
# Function to merge and do Fisher
# -----------------------------
def fisher_analysis(sv_final, pheno):
    sv_final["original_id"] = sv_final["original_id"].astype(str)
    merged = pheno.merge(
        sv_final[["original_id"]].drop_duplicates(),
        left_on="Sample",
        right_on="original_id",
        how="left"
    )

    merged["is_case"] = merged["original_dx"] != "control"
    merged["has_sv_lof"] = merged["original_id"].notna()

    n_cases = merged["is_case"].sum()
    n_controls = (~merged["is_case"]).sum()

    cases_with_sv = merged.loc[merged["is_case"], "has_sv_lof"].sum()
    controls_with_sv = merged.loc[~merged["is_case"], "has_sv_lof"].sum()

    # Fisher test
    table = [
        [cases_with_sv, n_cases - cases_with_sv],
        [controls_with_sv, n_controls - controls_with_sv]
    ]
    odds_ratio, p_value = fisher_exact(table)

    print(f"Cases: {n_cases} | Controls: {n_controls}")
    print(f"Cases w/ SV_LOF: {cases_with_sv} | Controls w/ SV_LOF: {controls_with_sv}")
    print(f"Odds Ratio = {odds_ratio:.4f} | P-value = {p_value:.4e}")
    print("-" * 50)

# -----------------------------
# Load SV once
# -----------------------------
df_sv = pd.read_csv(sv_file, sep="\t", compression="gzip")

# -----------------------------
# Run each filter
# -----------------------------
for run in runs:
    print(f"Processing {run['name']} SVs...")
    final_sv = process_sv(df_sv, AF=run["AF"], AC=run.get("AC"))
    out_file = f"dfci-ufc.aou.sv_lof.{run['name']}.tsv"
    final_sv.to_csv(out_file, sep="\t", index=False)
    print(f"Saved {out_file}")
    
    # Fisher analysis
    fisher_analysis(final_sv, pheno)
