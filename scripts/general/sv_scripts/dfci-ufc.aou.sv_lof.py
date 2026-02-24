import pandas as pd

# -----------------------------
# Load inputs
# -----------------------------
df = pd.read_csv("ufc.rare_svs.tsv.gz", sep="\t", compression="gzip")

# Load riaz gene list
with open("riaz_genes.list") as f:
    riaz_genes = set(line.strip() for line in f if line.strip())

# -----------------------------
# Filter AF == 0.001
# -----------------------------
df = df[(df["AF"] == 0.001)].copy()

# -----------------------------
# Combine SV gene columns
# -----------------------------
gene_cols = [
    "PREDICTED_INTRAGENIC_EXON_DUP",
    "PREDICTED_LOF",
    "PREDICTED_PARTIAL_EXON_DUP"
]
# Fill NA with empty string
df[gene_cols] = df[gene_cols].fillna("")

# Combine into one comma-delimited column
df["SV_LOF"] = (
    df["PREDICTED_INTRAGENIC_EXON_DUP"] + "," +
    df["PREDICTED_LOF"] + "," +
    df["PREDICTED_PARTIAL_EXON_DUP"]
)
# -----------------------------
# Tokenize genes safely
# -----------------------------
def tokenize_genes(x):
    if not x:
        return []
    return list({
        g.strip() for g in x.split(",")
        if g.strip()
    })

df["SV_LOF"] = df["SV_LOF"].apply(tokenize_genes)
# Keep rows with ≥1 exact match to riaz_genes
df = df[df["SV_LOF"].apply(lambda genes: any(g in riaz_genes for g in genes))]
df.to_csv("tmp.tsv",sep='\t')
# Keep only matching genes
df["SV_LOF"] = df["SV_LOF"].apply(
    lambda genes: [g for g in genes if g in riaz_genes]
)

# -----------------------------
# Explode genes
# -----------------------------
df = df.explode("SV_LOF")
#df.to_csv("tmp.tsv",sep='\t')
# -----------------------------
# Explode samples
# -----------------------------
df["SAMPLES"] = df["SAMPLES"].fillna("").apply(
    lambda x: [s.strip() for s in x.split(",") if s.strip()]
)
df = df.explode("SAMPLES")
df.to_csv("tmp.tsv",sep='\t')
# -----------------------------
# Final table
# -----------------------------
final = df[["ID", "AC", "SV_LOF","SAMPLES"]].copy()
final = final.rename(columns={"ID": "SV_ID"})
final = final.rename(columns={"SAMPLES": "original_id"})
final = final.drop_duplicates(subset=["SV_LOF", "original_id"])
# -----------------------------
# Write output
# -----------------------------
final.to_csv("dfci-ufc.aou.sv_lof.tsv", sep="\t", index=False)
print(f"Unique SV_ID: {final['SV_ID'].nunique()} | Unique SV_LOF: {final['SV_LOF'].nunique()} | Unique original_id: {final['original_id'].nunique()}")
print(final.groupby("SV_LOF")["SV_ID"].nunique().loc[lambda x: x > 1])


#### Part 2
### Final SV Stat stuff
import pandas as pd
from scipy.stats import fisher_exact

# -----------------------------
# Read phenotype file
# -----------------------------
pheno = pd.read_csv("dfci-ufc.aou.phenos.v2.tsv.gz", sep="\t", compression="gzip")
pheno["Sample"] = pheno["Sample"].astype(str)

# -----------------------------
# Exclude unwanted samples
# -----------------------------
exclude = pd.read_csv("samples_to_exclude.jan29.list", header=None)[0].astype(str)
pheno = pheno[~pheno["Sample"].isin(exclude)]

# -----------------------------
# Keep QC-passed samples only
# -----------------------------
qc_pass = pd.read_csv("all_samples_in_ufc_sv.qc_pass.list", header=None)[0].astype(str)
pheno = pheno[pheno["Sample"].isin(qc_pass)]

# -----------------------------
# Read SV LoF file
# -----------------------------
sv = pd.read_csv("dfci-ufc.aou.sv_lof.tsv", sep="\t")
sv["original_id"] = sv["original_id"].astype(str)

# -----------------------------
# Merge (Sample -> original_id)
# -----------------------------
merged = pheno.merge(
    sv[["original_id"]].drop_duplicates(),
    left_on="Sample",
    right_on="original_id",
    how="left"
)

# -----------------------------
# Define case/control
# -----------------------------
merged["is_case"] = merged["original_dx"] != "control"
merged["has_sv_lof"] = merged["original_id"].notna()

# -----------------------------
# Counts
# -----------------------------
n_cases = merged["is_case"].sum()
n_controls = (~merged["is_case"]).sum()

cases_with_sv = merged.loc[merged["is_case"], "has_sv_lof"].sum()
controls_with_sv = merged.loc[~merged["is_case"], "has_sv_lof"].sum()

# Percentages
pct_cases = cases_with_sv / n_cases * 100 if n_cases > 0 else 0
pct_controls = controls_with_sv / n_controls * 100 if n_controls > 0 else 0

# -----------------------------
# Fisher's Exact Test
# -----------------------------
table = [
    [cases_with_sv, n_cases - cases_with_sv],
    [controls_with_sv, n_controls - controls_with_sv]
]

odds_ratio, p_value = fisher_exact(table)

# -----------------------------
# Print results
# -----------------------------
print(f"Cases: {n_cases}")
print(f"Controls: {n_controls}")
print()
print(f"Cases w/ SV_LOF: {cases_with_sv} ({pct_cases:.2f}%)")
print(f"Controls w/ SV_LOF: {controls_with_sv} ({pct_controls:.2f}%)")
print()
print("Fisher's Exact Test:")
print(f"Odds Ratio = {odds_ratio:.4f}")
print(f"P-value = {p_value:.4e}")
