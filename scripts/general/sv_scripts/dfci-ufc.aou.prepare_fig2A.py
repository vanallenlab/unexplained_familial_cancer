import pandas as pd

# Input files
sv_file = "dfci-ufc.aou.sv_lof.rare.tsv"
phenos_file = "dfci-ufc.aou.phenos.v2.tsv.gz"

# Phenotype list
CANCERS = [
    "Basal_Cell_Carcinoma","Breast","Prostate","Squamous_Cell_Carcinoma","Melanoma",
    "Lung","Colorectal","Thyroid","Hematologic","Non-Hodgkin","Uterus","Bladder",
    "Sarcoma","Cervix","Ovary","Kidney","Brain","Neuroendocrine","control"
]

# Load SV LOF data
sv_df = pd.read_csv(sv_file, sep="\t")
# Make sure original_id exists
if 'original_id' not in sv_df.columns:
    raise ValueError("SV file must have 'original_id' column")

# Load phenotype data
phenos_df = pd.read_csv(phenos_file, sep="\t", compression='gzip')

# Merge on sample ID
merged = sv_df.merge(phenos_df, left_on='original_id', right_on='Sample', how='inner')

# Tokenize original_dx by ';' and keep only tokens in CANCERS
def filter_dx(dx):
    if pd.isna(dx):
        return []
    return [token for token in dx.split(';') if token in CANCERS]

# Create a new column with filtered tokens
merged['filtered_dx'] = merged['original_dx'].apply(filter_dx)

# Explode so each token becomes a separate row
merged_expanded = merged.explode('filtered_dx')

# Drop rows where no tokens were retained
merged_expanded = merged_expanded[merged_expanded['filtered_dx'].notna()]

# Create binary gene × phenotype matrix
matrix_binary = pd.crosstab(merged_expanded['SV_LOF'], merged_expanded['filtered_dx'])
matrix_binary = (matrix_binary > 0).astype(int)

# Save output
matrix_binary.to_csv("sv_gene_by_phenotype_binary.fig2a_prep.tsv", sep="\t")

print("\nGene × phenotype matrix (binary):")
print(matrix_binary.head())
