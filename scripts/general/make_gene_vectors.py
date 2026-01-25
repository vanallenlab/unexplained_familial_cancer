#!/usr/bin/env python3

import pandas as pd

# -------------------------
# Input / Output
# -------------------------
infile = "step_10_output.jan20_2025.tsv.gz"
outfile = "gene_by_patient_binary_matrix.tsv.gz"

# -------------------------
# Load data
# -------------------------
df = pd.read_csv(infile, sep="\t", compression="gzip")

# First column is gene_impact
gene_impact_col = df.columns[0]
patient_cols = df.columns[1:]

# 🚫 remove Tier0 rows
df = df[~df[gene_impact_col].str.contains("_Tier0_")]

# -------------------------
# Extract gene name
# Example: APC_Tier3_0001 -> APC
# -------------------------
df["gene"] = df[gene_impact_col].str.split("_").str[0]

# -------------------------
# Collapse Tier0–5 within gene
# For each gene, take max across all its Tier rows
# -------------------------
gene_patient_matrix = (
    df.groupby("gene")[patient_cols]
      .max()
      .T
)

# -------------------------
# Ensure binary (0/1)
# -------------------------
gene_patient_matrix = (gene_patient_matrix > 0).astype(int)

# -------------------------
# Save
# -------------------------
gene_patient_matrix.to_csv(outfile, sep="\t", compression="gzip")

# -------------------------
# Sanity checks
# -------------------------
print("Matrix shape (patients x genes):", gene_patient_matrix.shape)
print("Value counts:")
print(pd.Series(gene_patient_matrix.values.ravel()).value_counts())

