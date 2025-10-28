#!/usr/bin/env python3
import pandas as pd

# --- Input / Output ---
input_file = "ufc.cpg.variant_counts.tsv.gz"
output_file = "patient_criteria.tsv"

# --- Load data ---
df = pd.read_csv(input_file, sep="\t")

# Assume first column contains the variant criteria (like "ATM_REVEL_050_001")
criteria_col = df.columns[0]
df = df.set_index(criteria_col)

# --- Melt into long form and filter where value == 1 ---
df_long = (
    df.stack()
      .reset_index()
      .rename(columns={"level_1": "patient", 0: "value"})
)

# Keep only where the variant count == 1
df_filtered = df_long[df_long["value"] == 1][["patient", criteria_col]]

# --- Save output ---
df_filtered.to_csv(output_file, sep="\t", index=False)

print(f"✅ Saved {len(df_filtered)} rows to {output_file}")

