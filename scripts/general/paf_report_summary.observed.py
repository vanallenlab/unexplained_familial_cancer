#!/usr/bin/env python3
import pandas as pd
import glob
import os
import numpy as np

input_dir = "/Users/noah/Desktop/ufc_repository/results/paf_results/"
pattern = os.path.join(input_dir, "*adjusted*.tsv")

# ----------------------------
# Storage
# ----------------------------
baseline_r2 = []

selected_r2_full = []  # PGS or Tier depending on cancer

delta_prev_model = []  # PGS - previous row

novel_variance = []  # final - PGS/Tier

# cancers requiring Tier instead of PGS
tier_cancers = ["sarcoma", "neuroendocrine"]

# cancers for "novel variance"
novel_cancers = ["neuroendocrine", "thyroid", "kidney", "breast", "bladder", "cervix"]

# ----------------------------
# Helper
# ----------------------------
def get_cancer_type(filepath):
    fname = os.path.basename(filepath).lower()
    for c in tier_cancers + novel_cancers:
        if c in fname:
            return c
    return fname  # fallback

def report(values, name):
    if not values:
        print(f"No values for {name}")
        return
    print(f"{name}:")
    print(f"  range: {min(values):.4f} - {max(values):.4f}")
    print(f"  median: {np.median(values):.4f}")
    print()

def report_mean(values, name):
    if not values:
        print(f"No values for {name}")
        return
    print(f"{name}:")
    print(f"  range: {min(values):.4f} - {max(values):.4f}")
    print(f"  mean: {np.mean(values):.4f}")
    print()

# ----------------------------
# Main loop
# ----------------------------
for f in glob.glob(pattern):
    if "family" in f or "isolated" in f:
        continue
    
    df = pd.read_csv(f, sep="\t")
    df = df[df["prevalence_model"] == "observed"].reset_index(drop=True)
    
    if df.empty:
        continue
    
    cancer = get_cancer_type(f)
    
    # ----------------------------
    # Baseline (top row)
    # ----------------------------
    baseline = df.iloc[0]["R2_reduced"]
    baseline_r2.append(baseline)
    
    # ----------------------------
    # Select final row: PGS or Tier
    # ----------------------------
    if any(c in f.lower() for c in tier_cancers):
        subset = df[df["added_predictor"].str.contains("Tier", na=False)]
    else:
        subset = df[df["added_predictor"].str.contains("PGS", na=False)]
    
    if subset.empty:
        continue
    
    idx = subset.index[-1]  # take last matching row
    selected_r2 = df.loc[idx, "R2_full"]
    selected_r2_full.append(selected_r2 - baseline)
    
    # ----------------------------
    # Delta vs previous (PGS only, excluding tier cancers)
    # ----------------------------
    if not any(c in f.lower() for c in tier_cancers):
        if idx > 0 and not pd.isna(df.loc[idx - 1, "R2_full"]):
            prev = df.loc[idx - 1, "R2_full"]
        else:
            prev = baseline
        
        delta_prev_model.append(selected_r2 - prev)
    
    # ----------------------------
    # Novel variance
    # ----------------------------
    if any(c in f.lower() for c in novel_cancers):
        final_r2 = df.iloc[-1]["R2_full"]
        print(f.lower(),final_r2 - selected_r2)
        novel_variance.append(final_r2 - selected_r2)

# ----------------------------
# Reporting
# ----------------------------
report(baseline_r2, "Baseline R2_reduced")
report(selected_r2_full, "Added variance explained by PGS and ")
report_mean(delta_prev_model, "Delta vs previous model (PGS only)")
report(novel_variance, "Added variance from novel risk factors")

# ----------------------------
# Build summary table
# ----------------------------
summary_rows = []

def add_summary(values, name):
    if not values:
        summary_rows.append({
            "metric": name,
            "min": np.nan,
            "max": np.nan,
            "median": np.nan
        })
    else:
        summary_rows.append({
            "metric": name,
            "min": np.min(values),
            "max": np.max(values),
            "median": np.median(values)
        })

# add_summary(baseline_r2, "Baseline R2_reduced")
# add_summary(selected_r2_full, "Selected R2_full (PGS/Tier)")
# add_summary(delta_prev_model, "Delta vs previous model (PGS only)")
# add_summary(novel_variance, "Novel variance explained")

# summary_df = pd.DataFrame(summary_rows)

# # ----------------------------
# # Print nicely
# # ----------------------------
# pd.set_option("display.float_format", "{:.4f}".format)
# print("\nSummary Table:\n")
# print(summary_df)

# # ----------------------------
# # Save to file (optional)
# # ----------------------------
# out_path = os.path.join(input_dir, "summary_statistics.tsv")
# summary_df.to_csv(out_path, sep="\t", index=False)
# print(f"\nSaved summary table to: {out_path}")