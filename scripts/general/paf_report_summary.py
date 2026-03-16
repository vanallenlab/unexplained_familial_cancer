#!/usr/bin/env python3
import pandas as pd
import glob
import os
import numpy as np

# ----------------------------
# Input directory
# ----------------------------
input_dir = "/Users/noah/Desktop/ufc_repository/results/paf_results/"
pattern = os.path.join(input_dir, "*adjusted*.tsv")

# ----------------------------
# Initialize lists to store values with file info
# ----------------------------
top_r2_reduced = []
bottom_r2_full = []
bottom_fraction_explained = []

top_r2_reduced_files = []
bottom_r2_full_files = []
bottom_fraction_explained_files = []

# ----------------------------
# Loop over files
# ----------------------------
for f in glob.glob(pattern):
    if "family" in f or "isolated" in f:
        continue
    
    df = pd.read_csv(f, sep="\t")
    
    # Filter for prevalence_model == 'observed'
    df = df[df["prevalence_model"] == "observed"]
    #df = df[df["added_predictor"].str.contains("SV|Tier|PGS")]
    if df.empty:
        continue
    
    # Top row: R2_reduced
    top_r2_reduced.append(df.iloc[0]["R2_reduced"])
    top_r2_reduced_files.append(f)
    
    # Bottom row: R2_full and fraction_explained
    bottom_r2_full.append(df.iloc[-1]["R2_full"])
    bottom_r2_full_files.append(f)
    
    bottom_fraction_explained.append(df.iloc[-1]["fraction_explained"])
    bottom_fraction_explained_files.append(f)

# ----------------------------
# Helper function to report min, max and files
# ----------------------------
def report_range(values, files, name):
    if not values:
        print(f"No {name} values found.")
        return
    min_val = min(values)
    max_val = max(values)
    min_file = files[values.index(min_val)]
    max_file = files[values.index(max_val)]
    print(f"{name} range: {min_val:.4f} - {max_val:.4f}")
    print(f"  min in file: {min_file}")
    print(f"  max in file: {max_file}")

# ----------------------------
# Report ranges
# ----------------------------
report_range(top_r2_reduced, top_r2_reduced_files, "Top row R2_reduced")
report_range(bottom_r2_full, bottom_r2_full_files, "Bottom row R2_full")
report_range(bottom_fraction_explained, bottom_fraction_explained_files, "Bottom row fraction_explained")

# ----------------------------
# Report average fraction_explained and file
# ----------------------------
if bottom_fraction_explained:
    avg_fraction = np.mean(bottom_fraction_explained)
    # find closest value to average
    closest_idx = np.argmin(np.abs(np.array(bottom_fraction_explained) - avg_fraction))
    avg_file = bottom_fraction_explained_files[closest_idx]
    print(f"Average fraction_explained: {avg_fraction:.4f} (closest value from file: {avg_file})")
