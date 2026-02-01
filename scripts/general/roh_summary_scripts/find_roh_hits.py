#!/usr/bin/env python3

import pandas as pd
import os
import re
import glob

# ============================================================
# Configuration
# ============================================================
ROH_DIR = "/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/results"
PVALUE_THRESHOLD = 1e-5
OUTFILE = "roh_significant_hits.tsv"

# ============================================================
# Filename parser
# Example:
# basal_cell.sw_100000.roh_10000.short.tsv.gz
# ============================================================
FILENAME_RE = re.compile(
    r"(?P<cancer>.+?)\.sw_(?P<sw>\d+)\.roh_(?P<roh>\d+)\.short\.tsv\.gz$"
)

hits = []

# ============================================================
# Main loop
# ============================================================
for path in glob.glob(os.path.join(ROH_DIR, "*.tsv.gz")):
    fname = os.path.basename(path)
    m = FILENAME_RE.match(fname)
    if m is None:
        print(f"[SKIP] Could not parse filename: {fname}")
        continue

    cancer = m.group("cancer")
    sw_size = int(m.group("sw"))
    roh_size = int(m.group("roh"))

    try:
        df = pd.read_csv(path, sep="\t")
        # Strip whitespace and convert
        df = df[df['mean_case_value'] != 0]
    except Exception as e:
        print(f"[ERROR] Failed to read {fname}: {e}")
        continue

    if "pvalue_roh_score" not in df.columns:
        print(f"[SKIP] Missing pvalue_roh_score in {fname}")
        continue

    sig = df[df["pvalue_roh_score"] < PVALUE_THRESHOLD].copy()
    if sig.empty:
        continue

    sig["cancer_type"] = cancer
    sig["sliding_window"] = sw_size
    sig["roh_size"] = roh_size
    sig["source_file"] = fname
    

    hits.append(sig)

# ============================================================
# Output
# ============================================================
if not hits:
    print("No significant ROH hits found.")
    exit(0)

out = pd.concat(hits, ignore_index=True)

# Ensure pvalue column is numeric first
out['pvalue_roh_score'] = pd.to_numeric(out['pvalue_roh_score'], errors='coerce')

# Drop any rows where conversion failed (optional)
out = out.dropna(subset=['pvalue_roh_score'])

# Sort ascending: smallest p-values first
out = out.sort_values(by='pvalue_roh_score', ascending=True)

# Optional: reset index
out = out.reset_index(drop=True)

out.to_csv(OUTFILE, sep="\t", index=False)
print(f"Wrote {len(out)} significant ROH hits to {OUTFILE}")
