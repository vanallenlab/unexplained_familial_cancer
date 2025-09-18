#!/usr/bin/env python3
import pandas as pd
from collections import Counter, defaultdict
import itertools
import os

INPUT = "dfci-ufc.aou.phenos.v2.tsv"   # adjust path if needed
COL = "cancers_in_FDRs"
OUTDIR = "cancer_counts_output"
os.makedirs(OUTDIR, exist_ok=True)

# Read file (works with .tsv or .tsv.gz)
df = pd.read_csv(INPUT, sep="\t", dtype=str)
series = df[COL].dropna().astype(str)

# Counters
presence_counter = Counter()        # single-cancer presence counts (anywhere in row)
exact_counter = Counter()           # exact-set counts (order-normalized)
subset_counters = defaultdict(Counter)  # subset_counters[k][combo] = inclusive counts for size k
size_counter = Counter()            # how many rows have exactly k cancers

# Process rows
for s in series:
    # split, strip, remove empty tokens, and remove duplicates within the same row
    parts = [p.strip() for p in s.split(";") if p.strip() != ""]
    if not parts:
        continue
    # Unique within row (collapse duplicates), then sort for canonical ordering
    unique_sorted = sorted(dict.fromkeys(parts))
    sset = set(unique_sorted)
    k = len(sset)
    size_counter[k] += 1

    # exact combo key (sorted; semicolon-joined)
    exact_key = ";".join(unique_sorted)
    exact_counter[exact_key] += 1

    # presence counts: each cancer appears once for this row
    for c in sset:
        presence_counter[c] += 1

    # inclusive subset counts: for each subset size >=2 up to len(sset)
    for subset_size in range(2, k + 1):
        for comb in itertools.combinations(unique_sorted, subset_size):
            subset_key = ";".join(comb)
            subset_counters[subset_size][subset_key] += 1

# Write outputs

# 1) presence (single) counts
presence_df = (
    pd.DataFrame.from_records(list(presence_counter.items()), columns=["cancer", "count"])
      .sort_values("count", ascending=False)
      .reset_index(drop=True)
)
presence_df.to_csv(os.path.join(OUTDIR, "presence_counts.tsv"), sep="\t", index=False)

# 2) exact combo counts (any size)
exact_df = (
    pd.DataFrame.from_records(list(exact_counter.items()), columns=["combination", "count"])
      .sort_values("count", ascending=False)
      .reset_index(drop=True)
)
exact_df.to_csv(os.path.join(OUTDIR, "exact_combo_counts.tsv"), sep="\t", index=False)

# 3) subset (inclusive) counts for k = 2..max_k
max_k = max(size_counter.keys()) if size_counter else 0
for k in range(2, max_k + 1):
    counter_k = subset_counters.get(k, Counter())
    df_k = (
        pd.DataFrame.from_records(list(counter_k.items()), columns=[f"combination_k{k}", "count"])
          .sort_values("count", ascending=False)
          .reset_index(drop=True)
    )
    df_k.to_csv(os.path.join(OUTDIR, f"subset_counts_k{k}.tsv"), sep="\t", index=False)

# 4) counts by size (how many rows have exactly k cancers)
size_df = (
    pd.DataFrame.from_records(list(size_counter.items()), columns=["combo_size", "rows_with_exact_size"])
      .sort_values("combo_size")
      .reset_index(drop=True)
)
size_df.to_csv(os.path.join(OUTDIR, "counts_by_exact_size.tsv"), sep="\t", index=False)

# Print quick summary
print("Wrote outputs to", OUTDIR)
print("Top 10 presence (single) counts:")
print(presence_df.head(10).to_string(index=False))
print("\nTop 10 exact combos:")
print(exact_df.head(10).to_string(index=False))
if max_k >= 2:
    print(f"\nTop 10 inclusive pairs (k=2):")
    p2 = pd.read_csv(os.path.join(OUTDIR, "subset_counts_k2.tsv"), sep="\t")
    print(p2.head(10).to_string(index=False))

