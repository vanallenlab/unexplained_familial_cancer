#!/usr/bin/env python3

import numpy as np
import pandas as pd
from pathlib import Path

# ------------------------------------------------------------
# Hardcoded inputs (as requested)
# ------------------------------------------------------------
METADATA = "metadatas/kidney.metadata"
ROH_TSV  = "ufc_roh_chr22.tsv"                 # new ROH file
HAPLOTYPE = "chr22_32100001_32550000"          # new haplotype format
OUT_TSV = "kidney.chr22_32100001_32550000.roh_signal.tsv"

# ------------------------------------------------------------
# Helper: parse haplotype string
# Format: chr22_32100001_32550000
# ------------------------------------------------------------
def parse_haplotype(h):
    parts = h.split("_")
    if len(parts) != 3:
        raise ValueError(f"Bad haplotype format: {h} (expected chr22_32100001_32550000)")
    chrom = parts[0]
    start = int(parts[1])
    end = int(parts[2])
    if start > end:
        raise ValueError(f"Haplotype start > end: {start} > {end}")
    return chrom, start, end


def main():

    # -----------------------------
    # 1) Load metadata
    # -----------------------------
    meta = pd.read_csv(METADATA, sep="\t", dtype=str)

    if "original_id" not in meta.columns or "original_dx" not in meta.columns:
        raise SystemExit("metadata must contain columns: original_id, original_dx")

    meta["original_id"] = meta["original_id"].astype(str).str.strip()
    meta["original_dx"] = meta["original_dx"].fillna("").astype(str).str.strip()

    cases = set(meta.loc[meta["original_dx"] != "control", "original_id"])
    ctrls = set(meta.loc[meta["original_dx"] == "control", "original_id"])
    all_samples = set(meta["original_id"])

    n_cases = len(cases)
    n_ctrls = len(ctrls)

    if n_cases == 0 or n_ctrls == 0:
        raise SystemExit(f"Need both cases and controls. Found n_cases={n_cases}, n_ctrls={n_ctrls}")

    # -----------------------------
    # 2) Parse haplotype region
    # -----------------------------
    hap_chr, hap_start, hap_end = parse_haplotype(HAPLOTYPE)
    hap_len = hap_end - hap_start + 1

    print(f"Haplotype window: {hap_chr}:{hap_start}-{hap_end} (len={hap_len:,})")
    print(f"Cases={n_cases:,}  Controls={n_ctrls:,}")

    # -----------------------------
    # 3) Load ROH file
    # New format: chr, start, end, length, original_id
    # -----------------------------
    roh = pd.read_csv(
        ROH_TSV,
        sep="\t",
        names = ['chr','start','end','original_id','length'],
        dtype={"chr": str, "start": int, "end": int, "length": int, "original_id": str}
    )
    required_cols = {"chr", "start", "end", "length", "original_id"}
    missing = required_cols - set(roh.columns)
    if missing:
        raise SystemExit(f"ROH file missing columns: {sorted(missing)}")

    roh["original_id"] = roh["original_id"].astype(str).str.strip()
    print("overlap:", len(set(all_samples) & set(roh["original_id"].astype(str))))

    # Only keep ROHs from kidney metadata samples
    roh = roh[roh["original_id"].isin(all_samples)].copy()
    print(len(roh),"length")
    # Only keep ROHs on the correct chromosome
    roh = roh[roh["chr"].astype(str) == hap_chr].copy()
    print(len(roh))
    # Only keep ROHs overlapping the haplotype window
    roh = roh[(roh["start"] <= hap_end) & (roh["end"] >= hap_start)].copy()

    if roh.empty:
        raise SystemExit("No ROHs overlap this haplotype window after filtering.")

    print(f"ROHs overlapping window: {roh.shape[0]:,}")

    # -----------------------------
    # 4) Build per-base signal (fast)
    #
    # We do NOT iterate base-by-base for each ROH.
    # Instead:
    #   - use a difference array
    #   - add +inc at start_idx
    #   - add -inc at end_idx+1
    # then cumulative sum gives per-base signal
    # -----------------------------
    case_inc = 1.0 / n_cases
    ctrl_inc = 1.0 / n_ctrls

    case_diff = np.zeros(hap_len + 1, dtype=float)
    ctrl_diff = np.zeros(hap_len + 1, dtype=float)

    # classify each ROH row as case vs control
    # (if it’s neither, ignore it)
    for _, r in roh.iterrows():
        sid = r["original_id"]

        # clip to hap boundaries
        s = max(int(r["start"]), hap_start)
        e = min(int(r["end"]), hap_end)
        if s > e:
            continue

        s_idx = s - hap_start
        e_idx = e - hap_start

        if sid in cases:
            case_diff[s_idx] += case_inc
            case_diff[e_idx + 1] -= case_inc
        elif sid in ctrls:
            ctrl_diff[s_idx] += ctrl_inc
            ctrl_diff[e_idx + 1] -= ctrl_inc

    case_signal = np.cumsum(case_diff[:-1])
    ctrl_signal = np.cumsum(ctrl_diff[:-1])

    # -----------------------------
    # 5) Output table for plotting
    # -----------------------------
    out = pd.DataFrame({
        "chr": hap_chr,
        "pos": np.arange(hap_start, hap_end + 1, dtype=int),
        "case_roh_fraction": np.round(case_signal, 4),
        "control_roh_fraction": np.round(ctrl_signal, 4),
    })

    out.to_csv(OUT_TSV, sep="\t", index=False)
    print("Wrote:", OUT_TSV)


if __name__ == "__main__":
    main()

