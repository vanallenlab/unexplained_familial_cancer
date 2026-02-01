#!/usr/bin/env python3

import argparse
import pandas as pd
from sklearn.metrics import roc_auc_score


def clean_haplotype_name(h):
    # match your original replacement behavior
    return h.replace(":", "_").replace("-", "_")


def main():
    parser = argparse.ArgumentParser(
        description="Compute AUROC of percent_homozygosed vs cancer status for a given haplotype."
    )
    parser.add_argument("--metadata", required=True, help="Path to metadata TSV")
    parser.add_argument("--hap-matrix", required=True, help="Path to haplotype matrix TSV")
    parser.add_argument("--haplotype", required=True, help="Haplotype region")
    args = parser.parse_args()

    # Load inputs
    metadata = pd.read_csv(args.metadata, sep="\t", index_col=False, dtype=str)
    hap_matrix = pd.read_csv(args.hap_matrix, sep="\t", index_col=False)

    haplotype = clean_haplotype_name(args.haplotype)
    print("Using haplotype:", haplotype)

    # Subset to haplotype row
    hap_matrix = hap_matrix[hap_matrix["GeneName"] == haplotype].copy()
    if hap_matrix.empty:
        raise SystemExit(f"ERROR: haplotype {haplotype} not found in hap-matrix GeneName column.")

    # Drop first two columns and transpose
    hap_matrix = hap_matrix.iloc[:, 2:].T

    # Turn into 2-column DF: original_id + percent_homozygosed
    hap_matrix = (
        hap_matrix
        .reset_index()
        .rename(columns={
            "index": "original_id",
            hap_matrix.columns[0]: "percent_homozygosed"
        })
    )

    hap_matrix["original_id"] = hap_matrix["original_id"].astype(str)
    metadata["original_id"] = metadata["original_id"].astype(str)

    merged = metadata.merge(hap_matrix, on="original_id", how="left")

    # cancer status
    if "original_dx" not in merged.columns:
        raise SystemExit("ERROR: metadata missing required column: original_dx")

    merged["cancer_status"] = (merged["original_dx"] != "control").astype(int)

    roc_df = merged[["original_id", "cancer_status", "percent_homozygosed", "original_dx"]].dropna().copy()
    roc_df["percent_homozygosed"] = pd.to_numeric(roc_df["percent_homozygosed"], errors="coerce")
    roc_df = roc_df.dropna(subset=["percent_homozygosed"])

    if roc_df.empty:
        raise SystemExit("ERROR: no usable rows after merge/dropna.")

    # AUROC
    auroc = roc_auc_score(roc_df["cancer_status"], roc_df["percent_homozygosed"])
    print(f"AUROC (percent_homozygosed vs cancer_status): {auroc:.4f}")

    # Outputs
    out_tsv = f"{args.haplotype}.auroc.tsv"
    out_ids = f"{args.haplotype}.percent_homozygosed_ge_auroc.ids.txt"

    roc_df.to_csv(out_tsv, sep="\t", index=False)

    # Identify cancer samples with percent_homozygosed >= AUROC
    high_hom_df = roc_df[
        (roc_df["percent_homozygosed"] >= auroc) &
        (roc_df["original_dx"] != "control")
    ]

    high_hom_df[["original_id"]].to_csv(
        out_ids,
        sep="\t",
        index=False,
        header=False
    )

    print("Wrote:", out_tsv)
    print("Wrote:", out_ids)

#python3 analysis_1d.auroc.py --metadata metadatas/kidney.metadata --hap-matrix analysis_1b_output.jan9.sw_250000.roh_100000kb.tsv.gz.tsv.gz --haplotype chr22_32200001_32450000

if __name__ == "__main__":
    main()
