#!/usr/bin/env python3

import argparse
import pandas as pd
import re

# -----------------------------
# Argument parsing
# -----------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize Tiered variant burden in cases vs controls"
    )
    parser.add_argument(
        "--variant_matrix",
        required=True,
        help="TSV file with rows=gene_Tier_AF and columns=patient IDs"
    )
    parser.add_argument(
        "--exclude_samples",
        required=True,
        help="File with sample IDs to exclude (one per line)"
    )
    parser.add_argument(
        "--phenotypes",
        required=True,
        help="Phenotype file with Sample and original_dx columns"
    )
    parser.add_argument(
        "--out",
        default="tier_variant_case_control_summary.tsv",
        help="Output TSV"
    )
    return parser.parse_args()

# -----------------------------
# Main
# -----------------------------
def main():
    args = parse_args()

    # -----------------------------
    # Load variant matrix
    # -----------------------------
    df = pd.read_csv(args.variant_matrix, sep="\t", dtype=str)
    df = df.set_index(df.columns[0])  # gene_impact as index
    df.columns = df.columns.astype(str)

    # Convert entries to numeric
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)

    # -----------------------------
    # Load exclusions
    # -----------------------------
    exclude = (
        pd.read_csv(args.exclude_samples, header=None, dtype=str)[0]
        .str.strip()
        .astype(str)
    )

    df = df.drop(columns=[c for c in df.columns if c in set(exclude)], errors="ignore")

    # -----------------------------
    # Load phenotypes
    # -----------------------------
    pheno = pd.read_csv(args.phenotypes, sep="\t", dtype=str)
    pheno["Sample"] = pheno["Sample"].astype(str)
    pheno["is_control"] = pheno["original_dx"].str.lower().eq("control")

    pheno = pheno.set_index("Sample")

    # Keep only samples present in variant matrix
    common_samples = df.columns.intersection(pheno.index)
    df = df[common_samples]
    pheno = pheno.loc[common_samples]

    # -----------------------------
    # Parse Tier and AF from row names
    # -----------------------------
    meta = df.index.to_series().str.extract(
        r"Tier(?P<tier>[0-6])_(?P<af>001|0001)$"
    )

    df = df.loc[meta.dropna().index]
    meta = meta.dropna()

    df["Tier"] = meta["tier"].astype(int)
    df["AF"] = meta["af"]

    # -----------------------------
    # Analysis
    # -----------------------------
    results = []

    for af in ["001", "0001"]:
        for tier in range(0, 7):

            rows = df[(df["Tier"] == tier) & (df["AF"] == af)]
            if rows.empty:
                continue

            mat = rows.drop(columns=["Tier", "AF"])

            # Sample has ≥1 variant in this category
            has_variant = mat.sum(axis=0) > 0

            case_mask = ~pheno["is_control"]
            ctrl_mask = pheno["is_control"]

            results.append({
                "Tier": f"Tier{tier}",
                "AF": af,
                "Cases_with_variant": int((has_variant & case_mask).sum()),
                "Controls_with_variant": int((has_variant & ctrl_mask).sum()),
                "Total_cases": int(case_mask.sum()),
                "Total_controls": int(ctrl_mask.sum())
            })

    # -----------------------------
    # Output
    # -----------------------------
    out_df = pd.DataFrame(results)
    out_df = out_df.sort_values(["AF", "Tier"])

    out_df.to_csv(args.out, sep="\t", index=False)

    print(f"[DONE] Results written to {args.out}")

# -----------------------------
if __name__ == "__main__":
    main()

