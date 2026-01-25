#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact

# -----------------------
# Load data
# -----------------------
df = pd.read_csv(
    "/Users/noah/Desktop/ufc_repository/results/analysis_5_gsea_results/tier_variant_case_control_summary.tsv",
    sep="\t"
)

df["Tier"] = df["Tier"].astype(str)
df["AF"] = df["AF"].astype(float)

TIERS = ["Tier0", "Tier1", "Tier2", "Tier3", "Tier4"]
AFS = [0.001, 0.01]

df = df[df["Tier"].isin(TIERS) & df["AF"].isin(AFS)]

# -----------------------
# Statistical test
# -----------------------
results = []

for _, row in df.iterrows():
    a = row["Cases_with_variant"]
    b = row["Total_cases"] - a
    c = row["Controls_with_variant"]
    d = row["Total_controls"] - c

    table = np.array([[a, b], [c, d]])

    try:
        chi2, p_chi, _, exp = chi2_contingency(table)

        if (exp < 5).any():
            _, p = fisher_exact(table)
            test = "Fisher"
        else:
            p = p_chi
            test = "Chi-square"

    except Exception:
        p = np.nan
        test = "NA"

    results.append({
        "Tier": row["Tier"],
        "AF": row["AF"],
        "Cases_with_variant": a,
        "Cases_without_variant": b,
        "Controls_with_variant": c,
        "Controls_without_variant": d,
        "Test": test,
        "P_value": p
    })

# -----------------------
# Write output
# -----------------------
out_df = pd.DataFrame(results)

out_path = (
    "/Users/noah/Desktop/ufc_repository/results/"
    "analysis_5_gsea_results/"
    "tier_variant_case_control_stats.tsv"
)

out_df.to_csv(out_path, sep="\t", index=False)

print(f"Saved statistical results to:\n{out_path}")
