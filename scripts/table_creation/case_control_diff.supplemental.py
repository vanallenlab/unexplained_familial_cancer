#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, chi2_contingency
import scipy.stats as stats
import glob
import os

# Define major ancestry groups you care about
ANCESTRY_GROUPS = [
    "European",
    "African-American",
    "East-Asian",
    "South-Asian",
    "Latin-American-1",
    "Latin-American-2"
]

# Initialize results list
results = []

for file in glob.glob("*metadata"):
    cancer = os.path.basename(file).split(".")[0].capitalize()
    df = pd.read_csv(file, sep="\t")

    # Define cases and controls
    cases = df[df["original_dx"] != "control"]
    controls = df[df["original_dx"] == "control"]

    # Helper for missing %
    def missing_pct(series):
        return 100 * series.isna().sum() / len(series)

    # ---- AGE ----
    age_missing = missing_pct(df["age"])
    age_cases = cases["age"].dropna()
    age_controls = controls["age"].dropna()
    tstat, pval = stats.ttest_ind(age_cases, age_controls, equal_var=False)
    results.append([
        cancer, "Age", "Welch_t", round(age_missing, 2), pval,
        "", "", round(age_cases.mean(), 2), round(age_controls.mean(), 2)
    ])

    if "inferred_sex" in df.columns:
        # Contingency table of cases vs controls for female
        cont = pd.crosstab(df["original_dx"] != "control", df["inferred_sex"] == "female")

        if cont.shape == (2, 2):
            chi2, pval, _, _ = chi2_contingency(cont)

            case_female_pct = 100 * cont.loc[True, True] / cont.loc[True].sum()
            control_female_pct = 100 * cont.loc[False, True] / cont.loc[False].sum()

            results.append([
                cancer,
                "Inferred_Sex_Female",
                "Chi-squared",
                0,
                pval,
                round(case_female_pct, 2),
                round(control_female_pct, 2),
                "",""
            ])

    # ---- SMOKING HISTORY ----
    if "smoking_history" in df.columns:
        smoking_df = df.dropna(subset=["smoking_history"])
        missing = missing_pct(df["smoking_history"])
        # Convert to binary smoker/non-smoker
        smoking_df["smoking_history"] = smoking_df["smoking_history"].astype(int)
        table = pd.crosstab(smoking_df["original_dx"] != "control", smoking_df["smoking_history"])
        chi2, pval, _, _ = stats.chi2_contingency(table)
        case_pct = 100 * smoking_df.loc[smoking_df["original_dx"] != "control", "smoking_history"].mean()
        control_pct = 100 * smoking_df.loc[smoking_df["original_dx"] == "control", "smoking_history"].mean()
        results.append([
            cancer, "Smoking_History", "Chi-squared", round(missing, 2), pval,
            round(case_pct, 2), round(control_pct, 2), "", ""
        ])

    # ---- GRAFPOP ANCESTRY ----
    if "grafpop_ancestry" in df.columns:
        ancestry_df = df[~df["grafpop_ancestry"].isin(["Other", "Unknown"])]
        missing = missing_pct(df["grafpop_ancestry"])

        for ancestry in ANCESTRY_GROUPS:
            ancestry_df["is_group"] = ancestry_df["grafpop_ancestry"] == ancestry
            if ancestry_df["is_group"].sum() == 0:
                continue
            table = pd.crosstab(ancestry_df["original_dx"] != "control", ancestry_df["is_group"])
            if table.shape == (2, 2):
                chi2, pval, _, _ = stats.chi2_contingency(table)
                case_pct = 100 * ancestry_df.loc[
                    ancestry_df["original_dx"] != "control", "is_group"
                ].mean()
                control_pct = 100 * ancestry_df.loc[
                    ancestry_df["original_dx"] == "control", "is_group"
                ].mean()
                results.append([
                    cancer, f"Grafpop_Ancestry_{ancestry}", "Chi-squared",
                    round(missing, 2), pval,
                    round(case_pct, 2), round(control_pct, 2), "", ""
                ])

# ---- COMBINE RESULTS ----
out_df = pd.DataFrame(results, columns=[
    "Cancer", "Variable", "Test", "Missing_%", "p-value",
    "Case_%", "Control_%", "Case_Mean", "Control_Mean"
])
out_df.to_csv("covariate_differences.tsv", sep="\t", index=False)
print("✅ Wrote covariate_differences.tsv")

