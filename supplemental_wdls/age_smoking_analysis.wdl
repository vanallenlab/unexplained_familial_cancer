version 1.0

workflow AgeSmokingAnalysis {
  input {
    Array[File] cancer_tsvs
  }

  scatter (tsv_file in cancer_tsvs) {
    call AnalyzeCancerFile { input: tsv = tsv_file }
  }

  call MergeReports {
    input:
      reports = AnalyzeCancerFile.report_files
  }

  output {
    Array[File] individual_reports = AnalyzeCancerFile.report_files
    File merged_summary = MergeReports.merged_summary
  }
}

# ─────────────────────────────────────────────
# Task 1: Analyze each TSV file
# ─────────────────────────────────────────────
task AnalyzeCancerFile {
  input {
    File tsv
  }

  command <<<
  set -euxo pipefail
  python3 << CODE
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, chi2_contingency
import os

df = pd.read_csv("~{tsv}", sep="\t")

# --- Split into cases vs controls ---
cases = df[df["cancer"] != "control"]
controls = df[df["cancer"] == "control"]

# --- Compute age stats ---
age_cases = cases["age"].dropna()
age_controls = controls["age"].dropna()

t_stat, p_age = ttest_ind(age_cases, age_controls, equal_var=False)

age_report = {
    "mean_age_cases": np.mean(age_cases),
    "mean_age_controls": np.mean(age_controls),
    "t_pval": p_age
}

# --- Smoking history stats ---
def summarize_smoking(sub_df):
    n_total = len(sub_df)
    n_missing = sub_df["smoking_history"].isna().sum()
    n_smokers = (sub_df["smoking_history"].str.lower() == "smoker").sum()
    n_nonsmokers = (sub_df["smoking_history"].str.lower() == "non-smoker").sum()
    return {
        "missing_pct": n_missing / n_total * 100 if n_total > 0 else np.nan,
        "smoker_pct": n_smokers / n_total * 100 if n_total > 0 else np.nan,
        "nonsmoker_pct": n_nonsmokers / n_total * 100 if n_total > 0 else np.nan,
        "n_missing": n_missing,
        "n_smokers": n_smokers,
        "n_nonsmokers": n_nonsmokers,
        "n_total": n_total
    }

smoking_cases = summarize_smoking(cases)
smoking_controls = summarize_smoking(controls)

# --- Chi-squared test for smokers vs non-smokers ---
# Ignore missing data for chi-squared
table = np.array([
    [smoking_cases["n_smokers"], smoking_cases["n_nonsmokers"]],
    [smoking_controls["n_smokers"], smoking_controls["n_nonsmokers"]]
])
chi2, p_chi, _, _ = chi2_contingency(table)

smoking_report = {
    "pval_smoking": p_chi
}

# --- Combine all stats ---
report = {
    "file": os.path.basename("~{tsv}"),
    **age_report,
    **smoking_cases,
    **{f"control_{k}": v for k, v in smoking_controls.items()},
    **smoking_report
}

out_df = pd.DataFrame([report])
out_df.to_csv("report.tsv", sep="\t", index=False)
PYCODE
  >>>

  output {
    File report_file = "report.tsv"
  }

  runtime {
    docker: "python:3.11"
    memory: "4G"
    cpu: 1
  }
}

# ─────────────────────────────────────────────
# Task 2: Merge all reports into one summary table
# ─────────────────────────────────────────────
task MergeReports {
  input {
    Array[File] reports
  }

  command <<<
    set -euxo pipefail
    python3 <<'PYCODE'
import pandas as pd
import glob

files = [f.strip() for f in open("reports.list")] if False else ~{sep=" " reports}
dfs = [pd.read_csv(f, sep="\t") for f in ~{sep=" " reports}]
merged = pd.concat(dfs, ignore_index=True)
merged.to_csv("merged_summary.tsv", sep="\t", index=False)
PYCODE
  >>>

  output {
    File merged_summary = "merged_summary.tsv"
  }

  runtime {
    docker: "python:3.11"
    memory: "2G"
  }
}
