import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# -------------------------------
# 1. Load files
# -------------------------------
file_no_prs = "/Users/noah/Desktop/ufc_repository/results/epidemiological_results/patient_family_logistic_results.tsv"
file_with_prs = "/Users/noah/Desktop/ufc_repository/results/epidemiological_results/patient_family_logistic_with_prs.tsv"

df_no_prs = pd.read_csv(file_no_prs, sep="\t")
df_with_prs = pd.read_csv(file_with_prs, sep="\t")

# -------------------------------
# 2. Case-insensitive matching
# -------------------------------
df_no_prs["patient_cancer"] = df_no_prs["patient_cancer"].str.lower()
df_no_prs["family_cancer"] = df_no_prs["family_cancer"].str.lower()

df_with_prs["patient_cancer"] = df_with_prs["patient_cancer"].str.lower()
df_with_prs["family_cancer"] = df_with_prs["family_cancer"].str.lower()

# -------------------------------
# 3. Keep only patient_cancer == family_cancer
# -------------------------------
df_no_prs = df_no_prs[df_no_prs["patient_cancer"] == df_no_prs["family_cancer"]].copy()
df_with_prs = df_with_prs[df_with_prs["patient_cancer"] == df_with_prs["family_cancer"]].copy()

# -------------------------------
# 4. Merge datasets
# -------------------------------
merged = df_no_prs.merge(
    df_with_prs,
    on=["patient_cancer", "family_cancer"],
    suffixes=("_no_prs", "_with_prs")
)

# -------------------------------
# 5. Volcano plot: log2(OR) vs -log10(p)
# -------------------------------
plt.figure(figsize=(10, 8))
sns.set(style="whitegrid")

for _, row in merged.iterrows():
    # Arrow from OR/p before PRS to after PRS
    x_start = np.log2(row["odds_ratio_no_prs"])
    y_start = -np.log10(row["p_value_no_prs"])
    x_end = np.log2(row["odds_ratio_with_prs"])
    y_end = -np.log10(row["p_value_with_prs"])
    
    plt.arrow(
        x_start, y_start,
        x_end - x_start,
        y_end - y_start,
        length_includes_head=True,
        head_width=0.05,
        head_length=0.05,
        fc="blue",
        ec="blue",
        alpha=0.7
    )
    # Optional: label each arrow with cancer
    plt.text(x_end, y_end, row["patient_cancer"], fontsize=9, alpha=0.8)

plt.axvline(0, color="grey", linestyle="--")
plt.xlabel("log2(OR)")
plt.ylabel("-log10(p-value)")
plt.title("Family history OR before vs after including PRS")

plt.tight_layout()
plt.show()
