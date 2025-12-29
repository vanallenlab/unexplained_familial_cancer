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
for df in [df_no_prs, df_with_prs]:
    df["patient_cancer"] = df["patient_cancer"].str.lower()
    df["family_cancer"] = df["family_cancer"].str.lower()

# -------------------------------
# 3. Normalize family cancer
# -------------------------------
def normalize_family(row):
    pc = row["patient_cancer"]
    fc = row["family_cancer"]
    
    # Skin cancers
    if fc == "skin" and pc in ["melanoma", "basal_cell_carcinoma", "squamous_cell_carcinoma"]:
        return pc
    # Hematologic
    if fc == "blood_soft_tissue" and pc in ["hematologic", "non-hodgkin", "sarcoma"]:
        return pc
    return fc

for df in [df_no_prs, df_with_prs]:
    df["family_cancer"] = df.apply(normalize_family, axis=1)

# -------------------------------
# 4. Keep only patient_cancer == family_cancer
# -------------------------------
df_no_prs = df_no_prs[df_no_prs["patient_cancer"] == df_no_prs["family_cancer"]].copy()
df_with_prs = df_with_prs[df_with_prs["patient_cancer"] == df_with_prs["family_cancer"]].copy()

# -------------------------------
# 5. Merge datasets
# -------------------------------
merged = df_no_prs.merge(
    df_with_prs,
    on=["patient_cancer", "family_cancer"],
    suffixes=("_no_prs", "_with_prs")
)

# -------------------------------
# 6. Volcano plot with arrows + dots
# -------------------------------
plt.figure(figsize=(10, 8))
sns.set(style="whitegrid")

for _, row in merged.iterrows():
    x_start = np.log2(row["odds_ratio_no_prs"])
    y_start = -np.log10(row["p_value_no_prs"])
    x_end = np.log2(row["odds_ratio_with_prs"])
    y_end = -np.log10(row["p_value_with_prs"])
    
    # Arrow
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
    # Start and end dots
    plt.scatter(x_start, y_start, color="red", s=50, zorder=5)
    plt.scatter(x_end, y_end, color="green", s=50, zorder=5)
    
    # Label at end
    plt.text(x_end, y_end, row["patient_cancer"], fontsize=9, alpha=0.8)

plt.axvline(0, color="grey", linestyle="--")
plt.xlabel("log2(OR)")
plt.ylabel("-log10(p-value)")
plt.title("Family history OR before vs after including PRS")

plt.tight_layout()
plt.show()
