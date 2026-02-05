import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# 1. Load results
# ------------------------------------------------------------
df = pd.read_csv("/Users/noah/Desktop/ufc_repository/results/epidemiological_results/patient_family_logistic_results.tsv", sep="\t")
df = df[df['n_intersection'] >= 5]

# ------------------------------------------------------------
# 2. Prepare volcano plot values
# ------------------------------------------------------------
# Avoid log(0) for extremely small p-values
df["p_value_clipped"] = df["p_value"].clip(lower=1e-300)

df["log2_odds_ratio"] = np.log2(df["odds_ratio"])
df["neg_log10_p"] = -np.log10(df["p_value_clipped"])

# ------------------------------------------------------------
# 3. Define coloring
#    Black if same cancer, gray otherwise
# ------------------------------------------------------------

family_to_patient_refinement = {
    "skin": {
        "melanoma",
        "basal_cell_carcinoma",
        "squamous_cell_carcinoma",
    },
    "blood_soft_tissue": {
        "hematologic",
        "non-hodgkin",
        "sarcoma",
    },
    "bone": {
        "sarcoma"
    },
    "Skin": {
        "Melanoma",
        "Basal_Cell_Carcinoma",
        "Squamous_Cell_Carcinoma",
    },
    "Blood_Soft_Tissue": {
        "Hematologic",
        "Non-Hodgkin"
    },
    "Bone": {
        "Sarcoma"
    }
}


# Default behavior: if not in mapping, keep original label
def normalize_family_cancer(row):
    patient = row["patient_cancer"]
    family = row["family_cancer"]

    # Only refine family cancer if it is coarse
    if family in family_to_patient_refinement:
        if patient in family_to_patient_refinement[family]:
            return patient  # refine family history to patient subtype

    return family


df["patient_cancer_norm"] = df["patient_cancer"]
df["family_cancer_norm"] = df.apply(normalize_family_cancer, axis=1)

print(df.head())

df["color"] = np.where(
    df["patient_cancer_norm"] == df["family_cancer_norm"],
    "black",
    "gray"
)

# df["color"] = np.where(
#     df["patient_cancer"] == df["family_cancer"],
#     "black",
#     "gray"
# )

# ------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------
plt.figure(figsize=(2.5, 2.5))

plt.scatter(
    df["log2_odds_ratio"],
    df["neg_log10_p"],
    c=df["color"],
    alpha=0.7,
    edgecolors="none",
    s = 7
)

# ------------------------------------------------------------
# 4b. Label same-cancer enrichments (black points only)
# ------------------------------------------------------------
label_df = df[df["patient_cancer_norm"] == df["family_cancer_norm"]]

for _, row in label_df.iterrows():
    label = f"{row['patient_cancer']}"
    plt.text(
        row["log2_odds_ratio"],
        row["neg_log10_p"],
        label,
        fontsize=5,
        ha="left",
        va="bottom",
        color="black",
        alpha=0.9
    )
    
# ------------------------------------------------------------
# 5. Reference lines (optional but helpful)
# ------------------------------------------------------------
# Bonferroni or nominal line (example: p = 0.05)
plt.axhline(-np.log10(0.05), linestyle="--", color="black", linewidth=1)
plt.axhline(-np.log10(0.05/175), linestyle="--", color="black", linewidth=1)
plt.axvline(0, linestyle="--", color="black", linewidth=0.8)

# ------------------------------------------------------------
# 6. Labels & aesthetics
# ------------------------------------------------------------
plt.xlabel(r"$log_{2}(OR)$", fontsize=5)
plt.ylabel(r"$-\log_{10}{P}$", fontsize=5)
plt.title("Patient Cancer ~ Family Cancer + Sex ", fontsize=7)
#plt.legend(loc="upper left")

ax = plt.gca()

##### Remove top and right spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Keep left and bottom spines
ax.spines["left"].set_visible(True)
ax.spines["bottom"].set_visible(True)

# Optional: make axes slightly thicker
ax.spines["left"].set_linewidth(1.2)
ax.spines["bottom"].set_linewidth(1.2)

# Ensure ticks are only on left/bottom
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")
#####


plt.tight_layout(pad=0.2)

# ------------------------------------------------------------
# 7. Save
# ------------------------------------------------------------
plt.savefig("/Users/noah/Desktop/ufc_repository/results/epidemiological_results/familial_cancer_volcano.png", dpi=300)
plt.savefig("/Users/noah/Desktop/ufc_repository/results/epidemiological_results/familial_cancer_volcano.pdf")
#plt.show()
plt.close()
