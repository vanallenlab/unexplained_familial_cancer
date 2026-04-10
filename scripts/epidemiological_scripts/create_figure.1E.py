import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# 1. Load results
# ------------------------------------------------------------
df = pd.read_csv("/Users/noah/Desktop/ufc_repository/results/epidemiological_results/patient_family_logistic_results.tsv", sep="\t")
df = df[df['n_intersection'] >= 4]

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

# for _, row in label_df.iterrows():
#     label = f"{row['patient_cancer']}"
#     plt.text(
#         row["log2_odds_ratio"],
#         row["neg_log10_p"],
#         label,
#         fontsize=5,
#         ha="left",
#         va="bottom",
#         color="black",
#         alpha=0.9
#     )
    
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
plt.xlabel(r"$log_{2}$(Cancer Dx given Family History)", fontsize=5)
plt.ylabel(r"$-\log_{10}({P})$", fontsize=5)
plt.title("Patient Cancer ~ FDR Cancer + Sex ", fontsize=7)
#plt.legend(loc="upper left")

ax = plt.gca()
ax.text(-1.5,-np.log10(0.05)-0.3,"Nominal Significance",fontsize=5,fontfamily="Arial")
ax.text(-1.5,-np.log10(0.05/175)-0.3,"Bonferroni Significance",fontsize=5,fontfamily="Arial")

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

ax.text(1.1,6.45,"Breast",fontsize=5,fontfamily="Arial")
ax.text(0.95,5.95,"BCC",fontsize=5,fontfamily="Arial")
ax.text(1.0,5.25,"Prostate",fontsize=5,fontfamily="Arial")
ax.text(1.5,5.35,"Lung",fontsize=5,fontfamily="Arial")
ax.text(1.3,3.3,"Colorectal",fontsize=5,fontfamily="Arial")
ax.text(0.98,2.7,"Melanoma",fontsize=5,fontfamily="Arial")
ax.text(0.40,2.25,"SCC",fontsize=5,fontfamily="Arial")
ax.text(1.65,2.3,"Kidney",fontsize=5,fontfamily="Arial")
ax.text(1.1,2.22,"NHL",fontsize=5,fontfamily="Arial")
ax.text(0.34,1.8,"Hematologic",fontsize=5,fontfamily="Arial")
ax.text(1.22,1.03,"Ovary",fontsize=5,fontfamily="Arial")
ax.text(0.82,0.4,"Cervix",fontsize=5,fontfamily="Arial")
ax.text(0.53,0.07,"Bladder",fontsize=5,fontfamily="Arial")
ax.text(1.28,1.85,"Thyroid",fontsize=5,fontfamily="Arial")

ax.text(-1.1,5,"Breast-Skin",fontsize=5,fontfamily="Arial",color="gray")
plt.tight_layout(pad=0.0)

# ------------------------------------------------------------
# 7. Save
# ------------------------------------------------------------
plt.savefig("/Users/noah/Desktop/ufc_repository/results/epidemiological_results/familial_cancer_volcano.png", dpi=300)
plt.savefig("/Users/noah/Desktop/ufc_repository/results/epidemiological_results/familial_cancer_volcano.pdf")
#plt.show()
plt.close()
