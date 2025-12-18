import pandas as pd
import matplotlib.pyplot as plt
import yaml
import numpy as np

# --- Load input data ---
df = pd.read_csv("/Users/noah/Desktop/ufc_repository/results/demographics/demographics.tsv", sep="\t")

# --- Load colors from YAML ---
with open("/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml", "r") as f:
    colors = yaml.safe_load(f)

# --- Define fallback color for missing values ---
fallback_color = colors.get("Control", "#B0B0B0")

df["inferred_sex"] = df["inferred_sex"].str.capitalize()

# --- Map inferred_sex → color with fallback ---
df["color"] = df["inferred_sex"].map(colors).fillna(fallback_color)

# --- Define ancestry label mapping ---
label_map = {
    "male": "Male",
    "female": "Female",
    "other": "Other"
}

# --- Compute population percentages ---
pop_counts = df["inferred_sex"].value_counts(normalize=True) * 100
df["pct_label"] = df["inferred_sex"].apply(
    lambda x: f"{label_map.get(x, x)} ({pop_counts[x]:.0f}%)" if x in pop_counts else x
)

# --- Create plot ---
fig, ax = plt.subplots(figsize=(4, 3))

# Plot each ancestry group with black-edged markers
for ancestry, subdf in df.groupby("inferred_sex"):
    ax.scatter(
        subdf["chrX_ploidy"], subdf["chrY_ploidy"],
        label=df.loc[subdf.index[0], "pct_label"],
        c=subdf["color"],
        s=10,
        alpha=0.8,
        linewidths=0.3
    )

# --- Axis labels ---
ax.set_xlabel("chrX ploidy")
ax.set_ylabel("chrY ploidy")


# --- Keep only left and bottom spines ---
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_linewidth(1)
ax.spines["left"].set_linewidth(1)

# --- Legend: clean, without box ---
leg = ax.legend(
    loc="best",
    markerscale=2,
    frameon=False,
    fontsize=8
)

plt.tight_layout()
plt.savefig("/Users/noah/Desktop/ufc_repository/results/demographics/sex_ploidy.png",dpi=400)
#plt.show()
