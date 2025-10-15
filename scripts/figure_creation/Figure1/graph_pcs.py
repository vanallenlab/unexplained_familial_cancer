import pandas as pd
import matplotlib.pyplot as plt
import yaml
import numpy as np

# --- Load input data ---
df = pd.read_csv("/Users/noah/Desktop/ufc_repository/results/step_11_pcs/merged_pcs.tsv", sep="\t")

# --- Load colors from YAML ---
with open("/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml", "r") as f:
    colors = yaml.safe_load(f)

# --- Define ancestry label mapping ---
label_map = {
    "AFR": "African",
    "AMR": "American",
    "EAS": "E. Asian",
    "EUR": "European",
    "SAS": "S. Asian"
}

# --- Map ancestry → color and pretty label ---
df["color"] = df["intake_qc_pop"].map(colors)
df["label"] = df["intake_qc_pop"].map(label_map)

# --- Compute population percentages ---
pop_counts = df["intake_qc_pop"].value_counts(normalize=True) * 100
df["pct_label"] = df["intake_qc_pop"].apply(
    lambda x: f"{label_map.get(x, x)} ({pop_counts[x]:.0f}%)" if x in pop_counts else x
)

# --- Create plot ---
fig, ax = plt.subplots(figsize=(4, 3))

# Plot each ancestry group with black-edged markers
for ancestry, subdf in df.groupby("intake_qc_pop"):
    ax.scatter(
        subdf["PC1"], subdf["PC2"],
        label=df.loc[subdf.index[0], "pct_label"],
        c=subdf["color"],
        s=10,
        alpha=0.8,
        linewidths=0.3
    )

# --- Axis labels ---
ax.set_xlabel("Principal Component 1")
ax.set_ylabel("Principal Component 2")

# --- Simplify ticks: one every 0.4 starting from 0 ---
plt.xticks(np.arange(0, plt.xlim()[1], 0.04))
plt.yticks(np.arange(0, plt.ylim()[1], 0.04))

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
plt.savefig("/Users/noah/Desktop/ufc_repository/results/step_11_pcs/pc_figure.png")
#plt.show()
