import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# -----------------------------
# Settings
# -----------------------------
results_dir = "/Users/noah/Desktop/ufc_repository/results/paf_results"
#cancer_order = ['Ovary','Sarcoma','Hematologic','Uterus','Colorectal','Non-Hodgkin','Thyroid','Melanoma','Lung','Cervix','Basal_Cell','Squamous_Cell','Bladder','Breast','Prostate','Neuroendocrine','Kidney']
cancer_order = ['Kidney','Neuroendocrine','Prostate','Breast','Bladder','Squamous_Cell','Basal_Cell','Cervix','Lung','Melanoma','Thyroid','Non-Hodgkin','Colorectal','Uterus','Hematologic','Sarcoma','Ovary']

def prepare_data(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t")
    df = df.sort_values("R2_full")

    baseline = df["R2_reduced"].iloc[0]

    increments = []
    previous = baseline
    for val in df["R2_full"]:
        increments.append(val - previous)
        previous = val

    df["increment"] = increments

    total_r2 = df["R2_full"].max()

    return df, baseline, total_r2


# -----------------------------
# Load all cancers first
# -----------------------------
plot_data = {}
global_max = 0

for cancer in cancer_order:
    path = os.path.join(
        results_dir,
        f"{cancer.lower()}.attributable_fraction_results.tsv"
    )

    if os.path.exists(path):
        df, baseline, total_r2 = prepare_data(path)
        plot_data[cancer] = (df, baseline, total_r2)
        global_max = max(global_max, total_r2)
    else:
        print(f"Warning: file not found for {cancer}")


# -----------------------------
# Plot combined figure
# -----------------------------
fig_height = 4 #0.5 * len(plot_data)
print(fig_height)
fig, ax = plt.subplots(figsize=(7, fig_height))

bar_height = 0.3

for i, cancer in enumerate(cancer_order):

    if cancer not in plot_data:
        continue

    df, baseline, total_r2 = plot_data[cancer]

    spacing = 0.5
    y_pos = (len(plot_data) - i - 1) * spacing  # Top-to-bottom order

    # Baseline (black)
    ax.barh(
        y_pos,
        baseline,
        height=bar_height,
        color="gray",
        edgecolor="black",
        linewidth=1.5
    )

    left = baseline
    previous_color = None
    previous_r2 = df['R2_reduced'].iloc[0]

    for _, row in df.iterrows():

        ax.barh(
            y_pos,
            row["increment"],
            left=left,
            height=bar_height,
            color=row["color"],
            edgecolor="black"
        )

        # Add label if lightgreen
        if previous_r2 is not None and row["color"] == "lightgreen" and (row['R2_full'] - previous_r2 > 0.015) or (row['added_predictor'] == "ZBP1") :
            ax.text(
                left + row["increment"] / 2,  # center of segment
                y_pos - 0.02,
                row["added_predictor"],
                ha="center",
                va="center",
                fontweight="bold",
                fontstyle="italic",
                fontsize=5,
                fontfamily="Arial"
            )
        # Add label if purple
        if row['added_predictor'] == "FH_Tier2_001":
            print(row['R2_full'] - previous_r2)
        if previous_r2 is not None and row["color"] == "plum" and (row['R2_full'] - previous_r2 > 0.012):
            ax.text(
                left + row["increment"] / 2,  # center of segment
                y_pos - 0.02,
                row["added_predictor"].split('_')[0],
                ha="center",
                va="center",
                fontweight="bold",
                fontstyle="italic",
                fontsize=5,
                fontfamily="Arial"
            )

        # Divider
        if previous_color is not None:
            if row["color"] == previous_color:
                ax.vlines(left, y_pos - bar_height/2,
                          y_pos + bar_height/2,
                          colors="black", linestyles="--", linewidth=0.01)
            else:
                ax.vlines(left, y_pos - bar_height/2,
                          y_pos + bar_height/2,
                          colors="black", linestyles="-", linewidth=1.5)

        previous_color = row["color"]
        previous_r2 = row['R2_full']
        left += row["increment"]

    # Outline full bar
    ax.add_patch(
        plt.Rectangle(
            (0, y_pos - bar_height/2),
            total_r2,
            bar_height,
            fill=False,
            edgecolor="black",
            linewidth=2
        )
    )


# -----------------------------
# Formatting
# -----------------------------
ax.set_xlim(0, 0.3)
# Full-height grid lines (0 → 0.10)
for x in np.arange(0, 0.21, 0.05):
    ax.axvline(
        x,
        color="gray",
        linestyle=":",
        linewidth=0.8,
        zorder=0,
        ymin=0,
        ymax=1
    )

# Half-height grid lines (0.15 → 0.30)
for x in np.arange(0.25, 0.301, 0.05):
    ax.axvline(
        x,
        color="gray",
        linestyle=":",
        linewidth=0.8,
        zorder=0,
        ymin=0.24,   # start halfway up
        ymax=1
    )
ax.set_yticks([i * spacing for i in range(len(plot_data))])
ax.set_yticklabels(reversed(list(plot_data.keys())),fontsize=7)
ax.set_yticklabels(
    [ {"Neuroendocrine":"NETs",
       "Squamous_Cell":"SCC",
       "Basal_Cell":"BCC",
       "Non-Hodgkin":"NHL"}.get(k, k)
      for k in reversed(list(plot_data.keys())) ],
    fontsize=7
)
ax.set_xlabel("Proportion of Variance Explained by Identified Genetic Risk Factors",fontsize=7)
#ax.set_title("Variance Explained by Genetic Predictors")

# Make a legend
from matplotlib.patches import Patch

legend_elements = [
    Patch(facecolor='gray', edgecolor='black',
          label='Sex & Genetic Ancestry',linewidth=1.5),
    Patch(facecolor='plum', edgecolor='black',
          label='Damaging Variants in CPGs',linewidth=1.5),
    Patch(facecolor='yellow', edgecolor='black',
          label='Polygenic Risk',linewidth=1.5),
    Patch(facecolor='lightgreen', edgecolor='black',
          label='Nominated CPGs',linewidth=1.5),
    Patch(facecolor='lightblue', edgecolor='black',
          label='Runs of Homozygosity',linewidth=1.5),
]

ax.legend(
    handles=legend_elements,
    loc='lower right',
    frameon=False,
    fontsize=7
)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

output_path = os.path.join(results_dir, "combined_attributable_fraction.pdf")
plt.savefig(output_path, dpi=300,bbox_inches="tight",pad_inches=0,edgecolor="none")
plt.close()

print(f"Saved combined figure: {output_path}")