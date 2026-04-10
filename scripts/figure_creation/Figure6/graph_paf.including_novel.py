import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# Settings
# -----------------------------
results_dir = "/Users/noah/Desktop/ufc_repository/results/paf_results"

cancer_order = [
'Kidney','Neuroendocrine','Breast_Patient_and_Family_one','Breast','Cervix','Bladder','Thyroid'
]

# -----------------------------
# Prepare data
# -----------------------------
def prepare_data(tsv_path):

    df = pd.read_csv(tsv_path, sep="\t")

    df = df[df['prevalence_model'] == "observed"]
    #df = df[~df['added_predictor'].str.contains("chr")]
    # df = df[
    #     df['added_predictor'].str.contains("SV") |
    #     df['added_predictor'].str.contains("PGS") |
    #     df['added_predictor'].str.contains("Tier")
    # ]

    df["added_predictor"] = df["added_predictor"].astype(str)

    # Color assignment
    conditions = [
        df["added_predictor"] == "PGS",
        df["added_predictor"].str.contains("SV", na=False),
        df["added_predictor"].str.contains("Tier", na=False),
        df["added_predictor"].str.contains("chr", na=False)
    ]

    choices = [
        "yellow",   # PGS
        "orange",   # SV LoF
        "plum",      # damaging variants
        "lightblue"
    ]

    df["color"] = np.select(conditions, choices, default="lightgreen")

    # Ensure rows are ordered by R2_full
    df = df.sort_values("R2_full").reset_index(drop=True)

    # -----------------------------
    # Compute incremental R2
    # -----------------------------
    increments = []
    previous_r2 = df.loc[0, "R2_reduced"]

    for r2 in df["R2_full"]:
        increments.append(r2 - previous_r2)
        previous_r2 = r2

    df["increment"] = increments

    total_r2 = df["R2_full"].max() - df["R2_reduced"].iloc[0]

    return df, total_r2


# -----------------------------
# Load all cancers
# -----------------------------
plot_data = {}
global_max = 0

for cancer in cancer_order:

    path = os.path.join(
        results_dir,
        f"{cancer.lower()}.attributable_fraction_results.adjusted.tsv"
    )

    if os.path.exists(path):

        df, total_r2 = prepare_data(path)
        plot_data[cancer] = (df, total_r2)

        global_max = max(global_max, total_r2)

    else:
        print(f"Warning: file not found for {cancer}")


# -----------------------------
# Plot
# -----------------------------
fig_height = 4
fig, ax = plt.subplots(figsize=(7, fig_height))

bar_height = 0.3
spacing = 0.5

for i, cancer in enumerate(cancer_order):

    if cancer not in plot_data:
        continue

    df, total_r2 = plot_data[cancer]

    y_pos = (len(plot_data) - i - 1) * spacing

    left = 0
    previous_color = None

    for _, row in df.iterrows():

        width = row["increment"]

        ax.barh(
            y_pos,
            width,
            left=left,
            height=bar_height,
            color=row["color"],
            edgecolor="black"
        )

        # -----------------------------
        # Labeling
        # -----------------------------
        if row["color"] == "orange" and width > 0.008:

            label = row["added_predictor"].split("_")[0]

            ax.text(
                left + width/2,
                y_pos - 0.02,
                label,
                ha="center",
                va="center",
                fontweight="bold",
                fontstyle="italic",
                fontsize=5,
                fontfamily="Arial"
            )

        if row["color"] == "plum" and width > 0.008:

            label = row["added_predictor"].split("_")[0]

            ax.text(
                left + width/2,
                y_pos - 0.02,
                label,
                ha="center",
                va="center",
                fontweight="bold",
                fontstyle="italic",
                fontsize=5,
                fontfamily="Arial"
            )

        if row["color"] == "lightblue" and width > 0.015:
            roh_region = row["added_predictor"].split('_')[0] + ":" + row["added_predictor"].split('_')[1] + "-" + row["added_predictor"].split('_')[2]
            ax.text(
                left + width/2,
                y_pos - 0.02,
                roh_region,
                ha="center",
                va="center",
                fontweight="bold",
                fontstyle="italic",
                fontsize=5,
                fontfamily="Arial"
            )

        if row["color"] == "lightgreen" and width > 0.015:

            ax.text(
                left + width/2,
                y_pos - 0.02,
                row["added_predictor"],
                ha="center",
                va="center",
                fontweight="bold",
                fontstyle="italic",
                fontsize=5,
                fontfamily="Arial"
            )

        # Divider between predictors
        if previous_color is not None:

            if row["color"] == previous_color:
                ax.vlines(
                    left,
                    y_pos - bar_height/2,
                    y_pos + bar_height/2,
                    colors="black",
                    linestyles="--",
                    linewidth=0.01
                )
            else:
                ax.vlines(
                    left,
                    y_pos - bar_height/2,
                    y_pos + bar_height/2,
                    colors="black",
                    linestyles="-",
                    linewidth=1.5
                )

        previous_color = row["color"]
        left += width

    # Outline full genetic contribution
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
ax.set_xlim(0, 0.251)
# Set x-ticks every 0.1
ax.set_xticks(np.arange(0, 0.251, 0.05))
for x in np.arange(0, 0.151, 0.05):
    ax.axvline(
        x,
        color="gray",
        linestyle=":",
        linewidth=0.8,
        zorder=0
    )
for x in np.arange(0.15, 0.201, 0.05):
    ax.axvline(
        x,
        color="gray",
        linestyle=":",
        linewidth=0.8,
        zorder=0,
        ymin=0.25,
        ymax=1
    )

for x in np.arange(0.25, 0.251, 0.05):
    ax.axvline(
        x,
        color="gray",
        linestyle=":",
        linewidth=0.8,
        zorder=0
    )
ax.set_yticks([i * spacing for i in range(len(plot_data))])

ax.set_yticklabels(
    [ {"Neuroendocrine":"NETs",
       "Breast_Patient_and_Family_one":"Familial\nBreast",
       "Basal_Cell":"BCC",
       "Non-Hodgkin":"NHL"}.get(k, k)
      for k in reversed(list(plot_data.keys())) ],
    fontsize=7
)

ax.set_xlabel(
    "Proportion of Observed Variance Explained by Genetics",
    fontsize=7
)

# -----------------------------
# Legend
# -----------------------------
from matplotlib.patches import Patch

legend_elements = [

    Patch(facecolor='orange', edgecolor='black',
          label='LoF SVs in matched CPGs', linewidth=1.5),

    Patch(facecolor='plum', edgecolor='black',
          label='Damaging Variants in matched CPGs', linewidth=1.5),

    Patch(facecolor='yellow', edgecolor='black',
          label='Polygenic Risk', linewidth=1.5),
    Patch(facecolor='lightgreen', edgecolor='black',
          label='Nominated candidate CPGs', linewidth=1.5),
    Patch(facecolor='lightblue', edgecolor='black',
          label='Nominated candidate RoH loci', linewidth=1.5)
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

output_path = os.path.join(results_dir, "combined_attributable_fraction.novel_risk.pdf")

plt.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight",
    pad_inches=0,
    edgecolor="none"
)

plt.close()

print(f"Saved combined figure: {output_path}")