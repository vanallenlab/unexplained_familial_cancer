import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.patches import Patch


def parse_color_file(filepath: str) -> dict:
    """
    Parse a color mapping file and return a dictionary of name → hex color.

    Each line should look like:
      Breast: "#FFC0CB"  # pink

    The function ignores comments and blank lines.
    """
    color_map = {}

    with open(filepath, "r") as infile:
        for line in infile:
            # Strip whitespace and skip empty/comment lines
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Split key and value
            if ":" in line:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip().strip('"').strip("'")

                color_map[key.lower()] = value

    return color_map

cancer_color = parse_color_file("/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml")

or_table = pd.read_csv("/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.OR.tsv",sep='\t',index_col=False)
pval_table = pd.read_csv("/Users/noah/Desktop/ufc_repository/results/analysis_4b_all_by_all/FINAL_PRS.pvalue.tsv",sep='\t',index_col=False)


# Define the PGS identifiers per cancer
#"basal_cell":     ["PGS000790"],
cancer_configs = {
    "breast": ["PGS004688"],
    "cervix": ["PGS000784"],
    "colorectal": ["PGS000785"],
    "kidney": ["PGS000787"],
    "melanoma": ["PGS003382"],
    "non-hodgkin": ["PGS000791"],
    "prostate": ["PGS004694"],
    "thyroid": ["PGS000797"],
    "hematologic":["PGS000788"],
    "lung":["PGS000789"],
    "ovary":["PGS004249"],
    "uterus": ["PGS004244"],
    "brain": ["PGS003384"],
    "bladder":["PGS004687"]
}

# Set up plot style
#sns.set(style="whitegrid")
plt.rcParams.update({
    'axes.labelsize': 10,
    'axes.labelweight': 'bold',
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'axes.titlesize': 12,
    'axes.titleweight': 'bold',
})

# Directories
pgs_dir = Path("/Users/noah/Desktop/ufc_repository/results/analysis_4a_prs/")

# Set up the figure: 4 rows × 5 columns
fig, axes = plt.subplots(2, 7, figsize=(28, 4))
fig.subplots_adjust(hspace=0.8, wspace=0.5)
axes = axes.T.flatten()  # Fill left-to-right by cancer

plot_idx = 0
plot_handles_labels = None  # to store for shared legend

# Process each cancer type
for cancer, pgs_ids in cancer_configs.items():

    for pgs_id in pgs_ids:
        if pgs_id == "":
            plot_idx += 1
            continue

        # Store in tables using cancer and pgs_id
        or_value = float(or_table.loc[or_table['PGS_ID'] == pgs_id, cancer].values[0].split()[0])
        p_val = float(pval_table.loc[pval_table['PGS_ID'] == pgs_id, cancer])

        possible_files = [f for f in pgs_dir.glob(f"{cancer}.{pgs_id}*.anon_pgs")]
        if not possible_files:
            print(f"Missing: {cancer}.{pgs_id}*.pgs (excluding raw files)")
            continue

        pgs_path = possible_files[0]
        try:
            df = pd.read_csv(pgs_path, sep="\t")
        except Exception as e:
            print(f"Error reading {pgs_path.name}: {e}")
            continue

        if "case" not in df.columns or "PGS" not in df.columns:
            print(f"Skipping {pgs_path.name}, missing required columns.")
            continue

        ax = axes[plot_idx]

        # Define color mapping for case vs control
        palette = {
            1: cancer_color[cancer.lower()],  # case color depends on cancer type
            0: "#A9A9A9"                      # control = gray
        }

        plot = sns.kdeplot(
            data=df, x="PGS", hue="case", ax=ax,
            common_norm=False, fill=True, alpha=0.4, legend=False,palette = palette
        )
        # plot = sns.kdeplot(
        #     data=df, x="PGS",
        #     color=cancer_color[cancer.lower()],
        #     ax=ax,
        #     fill=True, alpha=0.4, common_norm=False, legend=False
        # )


        # Store legend handles/labels only once
        if plot_handles_labels is None:
            handles, labels = ax.get_legend_handles_labels()
            plot_handles_labels = (handles, labels)

        # Add text annotation inside the plot
        ax.text(
            0.95, 0.95,  # x, y in axes fraction coordinates
            f"OR = {or_value:.2f}\np = {p_val:.2e}",
            transform=ax.transAxes,  # use axes coordinates (0-1)
            ha='right', va='top',
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7)
        )
        # ax.text(
        #     0.05, 0.05,  # x=0.05, y=0.05 → bottom-left with a little padding
        #     f"OR = {or_value:.2f}\np = {p_val:.2e}",
        #     transform=ax.transAxes,  # axes coordinates
        #     ha='left', va='bottom',  # align text to the left and bottom
        #     fontsize=9,
        #     bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7)
        # )


        ax.set_title(f"{cancer.upper().replace("_"," ")} ({pgs_id})")
        ax.set_xlabel("PGS Score (Z-Score)")
        ax.set_ylabel("Density")
        ax.tick_params(width=1.5)
        plot_idx += 1

# Hide unused subplots
for i in range(plot_idx, len(axes)):
    fig.delaxes(axes[i])

# # Add one legend on the right
# if plot_handles_labels:
#     handles, labels = plot_handles_labels
#     fig.legend(
#         handles, labels,
#         loc='center right', title="Group",
#         fontsize=13, title_fontsize=14
#     )

# Define colors to match your plot
legend_elements = [
    Patch(facecolor=cancer_color["control"], linewidth=1.5, label="Control")
]

# Add one legend on the right
legend = fig.legend(
    handles=legend_elements,
    loc='center right',
    title="Group",
    fontsize=13,
    title_fontsize=12,
    frameon=True
)

# Customize legend frame and text
legend.get_frame().set_edgecolor("black")   # black border
legend.get_frame().set_linewidth(1)      # thicker frame
#legend.get_frame().set_boxstyle("round,pad=0.5")  # rounded frame

# Customize text colors and title style
for text in legend.get_texts():
    text.set_color("black")  # legend font color

legend.get_title().set_weight("bold")  # bold title
legend.get_title().set_color("black")  # title color

plt.tight_layout(rect=[0, 0, 0.95, 1])  # leave room for legend
plt.savefig("/Users/noah/Desktop/ufc_repository/results/analysis_4a_prs/Figure_4B.png", dpi=600)
#plt.show()