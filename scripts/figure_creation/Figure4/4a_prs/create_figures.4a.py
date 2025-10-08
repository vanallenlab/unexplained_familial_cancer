import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Define the PGS identifiers per cancer
# PGP000186, PGP000413, PGP000542, PGP000596
cancer_configs = {
    "breast":       ["PGS000783", "PGS003380", "PGS004242", "PGS004688"],
    "prostate":     ["PGS000795", "PGS003383", "PGS004251", "PGS004694"],
    "colorectal":   ["PGS000785", "PGS003386", "PGS004243", "PGS004689"],
    "cervical": ["PGS000784","PGS003389","",""],
    "melanoma":     ["PGS000790", "PGS003382", "PGS004247",""],
    "basal_cell":     ["PGS000790", "PGS003382", "PGS004247",""],
    "squamous_cell":     ["PGS000790", "PGS003382", "PGS004247",""],
    "thyroid": ["PGS000797","","",""],
    "nervous": ["","PGS003384","",""],
    "non-hodgkins": ["PGS000791","","PGS004248",""]
}

# Set up plot style
sns.set(style="whitegrid")
plt.rcParams.update({
    'axes.labelsize': 5,
    'axes.labelweight': 'bold',
    'xtick.labelsize': 5,
    'ytick.labelsize': 5,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'axes.titlesize': 7,
    'axes.titleweight': 'bold',
})

# Directories
pgs_dir = Path("/Users/noah/Desktop/ufc_repository/results/analysis_4a_prs/")

# Set up the figure: 4 rows × 5 columns
fig, axes = plt.subplots(4, 5, figsize=(22, 16))
fig.subplots_adjust(hspace=0.6, wspace=0.4)
axes = axes.T.flatten()  # Fill left-to-right by cancer

plot_idx = 0
plot_handles_labels = None  # to store for shared legend

# Process each cancer type
for cancer, pgs_ids in cancer_configs.items():

    for pgs_id in pgs_ids:
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
        plot = sns.kdeplot(
            data=merged, x="PGS", hue="case", ax=ax,
            common_norm=False, fill=True, alpha=0.4, legend=False
        )

        # Store legend handles/labels only once
        if plot_handles_labels is None:
            handles, labels = ax.get_legend_handles_labels()
            plot_handles_labels = (handles, labels)

        ax.set_title(f"{cancer.upper()} - {pgs_id}")
        ax.set_xlabel("PGS Score")
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

#plt.tight_layout(rect=[0, 0, 0.93, 1])  # leave room for legend
plt.savefig("pgs_case_control_density_grid.png", dpi=300)
plt.show()