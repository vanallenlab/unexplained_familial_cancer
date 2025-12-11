import gwaslab as gl
import os
import pandas as pd
import gwaslab as gl
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import gzip
from pathlib import Path
import matplotlib.colors as mcolors


def darken_color(color, factor=0.7):
    """Darken a color by multiplying (0–1 scale) by `factor` < 1."""
    c = mcolors.to_rgb(color)
    return tuple(max(0, min(1, channel * factor)) for channel in c)

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

# === Configuration ===
input_dir = Path("/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/")
output_dir = Path("/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/plots/")

output_dir.mkdir(parents=True, exist_ok=True)

fdr_q = 0.1  # 5% FDR

# === Benjamini-Hochberg FDR threshold function ===
def compute_bh_threshold(pvals, q=0.05):
    pvals = np.sort(pvals[~np.isnan(pvals)])
    pvals[pvals == 0] = 1e-300
    n = len(pvals)
    thresholds = (np.arange(1, n+1) / n) * q
    below = pvals <= thresholds
    if below.any():
        return pvals[below].max()
    else:
        return None

# === Process each HWAS file ===
for file_path in input_dir.glob("*.tsv"):
    if not file_path.name.endswith((".tsv", ".tsv.gz")):
        continue


    base_name = file_path.stem.replace(".tsv", "").replace(".gz", "")
    prefix = base_name.replace(" ", "_")
    print(f"Processing {file_path.name}")

    with open(file_path, "rt") as f:  # Open as text with gzip
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.strip().startswith("haplotype") and "beta_roh_score" in line:
            header_index = i
            break

    df = pd.read_csv(StringIO("".join(lines[header_index:])), sep="\t")



    # Parse haplotype to chrom/pos/end
    df[['chr', 'pos', 'end']] = df['haplotype'].astype(str).str.split('_', expand=True)
    df = df[df['mean_case_value'] > 0]
    df['chr'] = df['chr'].str.replace("chr", "", regex=False)
    #df['pos'] = df['pos'].astype(int)
    df['snp_id'] = df['chr'] + ":" + df['pos'] + "_X_X"
    df = df.drop_duplicates()

    # Calculate BH threshold
    bh_cutoff = compute_bh_threshold(df["pvalue_roh_score"], q=fdr_q)
    bh_cutoff = 2.5e-6
    if bh_cutoff is None:
        print(f"No significant hits at {fdr_q:.2%} FDR for {file_path.name}")
        bh_cutoff = 2.5e-6  # fallback for visualization

    # Create gwaslab Sumstats object
    # mysumstats = gl.Sumstats(df,
    #     snpid="gene",
    #     beta="beta_roh_score",
    #     chrom="chr",
    #     pos="pos",
    #     p="pvalue_roh_score",
    #     sep="\t"
    # )
    mysumstats = gl.Sumstats(df,
        snpid="haplotype",
        beta="beta_roh_score",
        chrom="chr",
        pos="pos",
        p="pvalue_roh_score",
        sep="\t"
    )

    # Plot Manhattan + QQ
    cancer = file_path.stem.split('_')[1].split('.')[0].replace("non-hodgkins","Non-Hodgkin Lymphoma").replace("bladder","Bladder").replace("thyroid","Thyroid")
    title = file_path.stem.replace(".tsv", "").replace("_", " ").title()

    try:
        # Plot MQQ
        mysumstats.plot_mqq(
            title=cancer,
            sig_level=2.5e-6,
            sig_level_lead=2.5e-6,
            colors=[cancer_color[cancer.lower()],darken_color(cancer_color[cancer.lower()], factor=0.6)],
            build="38"
        )

        # Save manually
        plt.savefig(output_dir / f"{prefix}_mqq.png", dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Failed to plot {file_path.name}: {e}")

    print(f"Saved plots for {file_path.name} to {output_dir}\n")
