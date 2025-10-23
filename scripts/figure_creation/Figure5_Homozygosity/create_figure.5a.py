import gwaslab as gl
import os
import pandas as pd
import gwaslab as gl
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import gzip
from pathlib import Path

# === Configuration ===
input_dir = Path("/Users/noah/Downloads/1C_RESULTS_10kb/")
output_dir = Path("/Users/noah/Downloads/1C_RESULTS_10kb/plots/")

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
for file_path in input_dir.glob("*.tsv.gz"):
    if not file_path.name.endswith((".tsv", ".tsv.gz")):
        continue


    base_name = file_path.stem.replace(".tsv", "").replace(".gz", "")
    prefix = base_name.replace(" ", "_")
    print(f"Processing {file_path.name}")

    with gzip.open(file_path, "rt") as f:  # Open as text with gzip
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.strip().startswith("haplotype") and "beta_roh_score" in line:
            header_index = i
            break

    df = pd.read_csv(StringIO("".join(lines[header_index:])), sep="\t")



    # Parse haplotype to chrom/pos/end
    df[['chr', 'pos', 'end','gene']] = df['haplotype'].astype(str).str.split('_', expand=True)
    df['snp_id'] = df['chr'] + ":" + df['pos'] + "_X_X"
    df = df.drop_duplicates()

    # Calculate BH threshold
    bh_cutoff = compute_bh_threshold(df["pvalue_roh_score"], q=fdr_q)
    if bh_cutoff is None:
        print(f"No significant hits at {fdr_q:.2%} FDR for {file_path.name}")
        bh_cutoff = 0.05/13000  # fallback for visualization

    # Create gwaslab Sumstats object
    mysumstats = gl.Sumstats(df,
        snpid="gene",
        beta="beta_roh_score",
        chrom="chr",
        pos="pos",
        p="pvalue_roh_score",
        sep="\t"
    )

    # Plot Manhattan + QQ
    title = file_path.stem.replace(".tsv", "").replace("_", " ").title()

    try:
        # Plot MQQ
        mysumstats.plot_mqq(
            title=title,
            sig_level=bh_cutoff,
            sig_level_lead=bh_cutoff,
            anno="GENENAME",
            build="38"
        )

        # Save manually
        plt.savefig(output_dir / f"{prefix}_mqq.png", dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Failed to plot {file_path.name}: {e}")

    print(f"Saved plots for {file_path.name} to {output_dir}\n")
