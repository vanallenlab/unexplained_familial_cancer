import gwaslab as gl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.colors as mcolors
import gzip
from io import StringIO

# ----- Helpers -----
def darken_color(color, factor=0.7):
    c = mcolors.to_rgb(color)
    return tuple(max(0, min(1, channel * factor)) for channel in c)

def parse_color_file(filepath: str) -> dict:
    color_map = {}
    with open(filepath, "r") as infile:
        for line in infile:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" in line:
                key, value = line.split(":", 1)
                color_map[key.strip().lower()] = value.strip().strip('"').strip("'")
    return color_map

def open_textmaybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def compute_bh_threshold(pvals, q=0.05):
    pvals = np.sort(pvals[~np.isnan(pvals)])
    pvals[pvals == 0] = 1e-300
    n = len(pvals)
    thresholds = (np.arange(1, n + 1) / n) * q
    below = pvals <= thresholds
    return pvals[below].max() if below.any() else None

# ----- Settings -----
cancer_color = parse_color_file("/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml")
output_dir = Path("/Users/noah/Desktop/ufc_repository/results/supplementary_figures/roh_figures/")
output_dir.mkdir(parents=True, exist_ok=True)
FDR_Q = 0.1

# ----- Files to process -----
roh_files = [
    "/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/roh_results/bladder.sw_100000.roh_500000.tsv.gz",
    "/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/roh_results/thyroid.sw_250000.roh_500000.tsv.gz",
]

# ----- Process each file -----
for file_path in roh_files:
    file_path = Path(file_path)
    print(f"Processing {file_path.name}")

    # Read lines
    with open_textmaybe_gz(file_path) as f:
        lines = f.readlines()

    # Find header line
    header_index = None
    for i, line in enumerate(lines):
        if line.strip().startswith("haplotype") and "beta_roh_score" in line:
            header_index = i
            break
    if header_index is None:
        print(f"  -> Could not find header in {file_path.name}, skipping")
        continue

    df = pd.read_csv(StringIO("".join(lines[header_index:])), sep="\t")

    # Parse haplotype
    df[["chr", "pos", "end"]] = df["haplotype"].astype(str).str.split("_", expand=True)
    df = df[df["mean_case_value"] > 0].copy()
    df["chr"] = df["chr"].str.replace("chr", "", regex=False)
    df["snp_id"] = df["chr"] + ":" + df["pos"] + "_X_X"
    df = df.drop_duplicates()
    if df.empty:
        print(f"  -> empty after filtering, skipping")
        continue

    # BH threshold
    bh_cutoff = compute_bh_threshold(df["pvalue_roh_score"], q=FDR_Q)
    if bh_cutoff is None:
        bh_cutoff = 2.5e-6
    else:
        bh_cutoff = 2.5e-6  # you hardcode it anyway

    # Create Sumstats
    mysumstats = gl.Sumstats(
        df,
        snpid="haplotype",
        beta="beta_roh_score",
        chrom="chr",
        pos="pos",
        p="pvalue_roh_score",
        sep="\t",
    )

    # Colors
    cancer = file_path.stem.split(".")[0]
    parts = file_path.stem.replace(".tsv.gz", "").split(".")
    sliding_window = int(parts[1].replace("sw_", ""))   # 100000
    roh_length = int(parts[2].replace("roh_", ""))      # 50000
    color1 = cancer_color.get(cancer.lower(), "#4C72B0")
    color2 = darken_color(color1, factor=0.6)

    # ----- Plot -----
    plt.close("all")
    out_mqq_pdf = output_dir / f"{file_path.name.replace(".tsv.gz", "")}_m.pdf"
    try:
        mysumstats.plot_mqq(
            sig_level=bh_cutoff,
            mode="m",
            title=f"{cancer.capitalize()} Runs of Homozygosity GWAS (Sliding Window = {int(sliding_window/1000)}kb; RoH Length ≥ {int(roh_length/1000)}kb)",
            title_pad=1,
            sig_level_lead=bh_cutoff,
            colors=[color1, color2],
            build="38",
            fontsize=7,
            anno_fontsize=5,
            title_fontsize=7,
            font_family="Arial",
            marker_size=(7, 7.6),
            fig_kwargs={"figsize": (6,2)},
            save=str(out_mqq_pdf),
            save_kwargs={"dpi":400, "facecolor":"white", "pad_inches":0},
        )
        print(f"Saved M vs Q plot: {out_mqq_pdf}")
    except Exception as e:
        print(f"Failed plotting M vs Q for {file_path.name.replace(".tsv.gz", "")}: {e}")

    # ----- QQ plot -----
    out_qq_pdf = output_dir / f"{file_path.name.replace(".tsv.gz", "")}_qq.pdf"
    plt.close("all")
    try:
        mysumstats.plot_mqq(
            sig_level=bh_cutoff,
            mode="qq",
            sig_level_lead=bh_cutoff,
            colors=[color1],
            build="38",
            fontsize=5,
            anno_fontsize=5,
            title_fontsize=5,
            font_family="Arial",
            marker_size=(7, 7.6),
            fig_kwargs={"figsize": (1.6,1.6)},
            save=str(out_qq_pdf),
            save_kwargs={"dpi":400, "facecolor":"white","bbox_inches": "tight","pad_inches": 0.02},
        )
        print(f"Saved QQ plot: {out_qq_pdf}")
    except Exception as e:
        print(f"Failed plotting QQ for {file_path.name}: {e}")