import gwaslab as gl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
from pathlib import Path
import matplotlib.colors as mcolors
import gzip


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
    """Open .tsv or .tsv.gz as text."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


cancer_color = parse_color_file("/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml")

input_dir = Path("/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/significant_regions/")
output_dir = Path("/Users/noah/Desktop/ufc_repository/results/analysis_1_roh/plots/")
output_dir.mkdir(parents=True, exist_ok=True)

FDR_Q = 0.1


def compute_bh_threshold(pvals, q=0.05):
    pvals = np.sort(pvals[~np.isnan(pvals)])
    pvals[pvals == 0] = 1e-300
    n = len(pvals)
    thresholds = (np.arange(1, n + 1) / n) * q
    below = pvals <= thresholds
    return pvals[below].max() if below.any() else None


# === Process each HWAS file ===
for file_path in input_dir.glob("*.tsv"):
    if not file_path.name.endswith((".tsv", ".tsv.gz")):
        continue

    base_name = file_path.name.replace(".tsv.gz", "").replace(".tsv", "")
    prefix = base_name.replace(" ", "_")

    print(f"Processing {file_path.name}")

    # read lines (handles .gz too)
    with open_textmaybe_gz(file_path) as f:
        lines = f.readlines()

    header_index = None
    for i, line in enumerate(lines):
        if line.strip().startswith("haplotype") and "beta_roh_score" in line:
            header_index = i
            break

    if header_index is None:
        print(f"  -> Could not find header line in {file_path.name}, skipping")
        continue

    df = pd.read_csv(StringIO("".join(lines[header_index:])), sep="\t")

    # Parse haplotype to chrom/pos/end
    df[["chr", "pos", "end"]] = df["haplotype"].astype(str).str.split("_", expand=True)
    df = df[df["mean_case_value"] > 0].copy()
    df["chr"] = df["chr"].str.replace("chr", "", regex=False)
    df["snp_id"] = df["chr"] + ":" + df["pos"] + "_X_X"
    df = df.drop_duplicates()

    if df.empty:
        print(f"  -> empty after filtering, skipping")
        continue

    # BH threshold (you override anyway)
    bh_cutoff = compute_bh_threshold(df["pvalue_roh_score"], q=FDR_Q)
    if bh_cutoff is None:
        bh_cutoff = 2.5e-6
    else:
        bh_cutoff = 2.5e-6  # you hardcode it

    # Create Sumstats object
    mysumstats = gl.Sumstats(
        df,
        snpid="haplotype",
        beta="beta_roh_score",
        chrom="chr",
        pos="pos",
        p="pvalue_roh_score",
        sep="\t",
    )

    cancer = file_path.stem.split(".")[0]
    color1 = cancer_color.get(cancer.lower(), "#4C72B0")
    color2 = darken_color(color1, factor=0.6)

    try:
        # IMPORTANT: close old figures so gwaslab doesn't reuse them
        plt.close("all")

        # gwaslab creates the figure internally
        out_pdf = str(output_dir / "test_mqq.pdf")
        
        # mysumstats.plot_mqq(
        #     sig_level=bh_cutoff,
        #     sig_level_lead=bh_cutoff,
        #     colors=[color1, color2],
        #     build="38",
        #     fontsize=5,          # tick labels
        #     anno_fontsize=5,     # annotation text
        #     title_fontsize=5,    # plot title
        #     font_family="Arial",
        #     marker_size=(1,1.6),
        #     save=out_pdf,
        #     save_kwargs={"dpi":400, "facecolor":"white"}  # pad_inches won't do much in v3.6.3
        # )

        out_pdf = str(output_dir / "kidney_mqq.pdf")
        out_qq_pdf = str(output_dir / "kidney_qq.pdf")
        mysumstats.plot_mqq(
            sig_level=bh_cutoff,
            mode="m",
            title = "Kidney Cancer Runs of Homozygosity GWAS",
            title_pad=1,
            sig_level_lead=bh_cutoff,
            colors=[color1, color2],
            build="38",
            fontsize=7,          # x/y tick labels
            anno_fontsize=5,     # gene/haplotype annotations
            title_fontsize=7,    # plot title
            font_family="Arial",
            marker_size=(7,7.6),
            fig_kwargs={"figsize": (8,2)},  # exact width x height in inches
            save=out_pdf,
            save_kwargs={"dpi":400, "facecolor":"white", "pad_inches":-2}  # remove extra space
        )
        mysumstats.plot_mqq(
            sig_level=bh_cutoff,
            mode="qq",
            sig_level_lead=bh_cutoff,
            colors=[color1],
            build="38",
            fontsize=7,          # x/y tick labels
            anno_fontsize=7,     # gene/haplotype annotations
            title_fontsize=7,    # plot title
            font_family="Arial",
            marker_size=(7,7.6),
            fig_kwargs={"figsize": (1.8,2)},  # exact width x height in inches
            save=out_qq_pdf,
            save_kwargs={"dpi":400, "facecolor":"white"}  # remove extra space
        )


    except Exception as e:
        print(f"Failed to plot {file_path.name}: {e}\n")
