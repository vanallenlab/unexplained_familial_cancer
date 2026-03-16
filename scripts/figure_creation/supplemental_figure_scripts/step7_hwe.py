#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import os

files = [
"/Users/noah/Desktop/ufc_repository/results/supplementary_figures/step_7_hwe/tp_prob_hwe_output_snp.tsv",
"/Users/noah/Desktop/ufc_repository/results/supplementary_figures/step_7_hwe/tp_prob_hwe_output_indel.tsv"
]

out_files = [
"/Users/noah/Desktop/ufc_repository/results/supplementary_figures/random_forest_figs/rf_hwe_snv.pdf",
"/Users/noah/Desktop/ufc_repository/results/supplementary_figures/random_forest_figs/rf_hwe_indel.pdf"
]

colors = ["#D2B48C","#8B5A2B"]

variant_types = ['SNVs','Indels']
# Colors
gray = "#4D4D4D"
euro = "#89CFF0"

# Nature-style settings
plt.rcParams.update({
"font.family": "Arial",
"font.size": 7,
"axes.linewidth": 0.6,
"xtick.major.width": 0.6,
"ytick.major.width": 0.6,
"xtick.direction": "out",
"ytick.direction": "out"
})

for i,file in enumerate(files):

    df = pd.read_csv(file, sep="\t")

    fig, ax = plt.subplots(figsize=(3.5, 2.5))

    # ----- All samples (gray) -----
    ax.plot(df["tp_prob"], df["median_p_hwe"], color=colors[i], linewidth=1.4, label="All samples")
    ax.plot(df["tp_prob"], df["median_p_hwe_1_percent"], color=colors[i], linestyle="--", linewidth=1.2, label="All samples (AF < 1%)")
    ax.plot(df["tp_prob"], df["median_p_hwe_5_percent"], color=colors[i], linestyle=":", linewidth=1.2, label="All samples (AF < 5%)")

    # ----- European subset (blue) -----
    ax.plot(df["tp_prob"], df["median_p_hwe_euro"], color=euro, linewidth=1.4, label="European")
    ax.plot(df["tp_prob"], df["median_p_hwe_1_percent_euro"], color=euro, linestyle="--", linewidth=1.2, label="European (AF < 1%)")
    ax.plot(df["tp_prob"], df["median_p_hwe_5_percent_euro"], color=euro, linestyle=":", linewidth=1.2, label="European (AF < 5%)")

    ax.set_xlabel("True positive probability threshold")
    ax.set_ylabel(f"Percentage of Biallelic {variant_types[i]}\nFailing HWE at Bonferonni Significance",fontsize=7)

    ax.legend(frameon=False)

    if variant_types[i] == "SNVs":
        ax.axvline(x=0.63, color="black", linestyle="--", linewidth=1)
        ax.text(0.63, ax.get_ylim()[1]*1, "SNV QC cutoff",
            rotation=90, va="top", ha="right", fontsize=5, color=colors[i])
    
    else:
        ax.axvline(x=0.76, color="black", linestyle="--", linewidth=1)
        ax.text(0.76, ax.get_ylim()[1]*1, "Indel QC cutoff",
            rotation=90, va="top", ha="right", fontsize=5, color=colors[i])

    # Nature-style axes
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    #ax.legend(loc="center left",frameon=False, fontsize=6)

    ax.legend(
        loc="center left",
        frameon=False,
        fontsize=5,
        ncol=2,
        handlelength=1.2,
        handletextpad=0.5,
        columnspacing=1.2
    )

    plt.tight_layout()

    out_file = file.replace(".tsv", ".pdf")
    plt.savefig(out_files[i])

    print(f"Saved: {out_files[i]}")