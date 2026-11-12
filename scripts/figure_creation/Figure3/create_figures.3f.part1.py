#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.colors as mcolors

# -----------------------
# Adjustable color
# -----------------------
CASE_COLOR = "#FF6F61"   # change this to whatever you want
CONTROL_COLOR = "#B0B0B0"

def darken_color(color, factor=0.7):
    rgb = mcolors.to_rgb(color)
    return tuple(channel * factor for channel in rgb)

# -----------------------
# Input file
# -----------------------
pgs_file = "/Users/noah/Desktop/ufc_repository/results/analysis_4c_results/4c_prs_results/breast.PGS004688.analysis_4c.tsv"
df = pd.read_csv(pgs_file, sep="\t")

def convert_OR_to_SD(odds_ratio):
    numerator = np.log10(odds_ratio)
    denomanator = np.log10(1.782)
    return numerator/denomanator
# -----------------------
# Gene markers (EDIT HERE)
# -----------------------
gene_lines = [
    ("ATM", convert_OR_to_SD(1.82)),
    ("BARD1",convert_OR_to_SD(1.37)),
    ("BRCA1",convert_OR_to_SD(7.62)),
    ("BRCA2",convert_OR_to_SD(5.23)),
    ("CDH1",convert_OR_to_SD(2.50)),
    ("CHEK2",convert_OR_to_SD(2.47)),
    ("NF1",convert_OR_to_SD(1.93)),
    ("PALB2",convert_OR_to_SD(3.83)),
    ("RAD51C",convert_OR_to_SD(1.2)),
    ("RAD51D",convert_OR_to_SD(1.72))
]

results = []
for gene,SD in gene_lines:
    percent_above_gene_case = len(df[(df['group'] != "Control") & (df['PGS'] >= SD)])/len(df[(df['group'] != "Control")])
    percent_above_gene_isolated = len(df[(df['group'] == "Isolated Breast") & (df['PGS'] >= SD)])/len(df[(df['group'] == "Isolated Breast")])
    percent_above_gene_familial = len(df[(df['group'] == "Familial Breast") & (df['PGS'] >= SD)])/len(df[(df['group'] == "Familial Breast")])
    percent_above_gene_control = len(df[(df['group'] == "Control") & (df['PGS'] >= SD)])/len(df[(df['group'] == "Control")])
    results.append([gene,SD,percent_above_gene_control,percent_above_gene_case,percent_above_gene_isolated,percent_above_gene_familial])

results_df = pd.DataFrame(results,columns=[
        "Gene",
        "SD_for_gene",
        "control_pct_above",
        "case_pct_above",
        "isolated_pct_above",
        "familial_pct_above"
    ])

results_df.to_csv("/Users/noah/Desktop/ufc_repository/results/misc_results/Fig12.tsv",sep='\t',index=False)