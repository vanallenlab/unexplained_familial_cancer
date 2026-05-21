#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.colors as mcolors
# https://www.sciencedirect.com/science/article/pii/S1542356520316645
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
pgs_file = "/Users/noah/Desktop/ufc_repository/results/analysis_4c_results/4c_prs_results/colorectal.PGS000785.analysis_4c.tsv"
df = pd.read_csv(pgs_file, sep="\t")

def convert_OR_to_SD(odds_ratio):
    numerator = np.log10(odds_ratio)
    denomanator = np.log10(1.493)
    return numerator/denomanator
# -----------------------
# Gene markers (EDIT HERE)
# -----------------------
gene_lines = [
    ("MSH2", convert_OR_to_SD(18.1)),
    ("APC",convert_OR_to_SD(49.4)),
    ("MLH1",convert_OR_to_SD(8.6)),
    ("MSH6",convert_OR_to_SD(4.9)),
    ("EPCAM",convert_OR_to_SD(1.9)),
    ("MUTYH",convert_OR_to_SD(1.3)),
    ("PMS2",convert_OR_to_SD(1.7)),
]

results = []
for gene,SD in gene_lines:
    percent_above_gene_case = len(df[(df['group'] != "Control") & (df['PGS'] >= SD)])/len(df[(df['group'] != "Control")])
    percent_above_gene_isolated = len(df[(df['group'] == "Isolated Colorectal") & (df['PGS'] >= SD)])/len(df[(df['group'] == "Isolated Colorectal")])
    percent_above_gene_familial = len(df[(df['group'] == "Familial Colorectal") & (df['PGS'] >= SD)])/len(df[(df['group'] == "Familial Colorectal")])
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

results_df.to_csv("/Users/noah/Desktop/ufc_repository/results/analysis_4f_prs_utility/4f.colorectal.tsv",sep='\t',index=False)