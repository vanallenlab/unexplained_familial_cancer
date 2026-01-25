#!/usr/bin/env python3

import sys, os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import ttest_ind

# -------------------------------
# Helpers
# -------------------------------
def darken_color(color, factor=0.7):
    c = mcolors.to_rgb(color)
    return tuple(max(0, min(1, ch * factor)) for ch in c)

def parse_color_file(filepath):
    cmap = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            k, v = line.split(":", 1)
            cmap[k.strip().lower()] = v.strip().strip('"').strip("'")
    return cmap

# -------------------------------
# I/O
# -------------------------------
input_file = sys.argv[1]
base = os.path.splitext(os.path.basename(input_file))[0]

output_root = "/Users/noah/Desktop/ufc_repository/results/analysis_4e_results"
png_file = f"{output_root}/pngs/{base}.results.png"
stats_file = f"{output_root}/stats/{base}.results.stats"

cancer_type = base.split(".")[0]
PGS_ID = base.split(".")[1]

cancer_color = parse_color_file(
    "/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml"
)

# -------------------------------
# Load data
# -------------------------------
df = pd.read_csv(input_file, sep="\t")

df["group"] = df["group"].str.replace("control", "Control")

group_order = [
    "Control",
    f"Not-Inherited {cancer_type.capitalize()}",
    f"Familial {cancer_type.capitalize()}",
]

df = df[df["group"].isin(group_order)].copy()

# -------------------------------
# Summary stats
# -------------------------------
summary = df.groupby("group")["PGS"].agg(["count", "mean", "std"])

# -------------------------------
# Statistics
# -------------------------------
pairs = [
    ("Control", f"Not-Inherited {cancer_type.capitalize()}"),
    ("Control", f"Familial {cancer_type.capitalize()}"),
    (f"Not-Inherited {cancer_type.capitalize()}",
     f"Familial {cancer_type.capitalize()}"),
]

def cohen_d(x1, x2):
    n1, n2 = len(x1), len(x2)
    pooled = np.sqrt(
        ((n1 - 1) * x1.var() + (n2 - 1) * x2.var()) / (n1 + n2 - 2)
    )
    return (x1.mean() - x2.mean()) / pooled

p_values = {}
or_values = {}

for g1, g2 in pairs:
    x1 = df.loc[df.group == g1, "PGS"]
    x2 = df.loc[df.group == g2, "PGS"]

    _, p_two = ttest_ind(x1, x2, equal_var=False)
    p_one = p_two / 2 if x1.mean() < x2.mean() else 1 - p_two / 2
    p_values[(g1, g2)] = p_one

    d = cohen_d(x1, x2)
    or_values[(g1, g2)] = np.exp(d * np.pi / np.sqrt(3))

# -------------------------------
# Labels
# -------------------------------
counts = df.group.value_counts().to_dict()

fam = next(g for g in group_order if g.startswith("Familial"))
disc = next(g for g in group_order if g.startswith("Not-Inherited"))

def paired_label(n, other):
    if n <= 20:
        return "≤20"
    if other <= 20:
        return "≥20"
    return "=" + str(n)

display_labels = []
for g in group_order:
    n = counts.get(g, 0)
    if g.startswith("Not-Inherited"):
        label = f"{cancer_type.replace('_',' ').title()} Family\nDiscordant Case"
        label_n = paired_label(n, counts.get(fam, 0))
    elif g.startswith("Familial"):
        label = f"Familial {cancer_type.replace('_',' ').title()}"
        label_n = paired_label(n, counts.get(disc, 0))
    else:
        label = "Control"
        label_n = "=" + str(n)

    display_labels.append(f"{label}\n(n{label_n})")

# -----
