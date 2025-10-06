#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import numpy as np
import sys, os

# -------------------------------
# 1. Input and output setup
# -------------------------------
input_file = sys.argv[1]
base = os.path.splitext(input_file)[0]
output_path = "../../../../results/analysis_4c_results/"
png_file = f"{output_path + base}.results.png"
stats_file = f"{output_path + base}.results.stats"

cancer_type = f"{base.split('.')[0]}"
PGS_ID = f"{base.split('.')[1]}"

# -------------------------------
# 2. Load and validate data
# -------------------------------
df = pd.read_csv(input_file, sep="\t")
df = df[df['group'].isin(['Control', f"Sporadic_{cancer_type.capitalize()}", f"Familial_{cancer_type.capitalize()}"])].copy()

print(f"Sporadic_{cancer_type.capitalize()}")
# -------------------------------
# 3. Summary stats per group
# -------------------------------
summary = df.groupby('group')['PGS'].agg(['count', 'mean', 'std'])
summary_str = summary.to_string()

# -------------------------------
# 4. Compute pairwise comparisons (independent t-tests)
# -------------------------------
pairs = [
    ('Control', 'Sporadic_Breast'),
    ('Control', 'Familial_Breast'),
    ('Sporadic_Breast', 'Familial_Breast')
]
p_values = {}
or_values = {}

def cohen_d(x1, x2):
    """Effect size (approximate OR proxy)"""
    nx1, nx2 = len(x1), len(x2)
    pooled_std = np.sqrt(((nx1 - 1) * x1.std()**2 + (nx2 - 1) * x2.std()**2) / (nx1 + nx2 - 2))
    return (x1.mean() - x2.mean()) / pooled_std

for g1, g2 in pairs:
    x1 = df.loc[df['group'] == g1, 'PGS']
    x2 = df.loc[df['group'] == g2, 'PGS']

    # t-test
    stat, p_val_two_tailed = ttest_ind(x1, x2, equal_var=False)

    # One-tailed test for mean(x1) > mean(x2)
    if np.mean(x1) < np.mean(x2):
        p_val_one_tailed = p_val_two_tailed / 2
    else:
        p_val_one_tailed = 1 - (p_val_two_tailed / 2)
    p_values[(g1, g2)] = p_val_one_tailed

    # Cohen’s d → approximate OR
    d = cohen_d(x1, x2)
    # Convert d to OR proxy (Chinn 2000)
    or_proxy = np.exp(d * np.pi / np.sqrt(3))
    or_values[(g1, g2)] = or_proxy

# -------------------------------
# 5. Plot with significance annotations
# -------------------------------

group_order = ['Control', 'Sporadic_Breast', 'Familial_Breast']

# --- Compute counts per group ---
counts = df['group'].value_counts().to_dict()

# --- Build custom labels with counts ---
group_labels = []
for g in group_order:
    n = counts.get(g, 0)
    label_count = f"<20" if n < 20 else str(n)
    group_labels.append(f"{g}\n(n={label_count})")

plt.figure(figsize=(8,6))


sns.boxplot(data=df, x='group', y='PGS', order=group_order, palette="Set2", showfliers=False)
sns.stripplot(data=df, x='group', y='PGS', order=group_order, color='black', size=3, alpha=0.5)

y_max = df['PGS'].max()
y_min = df['PGS'].min()
y_range = y_max - y_min
spacing = 0.15 * y_range

for i, (g1, g2) in enumerate(pairs):
    x1, x2 = group_order.index(g1), group_order.index(g2)
    y = y_max + (i+1) * spacing
    plt.plot([x1, x1, x2, x2], [y-0.01, y, y, y-0.01], lw=1.5, c='black')

    p_val = p_values[(g1, g2)]
    p_text = "p < 0.001" if p_val < 0.001 else f"p = {p_val:.3f}"
    plt.text((x1+x2)/2, y + 0.01, p_text, ha='center', va='bottom', fontsize=10)

plt.title(f"PGS ({PGS_ID}) distribution by {cancer_type.capitalize()} Family History")
plt.ylabel(f"PGS ({PGS_ID})")
plt.xlabel(f"{cancer_type.capitalize()} Cohort")

# --- Replace default x-axis labels with custom labels ---
plt.xticks(ticks=range(len(group_order)), labels=group_labels, rotation=20)
plt.ylim(y_min - 0.1*y_range, y_max + len(pairs)*spacing + 0.1*y_range)
plt.tight_layout()
plt.savefig(png_file, dpi=300)
plt.close()

# -------------------------------
# 6. Write results file
# -------------------------------
with open(stats_file, "w") as f:
    f.write("# Summary statistics (mean, std)\n")
    f.write(summary_str + "\n\n")

    f.write("# Pairwise comparisons (one tailed t-test p-value and Cohen’s d→OR proxy)\n")
    for (g1, g2), p_val in p_values.items():
        f.write(f"{g1} vs {g2}:\tp = {p_val:.4e},\tOR≈{or_values[(g1,g2)]:.3f}\n")

print(f"✅ Saved plot: {png_file}")
print(f"✅ Saved stats: {stats_file}")
