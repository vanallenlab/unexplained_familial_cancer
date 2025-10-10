#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import numpy as np
import sys, os, yaml
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

# -------------------------------
# 1. Input and output setup
# -------------------------------
input_file = sys.argv[1]
base = os.path.splitext(os.path.basename(input_file))[0]  # filename only
output_path = "/Users/noah/Desktop/ufc_repository/results/analysis_4e_results/"

png_file = f"{output_path}{base}.results.png"
stats_file = f"{output_path}{base}.results.stats"

# Extract cancer_type and PGS_ID
cancer_type = base.split('.')[0]  # e.g. "breast"
PGS_ID = base.split('.')[1]

# -------------------------------
# 2. Load and validate data
# -------------------------------
df = pd.read_csv(input_file, sep="\t")

# Filter to Familial and Not-Inherited
familial_label = f"Familial {cancer_type.capitalize()}"
not_inherited_label = f"Not-Inherited {cancer_type.capitalize()}"

df = df[df['group'].isin([familial_label, not_inherited_label])].copy()

# -------------------------------
# 3. Load color scheme
# -------------------------------
with open("/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml", "r") as f:
    colors = yaml.safe_load(f)

# Get familial color from cancer-specific entry
familial_color = colors.get(cancer_type.lower(), "#E69F00")  # fallback orange

# Use same color as control for Not-Inherited
control_color = colors.get("control", "#56B4E9")  # fallback blue

palette = {
    not_inherited_label: control_color,
    familial_label: familial_color
}

# -------------------------------
# 4. Summary stats per group
# -------------------------------
summary = df.groupby('group')['PGS'].agg(['count', 'mean', 'std'])
summary_str = summary.to_string()

# -------------------------------
# 5. Pairwise comparison
# -------------------------------
x1 = df.loc[df['group'] == not_inherited_label, 'PGS']
x2 = df.loc[df['group'] == familial_label, 'PGS']

# t-test (one-tailed)
stat, p_val_two_tailed = ttest_ind(x1, x2, equal_var=False)
if np.mean(x1) < np.mean(x2):
    p_val_one_tailed = p_val_two_tailed / 2
else:
    p_val_one_tailed = 1 - (p_val_two_tailed / 2)

# Cohen’s d → approximate OR proxy
def cohen_d(a, b):
    nx1, nx2 = len(a), len(b)
    pooled_std = np.sqrt(((nx1 - 1) * a.std()**2 + (nx2 - 1) * b.std()**2) / (nx1 + nx2 - 2))
    return (a.mean() - b.mean()) / pooled_std

d = cohen_d(x1, x2)
or_proxy = np.exp(d * np.pi / np.sqrt(3))

# -------------------------------
# 6. Plot
# -------------------------------
group_order = [not_inherited_label, familial_label]
counts = df['group'].value_counts().to_dict()

group_labels = []
for g in group_order:
    n = counts.get(g, 0)
    label_count = f"<20" if n < 20 else str(n)
    group_labels.append(f"{g}\n(n={label_count})")

plt.figure(figsize=(3,3))

sns.boxplot(
    data=df,
    x='group',
    y='PGS',
    order=group_order,
    palette={
        f"Not-Inherited {cancer_type.capitalize()}": cancer_color["control"],
        f"Familial {cancer_type.capitalize()}": darken_color(cancer_color[cancer_type], factor=0.6),
    },
    showfliers=False
)

#sns.boxplot(data=df, x='group', y='PGS', order=group_order, palette=palette, showfliers=False)
sns.stripplot(data=df, x='group', y='PGS', order=group_order, color='black', size=2, alpha=0.5)

y_max = df['PGS'].max()
y_min = df['PGS'].min()
y_range = y_max - y_min
spacing = 0.15 * y_range

# Add significance bar
x1_idx, x2_idx = 0, 1
y = y_max + spacing
plt.plot([x1_idx, x1_idx, x2_idx, x2_idx], [y-0.01, y, y, y-0.01], lw=1.0, c='black')
p_text = "p < 1e-15" if p_val_one_tailed < 1e-15 else f"p = {p_val_one_tailed:.2e}"
plt.text((x1_idx+x2_idx)/2, y + 0.01, p_text, ha='center', va='bottom', fontsize=6)

# Title & labels
plt.title(f"PGS ({PGS_ID}) in {cancer_type.capitalize()} by Inheritance", fontsize=7)
plt.ylabel(f"PGS ({PGS_ID}) Z-Score", fontsize=6)
plt.xlabel("")
plt.xticks(ticks=range(len(group_order)), labels=group_labels, rotation=15, fontsize=5)
plt.ylim(y_min - 0.1*y_range, y_max + 2*spacing)
plt.tight_layout()
plt.savefig(png_file, dpi=300)
plt.close()

# -------------------------------
# 7. Write results file
# -------------------------------
with open(stats_file, "w") as f:
    f.write("# Summary statistics (mean, std)\n")
    f.write(summary_str + "\n\n")
    f.write("# One-tailed t-test and Cohen’s d→OR proxy\n")
    f.write(f"{not_inherited_label} vs {familial_label}:\tp = {p_val_one_tailed:.4e},\tOR≈{or_proxy:.3f}\n")

print(f"✅ Saved plot: {png_file}")
print(f"✅ Saved stats: {stats_file}")
