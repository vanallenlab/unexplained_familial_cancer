#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import numpy as np
import sys, os
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
base = os.path.splitext(os.path.basename(input_file))[0]  # just filename without path or extension
output_path = "/Users/noah/Desktop/ufc_repository/results/analysis_4e_results/"

png_file = f"{output_path}pngs/{base}.results.png"
stats_file = f"{output_path}stats/{base}.results.stats"

# Extract cancer_type and PGS_ID
cancer_type = base.split('.')[0]  # "breast"
PGS_ID = base.split('.')[1]

# -------------------------------
# 2. Load and validate data
# -------------------------------
df = pd.read_csv(input_file, sep="\t")
df['group'] = df['group'].str.replace('control', 'Control')
df = df[df['group'].isin(['Control', f"Not-Inherited {cancer_type.capitalize()}", f"Familial {cancer_type.capitalize()}"])].copy()
print(f"Not-Inherited {cancer_type.capitalize()}")
# -------------------------------
# 3. Summary stats per group
# -------------------------------
summary = df.groupby('group')['PGS'].agg(['count', 'mean', 'std'])
summary_str = summary.to_string()

# -------------------------------
# 4. Compute pairwise comparisons (independent t-tests)
# -------------------------------
pairs = [
    ('Control', f"Not-Inherited {cancer_type.capitalize()}"),
    ('Control', f"Familial {cancer_type.capitalize()}"),
    (f"Not-Inherited {cancer_type.capitalize()}", f"Familial {cancer_type.capitalize()}")
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

group_order = ['Control', f"Not-Inherited {cancer_type.capitalize()}", f"Familial {cancer_type.capitalize()}"]

# --- Compute counts per group ---
counts = df['group'].value_counts().to_dict()

# --- Build custom labels with counts ---
group_labels = []
for g in group_order:
    n = counts.get(g, 0)
    label_count = f"<20" if n < 20 else str(n)
    group_labels.append(f"{g}\n(n={label_count})")

plt.figure(figsize=(4,3))

# ax = sns.boxplot(
#     data=df,
#     x='group',
#     y='PGS',
#     order=group_order,
#     palette={
#         f"Not-Inherited {cancer_type.capitalize()}": cancer_color[cancer_type],
#         f"Familial {cancer_type.capitalize()}": cancer_color[cancer_type],
#         'Control': cancer_color["control"]
#     },
#     showfliers=False
# )

# # --- Step 2: Add hatch (slashes) to the "Not-Inherited" box ---
# not_inherited_label = f"Not-Inherited {cancer_type.capitalize()}"

# for patch, label in zip(ax.artists, group_order):
#     if label == not_inherited_label:
#         patch.set_facecolor('white')          # white background
#         patch.set_edgecolor(cancer_color[cancer_type])  # colored border
#         patch.set_hatch("///")                # slashes
#         patch.set_linewidth(1.2)

# Palette must include all groups
palette_dict = {
    f"Not-Inherited {cancer_type.capitalize()}": cancer_color[cancer_type],  # initial placeholder
    f"Familial {cancer_type.capitalize()}": cancer_color[cancer_type],
    'Control': cancer_color["control"]
}

ax = sns.boxplot(
    data=df,
    x='group',
    y='PGS',
    order=group_order,
    palette=palette_dict,
    showfliers=False
)

# --- Step 2: Add hatch (slashes) to the "Not-Inherited" box ---
not_inherited_label = f"Not-Inherited {cancer_type.capitalize()}"

# Get index of the Not-Inherited box
xtick_labels = [t.get_text() for t in ax.get_xticklabels()]
x_pos = xtick_labels.index(not_inherited_label)

# Box patch
box = ax.patches[x_pos]
box.set_facecolor(cancer_color["control"])              # white background
box.set_edgecolor(cancer_color[cancer_type])  # this controls the hatch color
box.set_hatch("///")
box.set_linewidth(1.2)

# Now override the surrounding lines (the box edges) to gray
for line in ax.lines:
    # Each box has 6 lines (top, bottom, left, right, whiskers)
    # Check if the line corresponds to the box we want by x-position
    xdata = line.get_xdata()
    if np.all((xdata >= x_pos - 0.5) & (xdata <= x_pos + 0.5)):
        line.set_color("gray")          # box border color
        line.set_linewidth(1.0)

# for patch, label in zip(ax.artists, group_order):
#     if label == not_inherited_label:
#         patch.set_hatch("///")     # slashed fill
#         patch.set_edgecolor("black")
#         patch.set_linewidth(1.2)
        
#sns.boxplot(data=df, x='group', y='PGS', order=group_order, palette="Set2", showfliers=False)
#sns.stripplot(data=df, x='group', y='PGS', order=group_order, color='black', size=2, alpha=0.5)

# --- Overlay points only if group has >20 samples ---
for g in group_order:
    subset = df[df['group'] == g]
    if len(subset) > 20:
        sns.stripplot(
            data=subset, x='group', y='PGS',
            order=group_order, color='black', size=2, alpha=0.5
        )

y_max = df['PGS'].max()
y_min = df['PGS'].min()
y_range = y_max - y_min
spacing = 0.15 * y_range

#print(pairs)
for i, (g1, g2) in enumerate(pairs):
    #print(g1 + "-" + g2)
    x1, x2 = group_order.index(g1), group_order.index(g2)
    y = y_max + (i+1) * spacing
    plt.plot([x1, x1, x2, x2], [y-0.01, y, y, y-0.01], lw=1.0, c='black')

    p_val = p_values[(g1, g2)]
    p_text = "p < 1e-15" if p_val < 1e-15 else f"p = {p_val:.2e}"
    plt.text((x1+x2)/2, y + 0.01, p_text, ha='center', va='bottom', fontsize=5)

plt.title(f"PGS ({PGS_ID}) distribution by {cancer_type.capitalize()} Family History",fontsize=7)
plt.ylabel(f"PGS ({PGS_ID}) Z-Score",fontsize=5)
plt.xlabel("")

# --- Replace default x-axis labels with custom labels ---
plt.xticks(ticks=range(len(group_order)), labels=group_labels, rotation=20,fontsize=5)
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

#print(f"✅ Saved plot: {png_file}")
#print(f"✅ Saved stats: {stats_file}")
