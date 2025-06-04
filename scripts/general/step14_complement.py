import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
import numpy as np

# Load the data
df = pd.read_csv("af.tsv", sep='\t')

# Replace zeros with a small non-zero value to avoid log(0)
epsilon = 1e-5
df['group1_AF'] = df['group1_AF'].replace(0, epsilon)
df['group2_AF'] = df['group2_AF'].replace(0, epsilon)
df['gnomad_AF'] = df['gnomad_AF'].replace(0, epsilon)

# Define comparisons
comparisons = [
    ('group1_AF', 'group2_AF'),
    ('group1_AF', 'gnomad_AF'),
    ('group2_AF', 'gnomad_AF')
]

# Plot setup
sns.set(style="whitegrid")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for ax, (x_col, y_col) in zip(axes, comparisons):
    x = df[x_col]
    y = df[y_col]
    
    # Compute linear regression on log-log data
    slope, intercept, r_value, p_value, std_err = linregress(np.log10(x), np.log10(y))
    r2 = r_value**2

    # Plot
    ax.scatter(x, y, alpha=0.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title(f"{x_col} vs {y_col}")

    # Annotate R^2
    ax.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='left',
            fontsize=12, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="black"))

plt.tight_layout()
plt.savefig("./figures/allele_frequency_comparison.png", dpi=300)
plt.show()
