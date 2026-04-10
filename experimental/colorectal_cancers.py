import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2_contingency

# -----------------------------
# Data
# -----------------------------
# AoU
aou_total = 3082
aou_fdr = 234
aou_no = aou_total - aou_fdr

# Our cohort
cohort_total = 76
cohort_fdr = 33
cohort_no = cohort_total - cohort_fdr

# -----------------------------
# Chi-square test
# -----------------------------
table = np.array([[aou_fdr, aou_no],
                  [cohort_fdr, cohort_no]])
chi2, p_value, dof, expected = chi2_contingency(table)

# -----------------------------
# Plot
# -----------------------------
labels = ["All of Us (CDRv7)\nCRC", "UFC Cohort\nCRC"]
values = [aou_fdr / aou_total, cohort_fdr / cohort_total]

plt.rcParams["font.family"] = "Arial"  # Arial substitute

fig, ax = plt.subplots(figsize=(2, 1.2))  # Width 2", height 1"

# Horizontal bar plot
bars = ax.barh(labels, values, color="#17BECF", height=0.5)

# Annotate counts inside bars
for bar, count, total in zip(bars, [aou_fdr, cohort_fdr], [aou_total, cohort_total]):
    width = bar.get_width()
    ax.text(width + 0.01, bar.get_y() + bar.get_height()/2,
            f"{count}/{total}", va='center', fontsize=5)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add gridlines for x-axis
ax.xaxis.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
ax.set_axisbelow(True)

# Limit axes
# Set x-axis limits
ax.set_xlim(0, 0.5)

# Set x-axis ticks
ax.set_xticks([0,0.25, 0.5])

# Optional: format tick labels as percentages
ax.set_xticklabels([f"{x*100:.0f}%" for x in [0,0.25, 0.5]], fontsize=5)
ax.set_xlabel("Proportion of CRC cases\nwith ≥ 1 FDR reported with CRC", fontsize=5)
ax.tick_params(axis='both', labelsize=5)

ax.text(0.75, 0.15, r"$P < 10^{-5}$", ha='center', va='bottom', transform=ax.transAxes, fontsize=5)

plt.tight_layout()
plt.savefig("/Users/noah/Desktop/ufc_repository/experimental/colorectal_epidemeology.pdf", bbox_inches='tight',pad_inches=0, dpi=300)