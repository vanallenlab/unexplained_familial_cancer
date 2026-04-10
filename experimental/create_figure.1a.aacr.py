import matplotlib.pyplot as plt
import numpy as np

# Data
labels = [
    "All Malignant\nNeoplastic Disease",
    "≥1 Relative",
    "≥2 Relatives",
    "≥3 Relatives",
    "No PGV"
]

counts = np.array([34932, 11933, 7553, 2860, 1496])

# Convert to percentages
percents = counts / counts[0] * 100

# Colors
colors = ["#d9d9d9", "#bdbdbd", "#969696", "#636363", "#252525"]

# Plot
fig, ax = plt.subplots(figsize=(3.5, 4.5))
bars = ax.bar(range(len(labels)), percents, color=colors)

# Labels above bars
for i, bar in enumerate(bars):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 1,
        f"{counts[i]:,}",
        ha="center",
        va="bottom",
        fontsize=7
    )

# Axis formatting
ax.set_ylabel("Percent of Cancer Patients\nin All of Us", fontsize=9)
ax.set_xticks([])
ax.set_ylim(0, 105)

# Clean % ticks (0–100 every 25)
yticks = np.arange(0, 101, 25)
ax.set_yticks(yticks)
ax.set_yticklabels([f"{y}%" for y in yticks], fontsize=9)

# Style
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Reference line + label
ax.plot([0.44, 0.7], [100, 100], color='black', linewidth=1)
ax.text(0.76, 98, "All Cancer Patients with srWGS & EHR",
        va='bottom', ha='left', fontsize=8)

# UFC cohort label
ax.text(6, 2, "UFC Cohort", ha='center', va='bottom', fontsize=9)
ax.plot([4.44, 4.7], [4.3, 4.3], color='black', linewidth=1)

# Annotations
ax.plot([0.9, 0.9], [45, 75], color='black', lw=1)
ax.plot([0.9, 3.9], [75, 75], color='black', lw=1)
ax.text(0.9, 76, "1+ Relatives with Cancer", fontsize=8)

ax.plot([2, 2], [33, 60], color='black', lw=1)
ax.plot([2, 5], [60, 60], color='black', lw=1)
ax.text(2, 61, "2+ Relatives with Cancer", fontsize=8)

ax.plot([3, 3], [18, 45], color='black', lw=1)
ax.plot([3, 6], [45, 45], color='black', lw=1)
ax.text(3, 46, "3+ Relatives with Cancer", fontsize=8)

ax.plot([4, 4], [15, 30], color='black', lw=1)
ax.plot([4, 6], [30, 30], color='black', lw=1)
ax.text(4, 31, "No PGV in a CPG", fontsize=8)

# Bottom arrows
ax.annotate('', xy=(0.05, -0.1), xytext=(0.15, -0.1),
            xycoords='axes fraction',
            arrowprops=dict(arrowstyle='-', lw=1, color='black'))

ax.annotate('', xy=(0.85, -0.1), xytext=(1.05, -0.1),
            xycoords='axes fraction',
            arrowprops=dict(arrowstyle='<-', lw=1, color='black'))

ax.text(0.5, -0.1, "Cohort Selection Process",
        transform=ax.transAxes, ha='center', va='center', fontsize=9)

plt.tight_layout()
plt.savefig("Figure1a.aacr.pdf", bbox_inches="tight", pad_inches=0)