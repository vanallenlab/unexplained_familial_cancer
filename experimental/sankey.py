import matplotlib.pyplot as plt
import numpy as np

# data
breast_total, breast_dx, breast_fdr = 7083, 2439, 390
kidney_total, kidney_dx, kidney_fdr = 11933, 391, 40

labels = ["Breast", "Kidney"]

totals = np.array([breast_total, kidney_total])
dx = np.array([breast_dx, kidney_dx])
fdr = np.array([breast_fdr, kidney_fdr])

y = np.arange(len(labels))

plt.rcParams["font.family"] = "Arial"

fig, ax = plt.subplots(figsize=(2,1.2))

# base bars (total)
ax.barh(y, totals, color="#d9d9d9", height=0.5)

# diagnosed overlay
ax.barh(y, dx, color="#969696", height=0.5)

# ≥3 FDR overlay
ax.barh(y, fdr, color="#252525", height=0.5)

# labels
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=6)

ax.set_xlabel("All of Us cases", fontsize=6)

ax.tick_params(axis='x', labelsize=6)

# annotate
for i in range(len(labels)):
    ax.text(dx[i]+100, i, f"{dx[i]} dx", va="center", fontsize=6)
    ax.text(fdr[i]+50, i-0.18, f"{fdr[i]} ≥3 FDR", fontsize=6)

# clean style
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(pad=0.2)
plt.savefig("cancer_summary.pdf", bbox_inches="tight",pad_inches=0, dpi=300)