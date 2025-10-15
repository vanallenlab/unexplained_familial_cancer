import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# --- Hardcoded data ---
data = {
    "FDRs_with_cancer": [1, 2, 3, 4, 5, 6],
    "count": [4380, 4693, 2231, 585, 42, 2]
}
df = pd.DataFrame(data)

# --- Style setup ---
sns.set(style="whitegrid", context="talk", font_scale=1.2)

# --- Plot ---
plt.figure(figsize=(7, 5))
ax = sns.barplot(
    data=df,
    x="FDRs_with_cancer",
    y="count",
    color="steelblue",
    edgecolor="black",
)

# --- Add counts above bars ---
for i, row in df.iterrows():
    ax.text(i, row["count"] + max(df["count"]) * 0.02, f"{row['count']:,}",
            ha='center', va='bottom', fontsize=11)

# --- Labels & styling ---
ax.set_xlabel("Number of First-Degree Relatives with Cancer", fontsize=13, labelpad=10)
ax.set_ylabel("Number of Individuals", fontsize=13, labelpad=10)
ax.set_title("Distribution of Family Cancer Burden", fontsize=15, pad=15, weight='bold')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.yaxis.set_major_formatter(lambda x, _: f"{int(x):,}")

# --- Finalize ---
plt.tight_layout()
plt.savefig("family_cancer_distribution.png", dpi=300, bbox_inches="tight")
plt.show()
