import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import yaml

with open("../../yamls/cancer_colors.yaml", "r") as f:
    cancer_colors = yaml.safe_load(f)


# Load data
df = pd.read_csv("dfci-ufc.sample_meta.gatkhc_posthoc_outliers.tsv", sep="\t")

# Split semicolon-separated cancer types into multiple rows
df["cancer"] = df["cancer"].astype(str)
df = df.assign(cancer=df["cancer"].str.split(";")).explode("cancer")

# Capitalize the first letter of each cancer type
df["cancer"] = df["cancer"].str.strip().str.capitalize()

# Create 'pancancer' label for non-control
df["cancer_grouped"] = df["cancer"].where(df["cancer"].str.lower() == "control", "Pancancer")

# Compute summary stats
summary = (
    df.groupby("cancer")
    .age.agg(["median", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)])
    .rename(columns={"<lambda_0>": "q1", "<lambda_1>": "q3"})
    .reset_index()
)

# Add pancancer summary
pancancer_summary = df[df["cancer"].str.lower() != "control"].age.agg(
    {"median": "median", "q1": lambda x: x.quantile(0.25), "q3": lambda x: x.quantile(0.75)}
).to_frame().T
pancancer_summary["cancer"] = "Pancancer"
summary = pd.concat([summary, pancancer_summary], ignore_index=True)

# Define order: Control, Pancancer, then alphabetical
def custom_order(cancer_list):
    cancers = sorted(set(cancer_list) - {"Control", "Pancancer", "Other","Unknown"})
    return ["Control", "Pancancer"] + cancers + ["Other", "Unknown"]

summary["cancer"] = pd.Categorical(summary["cancer"], custom_order(summary["cancer"]))
summary = summary.sort_values("cancer", ascending=False)

# Define color palette
cancer_types = summary["cancer"].tolist()
palette = {
    "Control": "gray",
    "Pancancer": "black",
    "Breast": "pink",
    "Prostate": "skyblue"
}
# Assign random distinct colors for remaining cancer types
remaining = [c for c in cancer_types if c not in palette]
extra_colors = sns.color_palette("husl", len(remaining))
palette.update(dict(zip(remaining, extra_colors)))

# Save palette to YAML
with open("cancer_colors.yaml", "w") as f:
    yaml.dump({k: str(v) for k, v in palette.items()}, f)

# Plot
plt.figure(figsize=(9, 0.5 * len(summary)))
sns.set_style("whitegrid")

for idx, row in summary.iterrows():
    cancer = row["cancer"]
    plt.plot([row["q1"], row["q3"]], [idx, idx], color=palette[cancer], lw=2)
    plt.scatter(row["median"], idx, color=palette[cancer], zorder=3)

plt.yticks(ticks=range(len(summary)), labels=summary["cancer"], fontsize = 20)
plt.xlabel("Age", fontsize = 20)
plt.ylabel("")
plt.title("Distribution of Age by Cancer Type", fontsize = 20)
plt.grid(axis='x', linestyle='--', alpha=0.3)
plt.axvline(x=summary["median"].median(), color="black", linestyle=":", lw=1)  # optional global median
sns.despine(left=True, bottom=True)
plt.tight_layout()
plt.show()
