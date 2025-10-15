import pandas as pd
import matplotlib.pyplot as plt
import yaml
from matplotlib.colors import to_rgb, to_hex

# --- Example data ---
data = {
    "Pancancer": 10,
    "Breast": 12.5,
    "Prostate": 8.3,
    "Lung": 15.2,
    "Colorectal": 9.7,
    "Ovary": 14.1,
    "Pancreas": 10.2,
    "Melanoma": 5.5,
    "Control": 0
}

df = pd.DataFrame(list(data.items()), columns=["Cancer", "Percent"])

# --- Load colors from YAML ---
with open("/Users/noah/Desktop/ufc_repository/yamls/color_scheme.yaml", "r") as f:
    color_dict = yaml.safe_load(f)

# --- Function to darken color slightly ---
def darken_color(color, factor=0.8):
    r, g, b = to_rgb(color)
    r, g, b = r*factor, g*factor, b*factor
    return to_hex((r, g, b))

# --- Map cancer → color and darker edge ---
df["color"] = df["Cancer"].map(lambda x: color_dict.get(x, "#B0B0B0"))
df["edge_color"] = df["color"].apply(lambda x: darken_color(x))

# --- Sort by Percent ascending for horizontal bars ---
df = df.sort_values("Percent", ascending=True)

# --- Plot horizontal bar chart ---
fig, ax = plt.subplots(figsize=(7, 5))
bars = ax.barh(
    df["Cancer"],
    df["Percent"],
    color=df["color"],
    edgecolor=df["edge_color"],
    linewidth=1.2,
    alpha=0.5
)

# --- Add percentage labels at end of each bar ---
for bar in bars:
    width = bar.get_width()
    ax.text(
        width + 0.5,  # small offset
        bar.get_y() + bar.get_height()/2,
        f"{width:.1f}%",
        va="center",
        ha="left",
        fontsize=9
    )

# --- Axis labels ---
ax.set_xlabel("Percentage of Samples with TSG Deletion", fontsize=10)
ax.set_ylabel("")

# --- Bold y-axis labels ---
ax.set_yticklabels(df["Cancer"], fontweight="bold")

# --- Remove top and right spines ---
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# --- Add subtle grid lines on x-axis ---
ax.xaxis.grid(True, linestyle="--", alpha=0.3)

# --- Ticks style ---
ax.tick_params(axis='x', direction='out', length=4)

plt.tight_layout()
plt.show()
