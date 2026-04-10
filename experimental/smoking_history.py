import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

# ----------------------------
# Hardcode counts
# ----------------------------
our_smoker, our_non_smoker = 350, 650
aou_smoker, aou_non_smoker = 400, 600

# ----------------------------
# Percentages
# ----------------------------
our_pct = our_smoker / (our_smoker + our_non_smoker) * 100
aou_pct = aou_smoker / (aou_smoker + aou_non_smoker) * 100


def draw_cigarette(ax, x, height):

    width = 0.12

    filter_h = height * 0.18
    body_h = height * 0.75
    tip_h = height * 0.07

    # filter
    ax.add_patch(Rectangle((x - width/2, 0),
                           width, filter_h,
                           facecolor="#D8A36C",
                           edgecolor="black",
                           linewidth=0.6))

    # body
    ax.add_patch(Rectangle((x - width/2, filter_h),
                           width, body_h,
                           facecolor="white",
                           edgecolor="black",
                           linewidth=0.8))

    # burnt tip
    ax.add_patch(Rectangle((x - width/2, filter_h + body_h),
                           width, tip_h,
                           facecolor="#444444",
                           edgecolor="black",
                           linewidth=0.6))

    # smoke
    t = np.linspace(0, 3*np.pi, 200)
    smoke_x = x + 0.01 * np.sin(t)
    smoke_y = filter_h + body_h + tip_h + np.linspace(0, height * 0.25, 200)

    ax.plot(smoke_x, smoke_y, color="gray", linewidth=0.8, alpha=0.7)


# ----------------------------
# Plot
# ----------------------------

plt.rcParams["font.family"] = "Arial"

fig, ax = plt.subplots(figsize=(1, 1))

xpos = [0, 0.6]

draw_cigarette(ax, xpos[0], aou_pct)
draw_cigarette(ax, xpos[1], our_pct)

ax.set_xlim(-0.3, 0.9)
ax.set_ylim(0, 50)

ax.set_xticks(xpos)
ax.set_xticklabels(["AoU", "Our"], fontsize=6)

ax.set_ylabel("% smokers", fontsize=6)

ax.tick_params(axis='both', labelsize=6)

# clean style
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# percentage labels
ax.text(xpos[0], aou_pct + 2, f"{aou_pct:.0f}%", ha="center", fontsize=6)
ax.text(xpos[1], our_pct + 2, f"{our_pct:.0f}%", ha="center", fontsize=6)

plt.tight_layout(pad=0.2)
plt.savefig("smoking.pdf", bbox_inches="tight", pad_inches=0, dpi=300)