import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

# ----------------------------
# Hardcode counts here
# ----------------------------
our_smoker = 350
our_non_smoker = 650

aou_smoker = 400
aou_non_smoker = 600

# ----------------------------
# Convert to percentages
# ----------------------------
our_total = our_smoker + our_non_smoker
aou_total = aou_smoker + aou_non_smoker

our_pct = our_smoker / our_total * 100
aou_pct = aou_smoker / aou_total * 100


def draw_cigarette(ax, x, height):

    width = 0.18

    filter_h = height * 0.18
    body_h = height * 0.75
    tip_h = height * 0.07

    # filter
    ax.add_patch(Rectangle((x - width/2, 0),
                           width, filter_h,
                           facecolor="#D8A36C",
                           edgecolor="black"))

    # paper
    ax.add_patch(Rectangle((x - width/2, filter_h),
                           width, body_h,
                           facecolor="white",
                           edgecolor="black",
                           linewidth=1.5))

    # burnt tip
    ax.add_patch(Rectangle((x - width/2, filter_h + body_h),
                           width, tip_h,
                           facecolor="#444444",
                           edgecolor="black"))

    # smoke
    t = np.linspace(0, 3*np.pi, 200)
    smoke_x = x + 0.02*np.sin(t)
    smoke_y = filter_h + body_h + tip_h + np.linspace(0, height*0.4, 200)

    ax.plot(smoke_x, smoke_y, color="gray", linewidth=2, alpha=0.7)


# ----------------------------
# Plot
# ----------------------------

fig, ax = plt.subplots(figsize=(5,7))

draw_cigarette(ax, 0, aou_pct)
draw_cigarette(ax, 1, our_pct)

ax.set_xlim(-0.5, 1.5)
ax.set_ylim(0, 100)

ax.set_xticks([0,1])
ax.set_xticklabels(["All of Us", "Our Cohort"], fontsize=12)

ax.set_ylabel("% Smokers")

ax.set_title("Smoking Prevalence Comparison")

# cleaner plot
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# label percentages
ax.text(0, aou_pct + 3, f"{aou_pct:.1f}%", ha='center')
ax.text(1, our_pct + 3, f"{our_pct:.1f}%", ha='center')

plt.tight_layout()
plt.show()