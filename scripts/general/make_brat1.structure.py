import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ==========================================================
# Load GTF
# ==========================================================

gtf = pd.read_csv(
    "/Users/noah/Desktop/BRAT1.gtf",
    sep="\t",
    header=None,
    comment="#"
)

gtf = gtf.iloc[:, :5]

gtf.columns = [
    "chrom",
    "source",
    "feature",
    "start",
    "end"
]

# ==========================================================
# First three exons in transcript order
# ==========================================================

exons = (
    gtf[gtf.feature == "exon"]
    .sort_values("start", ascending=False)
    .reset_index(drop=True)
)

exons = exons.iloc[:3].copy()

# ==========================================================
# Display coordinates
# ==========================================================

INTRON_WIDTH = 120

display_x = 0
display_exons = []

for i, exon in exons.iterrows():

    exon_len = exon.end - exon.start + 1

    display_exons.append({
        "idx": i + 1,
        "g_start": exon.start,
        "g_end": exon.end,
        "x_start": display_x,
        "x_end": display_x + exon_len
    })

    display_x += exon_len

    if i < len(exons) - 1:
        display_x += INTRON_WIDTH

plot_end = display_x

# ==========================================================
# Figure
# ==========================================================

fig, ax = plt.subplots(figsize=(1.6, 0.4))

y = 0.5
exon_height = 0.15

# ==========================================================
# Draw introns
# ==========================================================

# for i in range(len(display_exons) - 1):

#     left = display_exons[i]["x_end"]
#     right = display_exons[i + 1]["x_start"]

#     ax.plot(
#         [left, right],
#         [y, y],
#         color="black",
#         lw=1
#     )

#     ax.text(
#         (left + right) / 2,
#         y,
#         "//",
#         fontsize=4,
#         ha="center",
#         va="center",
#         backgroundcolor="white"
#     )
for i in range(len(display_exons) - 1):

    left = display_exons[i]["x_end"]
    right = display_exons[i + 1]["x_start"]

    mid = (left + right) / 2

    gap = 7  # adjust as desired

    # left half
    ax.plot(
        [left, mid - gap],
        [y, y],
        color="black",
        lw=1
    )

    # right half
    ax.plot(
        [mid + gap, right],
        [y, y],
        color="black",
        lw=1
    )

    ax.text(
        mid,
        y,
        "//",
        fontsize=4,
        ha="center",
        va="center"
    )

# ==========================================================
# Draw exons
# ==========================================================

for exon in display_exons:

    color = "black"

    # Exon 1 = orange
    if exon["idx"] == 1:
        color = "#27acfd"

    ax.add_patch(
        Rectangle(
            (exon["x_start"],
             y - exon_height/2),
            exon["x_end"] - exon["x_start"],
            exon_height,
            facecolor=color,
            edgecolor="black",
            zorder=3
        )
    )

# ==========================================================
# Exon 2 UTR (16 bp)
# ==========================================================

exon2 = display_exons[1]

bp_scale = (
    (exon2["x_end"] - exon2["x_start"])
    /
    (exon2["g_end"] - exon2["g_start"] + 1)
)

utr_start = exon2["x_start"]

utr_width = 16 * bp_scale

ax.add_patch(
    Rectangle(
        (utr_start,
         y - exon_height/2),
        utr_width,
        exon_height,
        facecolor="#27acfd",
        edgecolor="none",
        zorder=4
    )
)

# ==========================================================
# First 5 bp of intron 1
# ==========================================================

intron1_start = display_exons[0]["x_end"]

ax.add_patch(
    Rectangle(
        (intron1_start+3,
         y - 0.05),
        5,
        0.10,
        facecolor="lightgray",
        edgecolor="none",
        zorder=5
    )
)

# ==========================================================
# Appearance
# ==========================================================

ax.set_xlim(-1, plot_end + 10)
ax.set_ylim(0.4, 0.6)

ax.set_xticks([])
ax.set_yticks([])

for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()
plt.savefig("/Users/noah/Desktop/ufc_repository/results/supp_figX/BRAT1.gene_structure.pdf")
