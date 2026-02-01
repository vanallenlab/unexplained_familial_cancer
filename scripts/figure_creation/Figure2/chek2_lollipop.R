library(trackViewer)
library(GenomicRanges)
library(grid)

# ---- Input variants ----
cases <- c("chr22:28734715-G-A","chr22:28725066-G-A","chr22:28696944-T-G",
           "chr22:28689207-G-T","chr22:28695232-A-G","chr22:28725069-C-T")
controls <- c("chr22:28695190-C-A","chr22:28696912-A-G","chr22:28699843-C-G")

parse_pos <- function(x) as.integer(sub("^chr[^:]+:([0-9]+).*", "\\1", x))
all_vars <- c(cases, controls)
pos <- parse_pos(all_vars)

# ---- Variant GRanges ----
sample.gr <- GRanges(
  seqnames = "chr22",
  ranges = IRanges(start = pos, width = 1)
)
sample.gr$group <- c(rep("Case", length(cases)), rep("Control", length(controls)))

# Colors and sizes
cols <- c(Case = "#FF6F61", Control = "#B0B0B0")
sample.gr$color <- cols[sample.gr$group]

sizes <- c(Case = 0.12, Control = 0.06)
sample.gr$node.label.cex <- sizes[sample.gr$group]

# ---- Hardcoded CHEK2 CDS blocks ----
features <- GRanges(
  seqnames = "chr22",
  ranges = IRanges(
    start = c(
      28734403, 28725243, 28724977, 28719395, 28711909, 28710006, 28703505,
      28699838, 28696901, 28695710, 28695127, 28694032, 28689135, 28687900
    ),
    end = c(
      28734721, 28725367, 28725124, 28719485, 28712017, 28710059, 28703566,
      28699937, 28696987, 28695873, 28695242, 28694117, 28689215, 28687986
    )
  ),
  strand = rep("-", 14)
)

# ---- Determine x-limits for plot ----
xmin <- min(start(features), start(sample.gr))
xmax <- max(end(features), end(sample.gr))

# ---- PDF output ----
pdf("/Users/noah/Desktop/ufc_repository/results/Figure_2/CHEK2_lollipop.pdf",
    width = 4.5, height = 2)  

lolliplot(
  sample.gr, features = features,
  ylab = "",
  showAxis = TRUE,
  legend = NULL,
  xaxis.gp = gpar(fontsize = 7),
  yaxis.gp = gpar(fontsize = 7),
  xlim = c(xmin, xmax),
  margin = -8,
  base_size = 1
)

# ---- Title ----
grid.text(expression(bold("Rare (AF < 0.1%) ") * bolditalic("CHEK2") * bold(" missense variants in breast cancer")),
          x = 0.55, y = 0.73,
          gp = gpar(fontsize = 7, fontface = "bold"))


dev.off()


