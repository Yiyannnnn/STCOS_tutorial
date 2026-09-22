#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: build_dlpfc_tutorial_figure.R reference.h5 coordinates.csv ",
    "simulated_counts.csv output.png"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

reference_path <- args[1]
coordinate_path <- args[2]
simulation_path <- args[3]
output_path <- args[4]
genes <- c("SLC17A7", "PCP4", "MOBP")

coordinates <- fread(coordinate_path, data.table = FALSE)
rownames(coordinates) <- coordinates$V1
coordinates$V1 <- NULL

reference <- Read10X_h5(reference_path)
simulation <- fread(
  simulation_path,
  select = c("V1", genes),
  data.table = FALSE
)
rownames(simulation) <- simulation$V1
simulation$V1 <- NULL

spots <- Reduce(
  intersect,
  list(rownames(coordinates), colnames(reference), rownames(simulation))
)
coordinates <- coordinates[spots, , drop = FALSE]
reference <- reference[genes, spots, drop = FALSE]
simulation <- simulation[spots, genes, drop = FALSE]

palette <- hcl.colors(100, "Inferno")
plot_map <- function(values, cap, title, show_y = FALSE) {
  transformed <- pmin(log1p(as.numeric(values)), cap)
  index <- cut(
    transformed,
    breaks = seq(0, cap, length.out = length(palette) + 1L),
    include.lowest = TRUE,
    labels = FALSE
  )
  plot(
    coordinates$col,
    -coordinates$row,
    pch = 16,
    cex = 0.72,
    col = palette[index],
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = title
  )
  box(col = "#CBD5E1")
  if (show_y) {
    mtext("Spatial position", side = 2, line = 0.5, cex = 0.75)
  }
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
png(output_path, width = 1800, height = 1600, res = 200, bg = "white")
par(
  mfrow = c(length(genes), 2),
  mar = c(1.1, 1.2, 2.4, 0.8),
  oma = c(1.2, 1.4, 3.8, 0.5),
  family = "sans"
)

for (gene in genes) {
  observed <- as.numeric(reference[gene, ])
  generated <- simulation[, gene]
  cap <- unname(quantile(log1p(c(observed, generated)), 0.99, na.rm = TRUE))
  if (!is.finite(cap) || cap <= 0) cap <- 1
  plot_map(observed, cap, paste(gene, "- observed"), show_y = TRUE)
  plot_map(generated, cap, paste(gene, "- ST-COS replay"))
}

mtext(
  "Reference-based output: DLPFC section 151673",
  side = 3,
  outer = TRUE,
  line = 2.0,
  font = 2,
  cex = 1.25
)
mtext(
  sprintf(
    "Common log1p count scale within each gene; %s matched spots",
    format(length(spots), big.mark = ",")
  ),
  side = 3,
  outer = TRUE,
  line = 0.6,
  cex = 0.82,
  col = "#475569"
)
dev.off()
