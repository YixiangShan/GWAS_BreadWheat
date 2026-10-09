suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(GGally)
  library(grid)
  library(ragg)
})

root <- normalizePath(".", mustWork = TRUE)
trait_file <- file.path(root, "data/phenotypes/YoGI_traits_303.tsv")
output_dir <- file.path(root, "figures/generated")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

phenotype <- fread(trait_file)
stopifnot(nrow(phenotype) == 303L)
stopifnot(uniqueN(phenotype$Taxa) == 303L)
stopifnot(!"YoGI_342" %in% phenotype$Taxa)

traits <- data.frame(
  `Plant height` = phenotype$Plant_height_node,
  `Peduncle length` = phenotype$Peduncle_length,
  `Leaf length` = phenotype$Leaf_length,
  `Spike number` = phenotype$Spike_number,
  `Spike length` = phenotype$Spike_length,
  `Spikelet number` = phenotype$Spikelet_number,
  `Plant height with spike` = phenotype$Plant_height_spike,
  `Awn` = phenotype$Awn_check,
  check.names = FALSE
)
stopifnot(complete.cases(traits))

pca_result <- prcomp(traits, center = TRUE, scale. = TRUE)
variance <- summary(pca_result)$importance[2, ] * 100
pc1_var <- round(variance[1], 1)
pc2_var <- round(variance[2], 1)
pc2_sign <- if (pca_result$rotation["Awn", 2] < 0) -1 else 1

pca_scores <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2] * pc2_sign,
  Taxa = phenotype$Taxa
)

loadings <- data.frame(
  PC1 = pca_result$rotation[, 1] * 5,
  PC2 = pca_result$rotation[, 2] * pc2_sign * 5,
  Trait = rownames(pca_result$rotation),
  row.names = NULL
)

p_pca <- ggplot(pca_scores, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#D7D7D7") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#D7D7D7") +
  geom_point(size = 1.15, alpha = 0.68, colour = "#6B9CC4") +
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(1.7, "mm"), type = "closed"),
    linewidth = 0.45,
    colour = "#D62728"
  ) +
  geom_text_repel(
    data = loadings,
    aes(PC1, PC2, label = Trait),
    inherit.aes = FALSE,
    colour = "#A7191F",
    size = 2.4,
    fontface = "bold",
    box.padding = 0.25,
    point.padding = 0.15,
    segment.colour = "#808080",
    segment.size = 0.25,
    seed = 20260901,
    max.overlaps = Inf
  ) +
  coord_equal() +
  labs(
    x = sprintf("PC1 (%.1f%%)", pc1_var),
    y = sprintf("PC2 (%.1f%%)", pc2_var)
  ) +
  theme_minimal(base_size = 7.5, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7, colour = "black"),
    plot.margin = margin(4, 6, 2, 7)
  )

p_corr <- ggpairs(
  traits,
  lower = list(continuous = wrap("points", alpha = 0.45, size = 0.45, colour = "#202020")),
  diag = list(continuous = wrap("densityDiag", fill = "#8A84D6", alpha = 0.55)),
  upper = list(continuous = wrap("cor", size = 2.6, colour = "#D62728")),
  columnLabels = c(
    "Plant height", "Peduncle length", "Leaf length", "Spike number",
    "Spike length", "Spikelet number", "Plant height\nwith spike", "Awn (1 or 0)"
  ),
  progress = FALSE
) +
  theme_minimal(base_size = 5.2, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.18, colour = "#E6E6E6"),
    strip.background = element_blank(),
    strip.text = element_text(size = 5.2, colour = "#202020"),
    axis.text = element_text(size = 4.5, colour = "#202020"),
    axis.title = element_text(size = 5.5),
    plot.margin = margin(1, 4, 3, 4)
  )

draw_figure <- function() {
  grid.newpage()
  layout <- grid.layout(nrow = 2, ncol = 1, heights = unit(c(0.43, 0.57), "null"))
  pushViewport(viewport(layout = layout))

  vp_a <- viewport(layout.pos.row = 1, layout.pos.col = 1, name = "panel_a")
  print(p_pca, vp = vp_a)
  pushViewport(vp_a)
  grid.text("A", x = unit(0.012, "npc"), y = unit(0.985, "npc"),
            just = c("left", "top"), gp = gpar(fontfamily = "Arial", fontsize = 9, fontface = "bold"))
  popViewport()

  vp_b <- viewport(layout.pos.row = 2, layout.pos.col = 1, name = "panel_b")
  print(p_corr, vp = vp_b)
  pushViewport(vp_b)
  grid.text("B", x = unit(0.012, "npc"), y = unit(0.985, "npc"),
            just = c("left", "top"), gp = gpar(fontfamily = "Arial", fontsize = 9, fontface = "bold"))
  popViewport(2)
}

width_in <- 183 / 25.4
height_in <- 245 / 25.4
stem <- file.path(output_dir, "Fig2_Phenotype_PCA_Correlation_n303")

cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial")
draw_figure()
dev.off()

agg_tiff(
  paste0(stem, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  compression = "lzw"
)
draw_figure()
dev.off()

agg_png(
  paste0(stem, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 300
)
draw_figure()
dev.off()

fwrite(cbind(Taxa = phenotype$Taxa, traits), file.path(output_dir, "Fig2_source_data_n303.csv"))

correlations <- as.data.table(as.table(cor(traits, use = "pairwise.complete.obs")))
setnames(correlations, c("Trait_1", "Trait_2", "Pearson_r"))
fwrite(correlations, file.path(output_dir, "Fig2_Pearson_correlations_n303.csv"))

pca_summary <- data.table(
  Metric = c("n_accessions", "excluded_accession", "PC1_percent", "PC2_percent"),
  Value = c("303", "YoGI_342", sprintf("%.1f", pc1_var), sprintf("%.1f", pc2_var))
)
fwrite(pca_summary, file.path(output_dir, "Fig2_PCA_summary_n303.csv"))

writeLines(
  c(
    "Figure conclusion: 303 accessions with matched genotype and phenotype data retain broad morphological variation.",
    "Panel A: standardized PCA biplot of eight traits.",
    "Panel B: distributions, pairwise scatterplots, and Pearson correlations for the same eight traits.",
    "Excluded accession: YoGI_342 (phenotype available, genotype unavailable).",
    "PC2 sign was oriented to keep the awn loading positive, matching the original figure; PCA axis signs are arbitrary.",
    sprintf("PC1 = %.1f%%; PC2 = %.1f%%.", pc1_var, pc2_var)
  ),
  file.path(output_dir, "Fig2_QA_notes_n303.txt")
)

cat(sprintf("n = %d\nPC1 = %.1f%%\nPC2 = %.1f%%\nOutput: %s\n", nrow(phenotype), pc1_var, pc2_var, output_dir))
