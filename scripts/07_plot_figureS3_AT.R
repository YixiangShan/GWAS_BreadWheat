library(data.table)
library(ggplot2)
library(patchwork)

root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts/lib/read_compressed.R"))
input_dir <- "data/associative_transcriptomics"
output_dir <- "figures/generated"

traits <- data.frame(
  file = c(
    "Awn_AT_results.tsv.gz",
    "Plant_height_node_AT_results.tsv.gz",
    "Spike_number_AT_results.tsv.gz",
    "Spike_length_AT_results.tsv.gz",
    "Spikelet_number_AT_results.tsv.gz",
    "Plant_height_spike_AT_results.tsv.gz",
    "Peduncle_length_AT_results.tsv.gz",
    "Leaf_length_AT_results.tsv.gz"
  ),
  label = c(
    "Awn presence",
    "Plant height without spike",
    "Spike number",
    "Spike length",
    "Spikelet number",
    "Plant height with spike",
    "Peduncle length",
    "Leaf length"
  ),
  stringsAsFactors = FALSE
)

palette <- c("1" = "#202020", "2" = "#C43C39")

make_panel <- function(filename, trait_label, panel_label) {
  dat <- read_compressed_table(file.path(input_dir, filename))
  stopifnot(nrow(dat) == 46248)

  dat[, chromosome_group := ifelse(as.integer(Graph) %% 2L == 1L, "1", "2")]
  threshold <- -log10(0.05 / nrow(dat))
  separators <- c(
    max(dat[as.integer(Graph) <= 7L, Sort]),
    max(dat[as.integer(Graph) <= 14L, Sort])
  )
  genome_centres <- c(
    mean(range(dat[as.integer(Graph) <= 7L, Sort])),
    mean(range(dat[as.integer(Graph) >= 8L & as.integer(Graph) <= 14L, Sort])),
    mean(range(dat[as.integer(Graph) >= 15L, Sort]))
  )

  ggplot(dat, aes(x = Sort, y = log10PA, colour = chromosome_group)) +
    geom_point(size = 0.32, alpha = 0.82) +
    geom_vline(xintercept = separators, linewidth = 0.25,
               linetype = "dashed", colour = "#777777") +
    geom_hline(yintercept = threshold, linewidth = 0.35,
               linetype = "dashed", colour = "#315EF5") +
    scale_colour_manual(values = palette, guide = "none") +
    scale_x_continuous(breaks = genome_centres, labels = c("A", "B", "D"),
                       expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(ylim = c(0, threshold + 1)) +
    labs(
      x = NULL,
      y = expression(-log[10](italic(P))),
      title = paste0(panel_label, ". ", trait_label)
    ) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.text = element_text(size = 6.5, colour = "black"),
      axis.title.y = element_text(size = 7),
      plot.title = element_text(size = 7.5, face = "bold", hjust = 0),
      plot.margin = margin(3, 3, 3, 3)
    )
}

plots <- lapply(seq_len(nrow(traits)), function(i) {
  make_panel(traits$file[i], traits$label[i], LETTERS[i])
})

figure <- wrap_plots(plots, ncol = 2) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
base <- file.path(output_dir, "FigureS3_AT_Manhattan_n303")

ragg::agg_png(paste0(base, ".png"), width = 183, height = 230,
              units = "mm", res = 600, background = "white")
print(figure)
dev.off()

ragg::agg_jpeg(paste0(base, ".jpeg"), width = 183, height = 230,
               units = "mm", res = 300, quality = 95,
               background = "white")
print(figure)
dev.off()

ragg::agg_tiff(paste0(base, ".tiff"), width = 183, height = 230,
               units = "mm", res = 600, compression = "lzw",
               background = "white")
print(figure)
dev.off()

grDevices::cairo_pdf(paste0(base, ".pdf"), width = 183 / 25.4,
                     height = 230 / 25.4, family = "Arial")
print(figure)
dev.off()

cat("Saved Fig. S3 to", output_dir, "\n")
