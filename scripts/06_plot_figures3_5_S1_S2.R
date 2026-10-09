#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

root <- normalizePath(".")
gwas_dir <- file.path(root, "data/gwas")
out_dir <- file.path(root, "figures/generated")
geno_file <- file.path(root, "data/genotypes/YoGI_genotypes_MAF0.05_303.hmp.txt.gz")
trait_file <- file.path(root, "data/phenotypes/YoGI_traits_303.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "scripts/lib/read_compressed.R"))

trait_labels <- c(
  Awn = "Awn presence",
  Plant_height_node = "Plant height (without spike)",
  Plant_height_spike = "Plant height (with spike)",
  Spike_number = "Spike number",
  Spike_length = "Spike length",
  Spikelet_number = "Spikelet number",
  Peduncle_length = "Peduncle length",
  Leaf_length = "Leaf length"
)

result_file <- function(trait) {
  file.path(gwas_dir, sprintf("%s_GAPIT_MLM.csv.gz", trait))
}

read_gwas <- function(trait) {
  x <- read_compressed_table(result_file(trait))
  setnames(x, c("Chr", "Pos", "P.value"), c("Chromosome", "Position", "P"))
  x[, `:=`(
    Trait = trait,
    Trait_label = unname(trait_labels[trait]),
    Chromosome = as.integer(Chromosome),
    Position = as.numeric(Position),
    P = as.numeric(P),
    MAF = as.numeric(MAF),
    neglog10P = -log10(as.numeric(P))
  )]
  x
}

traits_to_read <- names(trait_labels)
gwas <- rbindlist(lapply(traits_to_read, read_gwas), use.names = TRUE, fill = TRUE)
marker_count <- uniqueN(gwas[Trait == "Awn", SNP])
stopifnot(marker_count == 32221L)
bonferroni_p <- 0.05 / marker_count
bonferroni_logp <- -log10(bonferroni_p)

chromosome_labels <- c(paste0(1:7, "A"), paste0(1:7, "B"), paste0(1:7, "D"))
chr_info <- unique(gwas[Trait == "Awn", .(Chromosome, Position)])
chr_offsets <- chr_info[, .(chr_len = max(Position, na.rm = TRUE)), by = Chromosome]
setorder(chr_offsets, Chromosome)
chr_offsets[, offset := shift(cumsum(chr_len), fill = 0)]
chr_offsets[, center := offset + chr_len / 2]
chr_offsets[, label := chromosome_labels[Chromosome]]

gwas <- merge(gwas, chr_offsets[, .(Chromosome, offset)], by = "Chromosome")
gwas[, `:=`(
  cum_pos = Position + offset,
  chr_factor = factor(Chromosome),
  significant = P < bonferroni_p
)]

chr_colors <- rep(c("#D95F5F", "#E7B353", "#74A57F", "#5E8DB8"), length.out = 21)
names(chr_colors) <- as.character(1:21)

theme_publication <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", color = "#1F1F1F"),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(color = "#303030"),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.margin = margin(7, 9, 7, 7)
    )
}

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  median(qchisq(1 - p, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
}

make_manhattan <- function(trait, label_hits = c("significant", "top", "none")) {
  label_hits <- match.arg(label_hits)
  dt <- copy(gwas[Trait == trait & is.finite(P)])
  top <- dt[which.min(P)]
  labels <- switch(
    label_hits,
    significant = dt[significant == TRUE],
    top = top,
    none = dt[0]
  )

  ymax <- max(c(dt$neglog10P, bonferroni_logp), na.rm = TRUE)
  p <- ggplot(dt, aes(cum_pos, neglog10P, color = chr_factor)) +
    geom_point(size = 0.55, alpha = 0.78, na.rm = TRUE) +
    geom_hline(
      yintercept = bonferroni_logp,
      color = "#B7242E",
      linewidth = 0.55
    ) +
    scale_color_manual(values = chr_colors, guide = "none") +
    scale_x_continuous(
      breaks = chr_offsets$center,
      labels = chr_offsets$label,
      expand = expansion(mult = c(0.004, 0.008))
    ) +
    scale_y_continuous(
      limits = c(0, ymax * 1.13),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      title = unname(trait_labels[trait]),
      x = "Chromosome",
      y = expression(-log[10](P))
    ) +
    theme_publication(11) +
    theme(
      axis.text.x = element_text(size = 8),
      axis.ticks.x = element_line(linewidth = 0.25)
    )

  if (nrow(labels) > 0) {
    labels[, label := if (label_hits == "top") paste0(SNP, "\n(top signal; not significant)") else SNP]
    p <- p +
      geom_point(
        data = labels,
        aes(cum_pos, neglog10P),
        inherit.aes = FALSE,
        shape = 21,
        size = 2.3,
        stroke = 0.7,
        color = "#79242A",
        fill = "#F5C7CA"
      ) +
      geom_text_repel(
        data = labels,
        aes(cum_pos, neglog10P, label = label),
        inherit.aes = FALSE,
        size = if (nrow(labels) > 2) 3.0 else 3.2,
        family = "Arial",
        min.segment.length = 0,
        segment.color = "#606060",
        segment.size = 0.25,
        box.padding = 0.35,
        point.padding = 0.25,
        max.overlaps = Inf,
        seed = 20260901
      )
  }
  p
}

make_qq <- function(trait, highlight = c("significant", "top")) {
  highlight <- match.arg(highlight)
  dt <- copy(gwas[Trait == trait & is.finite(P) & P > 0 & P <= 1])
  setorder(dt, P)
  n <- nrow(dt)
  dt[, rank := seq_len(.N)]
  dt[, expected_p := (rank - 0.5) / n]
  dt[, `:=`(
    expected = -log10(expected_p),
    observed = -log10(P),
    envelope_low = -log10(qbeta(0.975, rank, n - rank + 1)),
    envelope_high = -log10(qbeta(0.025, rank, n - rank + 1))
  )]
  setorder(dt, expected)
  hi <- if (highlight == "significant") dt[significant == TRUE] else dt[which.max(observed)]
  xmax <- max(dt$expected, na.rm = TRUE) * 1.04
  ymax <- max(dt$observed, na.rm = TRUE) * 1.04
  lambda <- lambda_gc(dt$P)

  ggplot(dt, aes(expected, observed)) +
    geom_ribbon(
      aes(ymin = envelope_low, ymax = envelope_high),
      fill = "#D9D9D9",
      color = NA
    ) +
    geom_abline(intercept = 0, slope = 1, color = "#C8333D", linewidth = 0.6) +
    geom_point(size = 0.65, shape = 1, stroke = 0.45, color = "#1457C5") +
    geom_point(
      data = hi,
      aes(expected, observed),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.0,
      stroke = 0.65,
      color = "#79242A",
      fill = "#F5C7CA"
    ) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
    labs(
      title = sprintf("Q-Q plot (lambda GC = %.3f)", lambda),
      x = expression(Expected~~-log[10](P)),
      y = expression(Observed~~-log[10](P))
    ) +
    theme_publication(11)
}

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 400, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".tiff")), plot, width = width, height = height, dpi = 400, compression = "lzw", bg = "white")
}

make_two_panel <- function(trait, stem, hit_mode) {
  manhattan <- make_manhattan(trait, hit_mode)
  qq <- make_qq(trait, if (hit_mode == "significant") "significant" else "top")
  combined <- manhattan + qq +
    plot_layout(widths = c(3.35, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(family = "Arial", face = "bold", size = 16))
  save_plot(combined, stem, 12, 6.25)
}

make_two_panel("Awn", "Fig3_Awn_MAF0.05_Manhattan_QQ", "significant")
make_two_panel("Plant_height_node", "Fig4_Plant_height_MAF0.05_Manhattan_QQ", "top")
make_two_panel("Spike_number", "Fig5_Spike_number_MAF0.05_Manhattan_QQ", "top")

supp_order <- c("Spike_length", "Spikelet_number", "Plant_height_spike", "Peduncle_length", "Leaf_length")
supp_labels <- c(
  Spike_length = "A. Spike length",
  Spikelet_number = "B. Spikelet number per spike",
  Plant_height_spike = "C. Plant height with spike",
  Peduncle_length = "D. Peduncle length",
  Leaf_length = "E. Leaf length"
)
supp <- copy(gwas[Trait %in% supp_order])
supp[, Panel := factor(supp_labels[Trait], levels = unname(supp_labels[supp_order]))]

p_s4 <- ggplot(supp, aes(cum_pos, neglog10P, color = chr_factor)) +
  geom_point(size = 0.36, alpha = 0.76, na.rm = TRUE) +
  geom_hline(yintercept = bonferroni_logp, color = "#B7242E", linewidth = 0.45) +
  facet_wrap(~ Panel, ncol = 1, scales = "fixed") +
  scale_color_manual(values = chr_colors, guide = "none") +
  scale_x_continuous(
    breaks = chr_offsets$center,
    labels = chr_offsets$label,
    expand = expansion(mult = c(0.004, 0.008))
  ) +
  scale_y_continuous(
    limits = c(0, max(bonferroni_logp * 1.08, max(supp$neglog10P, na.rm = TRUE) * 1.08)),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_publication(10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 10.5),
    axis.text.x = element_text(size = 7),
    panel.spacing = unit(0.8, "lines")
  )
save_plot(p_s4, "FigureS2_Other_five_traits_MAF0.05_Manhattan", 10.5, 11.5)

iupac_label <- c(
  A = "A", C = "C", G = "G", T = "T",
  R = "A/G", Y = "C/T", S = "C/G", W = "A/T", K = "G/T", M = "A/C",
  B = "C/G/T", D = "A/G/T", H = "A/C/T", V = "A/C/G", N = NA_character_
)
geno <- read_compressed_table(geno_file, check.names = FALSE)
phenotype <- fread(trait_file)
awn_hits <- gwas[Trait == "Awn" & significant == TRUE][order(P), SNP]
awn_geno <- melt(
  geno[rs %in% awn_hits],
  id.vars = names(geno)[1:11],
  variable.name = "Taxa",
  value.name = "Genotype_code"
)
awn_geno <- merge(awn_geno, phenotype[, .(Taxa, Awn)], by = "Taxa")
awn_geno[, Genotype := iupac_label[toupper(Genotype_code)]]
awn_geno <- awn_geno[!is.na(Genotype) & !is.na(Awn)]
awn_prop <- awn_geno[, .(n = .N, proportion = mean(Awn == 1)), by = .(rs, Genotype)]
awn_prop[, rs := factor(rs, levels = awn_hits)]

p_s1 <- ggplot(awn_prop, aes(Genotype, proportion)) +
  geom_col(width = 0.7, fill = "#5A6773") +
  geom_text(aes(label = paste0("n = ", n)), vjust = -0.28, size = 3.1, family = "Arial") +
  facet_wrap(~ rs, ncol = 2, scales = "free_x") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.08),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(x = "Genotype", y = "Proportion of awned accessions") +
  theme_bw(base_size = 11) +
  theme(
    text = element_text(family = "Arial", color = "#1F1F1F"),
    strip.background = element_rect(fill = "#E2E5E7", color = "#4A4A4A"),
    strip.text = element_text(size = 9.5),
    panel.grid.minor = element_blank()
  )
save_plot(p_s1, "FigS1_Awn_by_genotype_MAF0.05", 10, 6.2)

summary_table <- gwas[, {
  top <- .SD[which.min(P)]
  .(
    markers_tested = .N,
    top_SNP = top$SNP,
    top_P = top$P,
    top_neglog10P = top$neglog10P,
    top_MAF = top$MAF,
    significant_SNPs = sum(P < bonferroni_p, na.rm = TRUE),
    lambda_GC = lambda_gc(P)
  )
}, by = .(Trait, Trait_label)]
summary_table[, trait_order := match(Trait, traits_to_read)]
setorder(summary_table, trait_order)
summary_table[, trait_order := NULL]
fwrite(summary_table, file.path(out_dir, "MAF0.05_GWAS_figure_summary.csv"))
fwrite(awn_prop, file.path(out_dir, "FigureS1_source_data.csv"))

cat("Output directory:", out_dir, "\n")
cat("Markers:", marker_count, "\n")
cat("Bonferroni P:", format(bonferroni_p, scientific = TRUE, digits = 6), "\n")
cat("Bonferroni -log10(P):", sprintf("%.4f", bonferroni_logp), "\n")
print(summary_table)
