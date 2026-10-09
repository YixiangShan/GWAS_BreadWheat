#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(scales)
  library(tidyr)
})

root <- normalizePath(".", mustWork = TRUE)
input_dir <- file.path(root, "data/population_structure")
output_dir <- file.path(root, "figures/generated")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

individual_data <- fread(file.path(input_dir, "Figure1_individual_ancestry_K5.csv"))
mean_ancestry <- fread(file.path(input_dir, "Figure1_mean_ancestry_by_continent.csv"))
assignment_data <- fread(file.path(input_dir, "Figure1_assignment_proportions_Q0.70.csv"))
cv_data <- fread(file.path(input_dir, "ADMIXTURE_CV_K1_K6.csv"))

stopifnot(uniqueN(individual_data$IID) == 303L)
stopifnot(all(c("A", "B", "D") %in% individual_data$subgenome))

ancestry_colors <- c(
  Q1 = "#C94A44", Q2 = "#3A8F61", Q3 = "#477DB3",
  Q4 = "#8A62A8", Q5 = "#E59A35"
)
assignment_colors <- c(ancestry_colors, Admixed = "#B8B8B8")
continent_levels <- c("AF", "AS", "EU", "NAM", "OC", "SAM", "N/A")

individual_data[, `:=`(
  subgenome = factor(subgenome, levels = c("A", "B", "D")),
  Continent_short = factor(Continent_short, levels = continent_levels)
)]
individual_long <- individual_data %>%
  pivot_longer(starts_with("Q"), names_to = "Ancestry", values_to = "Fraction") %>%
  mutate(Ancestry = factor(Ancestry, levels = names(ancestry_colors)))

mean_ancestry <- mean_ancestry %>%
  mutate(
    subgenome = factor(subgenome, levels = c("A", "B", "D")),
    Continent_short = factor(Continent_short, levels = continent_levels),
    Ancestry = factor(Ancestry, levels = names(ancestry_colors))
  )

assignment_data <- assignment_data %>%
  mutate(
    subgenome = factor(subgenome, levels = c("A", "B", "D")),
    Continent_short = factor(Continent_short, levels = continent_levels),
    Assignment = factor(Assignment, levels = c(names(ancestry_colors), "Admixed"))
  )

theme_population <- theme_classic(base_size = 6.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.text = element_text(size = 5.8, colour = "black"),
    axis.title = element_text(size = 6.5, colour = "black"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5.5),
    strip.background = element_rect(fill = "#F0F0F0", colour = "#B5B5B5", linewidth = 0.25),
    strip.text = element_text(size = 5.8, face = "bold"),
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.004, 0.985),
    plot.margin = margin(2, 2, 2, 2)
  )

make_admixture_panel <- function(subgenome_name, tag, show_legend = FALSE) {
  ggplot(filter(individual_long, subgenome == subgenome_name),
         aes(within_continent_index, Fraction, fill = Ancestry)) +
    geom_col(width = 1, colour = NA) +
    facet_grid(~ Continent_short, scales = "free_x", space = "free_x", drop = FALSE) +
    scale_fill_manual(values = ancestry_colors, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = expansion(mult = 0)) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    labs(x = NULL, y = "Ancestry proportion", fill = NULL, tag = tag) +
    theme_population +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.key.height = grid::unit(3, "mm"),
      legend.key.width = grid::unit(3, "mm"),
      panel.spacing.x = grid::unit(0.7, "mm")
    )
}

p_a <- make_admixture_panel("A", "A", TRUE)
p_b <- make_admixture_panel("B", "B")
p_c <- make_admixture_panel("D", "C")

p_d <- ggplot(mean_ancestry, aes(Ancestry, Continent_short, fill = Mean_Q)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = percent(Mean_Q, accuracy = 1)), size = 1.65,
            family = "Arial", colour = "#222222") +
  facet_grid(~ subgenome, drop = FALSE) +
  scale_fill_gradient(low = "#F4F6F8", high = "#2F6F9F", limits = c(0, 1), name = "Mean Q") +
  scale_y_discrete(limits = rev(continent_levels), drop = FALSE) +
  labs(x = "Ancestry component", y = "Continent", tag = "D") +
  theme_population +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    panel.border = element_rect(colour = "#B5B5B5", fill = NA, linewidth = 0.25),
    legend.position = "bottom", legend.direction = "horizontal"
  ) +
  guides(fill = guide_colourbar(title.position = "top",
                                barwidth = grid::unit(24, "mm"),
                                barheight = grid::unit(2.2, "mm")))

p_e <- ggplot(assignment_data, aes(Continent_short, Proportion, fill = Assignment)) +
  geom_col(width = 0.78, colour = NA) +
  facet_grid(~ subgenome, drop = FALSE) +
  scale_fill_manual(values = assignment_colors, drop = FALSE, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = 0)) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  labs(x = "Continent", y = "Assigned proportion", tag = "E") +
  theme_population +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.width = grid::unit(2.5, "mm")
  )

top_panels <- wrap_plots(list(p_a, p_b, p_c), ncol = 1)
bottom_panels <- wrap_plots(list(p_d, p_e), ncol = 2, widths = c(1.15, 1))
figure <- wrap_plots(list(top_panels, bottom_panels), ncol = 1, heights = c(2.65, 1.35)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

continent_note <- paste0(
  "Continent abbreviations: AF, Africa; AS, Asia; EU, Europe; ",
  "NAM, North America; OC, Oceania; SAM, South America; and N/A, unknown origin."
)
figure <- figure + plot_annotation(
  caption = continent_note,
  theme = theme(plot.caption = element_text(family = "Arial", size = 5.5,
                                             colour = "#222222", hjust = 0,
                                             margin = margin(t = 3)))
)

base <- file.path(output_dir, "Figure1_population_structure")
ggsave(paste0(base, ".png"), figure, width = 183, height = 135,
       units = "mm", dpi = 300, device = agg_png, bg = "white")
ggsave(paste0(base, ".pdf"), figure, width = 183, height = 135,
       units = "mm", device = cairo_pdf, bg = "white")

p_cv <- ggplot(cv_data, aes(K, CV_error, colour = subgenome)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.7) +
  geom_vline(xintercept = 5, linetype = "dashed", linewidth = 0.35, colour = "#666666") +
  scale_x_continuous(breaks = 1:6) +
  scale_colour_manual(values = c(A = "#C94A44", B = "#3A8F61", D = "#477DB3")) +
  labs(x = "Number of ancestral populations (K)", y = "10-fold CV error", colour = "Subgenome") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "Figure1_ADMIXTURE_CV.png"), p_cv,
       width = 100, height = 70, units = "mm", dpi = 300,
       device = agg_png, bg = "white")

message("Saved Figure 1 and ADMIXTURE CV plot to ", output_dir)
