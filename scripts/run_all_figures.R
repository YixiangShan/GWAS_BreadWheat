#!/usr/bin/env Rscript

scripts <- c(
  "scripts/04_plot_figure1.R",
  "scripts/05_plot_figure2.R",
  "scripts/06_plot_figures3_5_S1_S2.R",
  "scripts/07_plot_figureS3_AT.R"
)

rscript <- file.path(R.home("bin"), "Rscript")
for (script in scripts) {
  message("Running ", script)
  status <- system2(rscript, script)
  if (!identical(status, 0L)) {
    stop("Script failed: ", script)
  }
}

message("All manuscript figures were generated successfully.")
