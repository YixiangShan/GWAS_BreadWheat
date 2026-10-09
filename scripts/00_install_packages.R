#!/usr/bin/env Rscript

cran_packages <- c(
  "data.table", "dplyr", "GGally", "ggplot2", "ggrepel",
  "patchwork", "ragg", "R.utils", "scales", "tidyr"
)

missing <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

message("Required plotting packages are installed.")
