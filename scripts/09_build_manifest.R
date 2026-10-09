#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(".", mustWork = TRUE)
files <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE)
files <- files[file.info(files)$isdir %in% FALSE]
relative <- substring(files, nchar(root) + 2L)

keep <- !grepl("^([.]git/|results/|figures/generated/)", relative) &
  relative != "docs/FILE_MANIFEST.tsv"
files <- files[keep]
relative <- relative[keep]

sha256 <- vapply(files, function(path) {
  output <- system2("shasum", c("-a", "256", path), stdout = TRUE)
  sub("[[:space:]].*$", "", output[[1]])
}, character(1))

manifest <- data.table(
  path = relative,
  bytes = as.numeric(file.info(files)$size),
  sha256 = sha256
)
setorder(manifest, path)
fwrite(manifest, file.path(root, "docs/FILE_MANIFEST.tsv"), sep = "\t")

message("Wrote ", nrow(manifest), " entries to docs/FILE_MANIFEST.tsv")
