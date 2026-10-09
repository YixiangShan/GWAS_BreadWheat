#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- normalizePath(".", mustWork = TRUE)
geno_file <- file.path(project_dir, "data/genotypes/YoGI_genotypes_MAF0.05_303.hmp.txt.gz")
trait_file <- file.path(project_dir, "data/phenotypes/YoGI_traits_303.tsv")
out_dir <- file.path(project_dir, "results/gwas")
gapit_local <- file.path(project_dir, "scripts/vendor/gapit_functions.txt")
source(file.path(project_dir, "scripts/lib/read_compressed.R"))

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

if (!exists("GAPIT", mode = "function")) {
  if (file.exists(gapit_local)) {
    source(gapit_local)
  } else if (requireNamespace("GAPIT3", quietly = TRUE)) {
    suppressPackageStartupMessages(library(GAPIT3))
  } else {
    source("https://raw.githubusercontent.com/jiabowang/GAPIT/refs/heads/master/gapit_functions.txt")
  }
}

traits <- fread(trait_file, data.table = FALSE)
geno <- read_compressed_table(
  geno_file,
  header = FALSE,
  data.table = FALSE,
  check.names = FALSE
)
trait_names <- c(
  "Awn", "Plant_height_node", "Plant_height_spike", "Spike_number",
  "Spike_length", "Spikelet_number", "Peduncle_length", "Leaf_length"
)
stopifnot(all(trait_names %in% names(traits)))

old_wd <- getwd()
setwd(out_dir)
on.exit(setwd(old_wd), add = TRUE)

for (trait in trait_names) {
  message("Running MLM GWAS with filtered-SNP PCA for: ", trait)
  y <- traits[, c("Taxa", trait)]
  colnames(y) <- c("Taxa", trait)
  GAPIT(
    Y = y,
    G = geno,
    PCA.total = 3,
    model = "MLM",
    Multiple_analysis = FALSE
  )
}

message("GWAS output directory: ", out_dir)
