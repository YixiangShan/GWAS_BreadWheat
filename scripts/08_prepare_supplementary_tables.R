#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts/lib/read_compressed.R"))
input_dir <- file.path(root, "data/gwas")
output_dir <- file.path(root, "data/supplementary_tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

traits <- c(
  "Awn", "Plant_height_node", "Plant_height_spike", "Spike_number",
  "Spike_length", "Spikelet_number", "Peduncle_length", "Leaf_length"
)

gwas <- rbindlist(lapply(traits, function(trait) {
  x <- read_compressed_table(file.path(input_dir, paste0(trait, "_GAPIT_MLM.csv.gz")))
  setnames(x, c("Chr", "Pos", "P.value"), c("Chromosome", "Position", "P"))
  x[, Trait := trait]
  x
}), use.names = TRUE, fill = TRUE)

threshold <- 0.05 / uniqueN(gwas[Trait == "Awn", SNP])
columns <- c("Trait", "SNP", "Chromosome", "Position", "P", "MAF", "nobs", "Effect")
table_s1 <- gwas[Trait == "Awn" & P < threshold, ..columns][order(P)]
table_s2 <- gwas[Trait != "Awn", head(.SD[order(P)], 3L), by = Trait, .SDcols = columns[-1L]]

stopifnot(nrow(table_s1) == 4L, nrow(table_s2) == 21L)
fwrite(table_s1, file.path(output_dir, "Table_S1_significant_awn_SNPs.csv"))
fwrite(table_s2, file.path(output_dir, "Table_S2_top3_nonsignificant_SNPs.csv"))

message("Prepared Tables S1 and S2.")
