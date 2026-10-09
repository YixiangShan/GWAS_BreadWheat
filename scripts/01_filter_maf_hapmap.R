#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("scripts/lib/read_compressed.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(
    paste(
      "Usage: Rscript scripts/01_filter_maf_hapmap.R <unfiltered_hapmap>",
      "[filtered_hapmap] [maf_summary] [threshold]",
      "\nThe large unfiltered matrix is not included in this repository."
    )
  )
}

input_file <- args[[1]]
output_file <- if (length(args) >= 2) args[[2]] else "data/genotypes/YoGI_genotypes_MAF0.05_303.hmp.txt"
summary_file <- if (length(args) >= 3) args[[3]] else sub("\\.txt$", "_MAF_summary.csv", output_file)
threshold <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.05

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}
if (!is.finite(threshold) || threshold < 0 || threshold > 0.5) {
  stop("MAF threshold must be between 0 and 0.5")
}

iupac_counts <- list(
  A = c(A = 2, C = 0, G = 0, T = 0),
  C = c(A = 0, C = 2, G = 0, T = 0),
  G = c(A = 0, C = 0, G = 2, T = 0),
  T = c(A = 0, C = 0, G = 0, T = 2),
  R = c(A = 1, C = 0, G = 1, T = 0),
  Y = c(A = 0, C = 1, G = 0, T = 1),
  S = c(A = 0, C = 1, G = 1, T = 0),
  W = c(A = 1, C = 0, G = 0, T = 1),
  K = c(A = 0, C = 0, G = 1, T = 1),
  M = c(A = 1, C = 1, G = 0, T = 0),
  B = c(A = 0, C = 1, G = 1, T = 1),
  D = c(A = 1, C = 0, G = 1, T = 1),
  H = c(A = 1, C = 1, G = 0, T = 1),
  V = c(A = 1, C = 1, G = 1, T = 0)
)

parse_genotype_count <- function(code) {
  code <- toupper(trimws(as.character(code)))
  if (is.na(code) || code == "" || code %in% c("N", "NN", "NA", "-", "--", "?")) {
    return(c(A = 0, C = 0, G = 0, T = 0))
  }
  if (!is.null(iupac_counts[[code]])) {
    return(iupac_counts[[code]])
  }
  alleles <- unlist(strsplit(gsub("[/|]", "", code), "", fixed = FALSE), use.names = FALSE)
  alleles <- alleles[alleles %in% c("A", "C", "G", "T")]
  counts <- c(A = 0, C = 0, G = 0, T = 0)
  if (length(alleles) > 0) {
    tab <- table(alleles)
    counts[names(tab)] <- as.numeric(tab)
  }
  counts
}

summarise_marker <- function(genotypes) {
  tab <- table(toupper(trimws(as.character(genotypes))), useNA = "no")
  tab <- tab[!names(tab) %in% c("", "N", "NN", "NA", "-", "--", "?")]

  allele_counts <- c(A = 0, C = 0, G = 0, T = 0)
  heterozygote_calls <- 0
  called_genotypes <- 0

  if (length(tab) > 0) {
    for (code in names(tab)) {
      counts <- parse_genotype_count(code)
      allele_counts <- allele_counts + counts * as.numeric(tab[[code]])
      if (sum(counts > 0) > 1) {
        heterozygote_calls <- heterozygote_calls + as.numeric(tab[[code]])
      }
      called_genotypes <- called_genotypes + as.numeric(tab[[code]])
    }
  }

  observed <- allele_counts[allele_counts > 0]
  total_alleles <- sum(observed)
  call_rate <- called_genotypes / length(genotypes)

  if (length(observed) < 2 || total_alleles == 0) {
    return(c(
      MAF = 0,
      call_rate = call_rate,
      heterozygosity = ifelse(called_genotypes > 0, heterozygote_calls / called_genotypes, NA_real_),
      n_called = called_genotypes,
      major_allele = ifelse(length(observed) == 1, names(observed)[1], NA),
      minor_allele = NA
    ))
  }

  freqs <- sort(observed / total_alleles, decreasing = TRUE)
  c(
    MAF = unname(tail(freqs, 1)),
    call_rate = call_rate,
    heterozygosity = heterozygote_calls / called_genotypes,
    n_called = called_genotypes,
    major_allele = names(freqs)[1],
    minor_allele = names(tail(freqs, 1))
  )
}

message("Reading: ", input_file)
geno <- read_compressed_table(
  input_file,
  data.table = TRUE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
if (ncol(geno) <= 11) {
  stop("Expected HapMap/GAPIT format with 11 metadata columns followed by sample genotype columns.")
}

sample_cols <- 12:ncol(geno)
message("Markers: ", nrow(geno))
message("Samples: ", length(sample_cols))
message("Calculating MAF...")

maf_mat <- t(apply(as.matrix(geno[, ..sample_cols]), 1, summarise_marker))
maf_dt <- as.data.table(maf_mat)
numeric_cols <- c("MAF", "call_rate", "heterozygosity", "n_called")
maf_dt[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]

marker_id_col <- names(geno)[1]
chrom_col <- names(geno)[3]
pos_col <- names(geno)[4]
summary_dt <- data.table(
  SNP = geno[[marker_id_col]],
  Chromosome = geno[[chrom_col]],
  Position = geno[[pos_col]]
)
summary_dt <- cbind(summary_dt, maf_dt)
summary_dt[, keep_MAF := MAF >= threshold]

kept <- summary_dt[["keep_MAF"]]
filtered <- geno[kept]

message("Keeping ", nrow(filtered), " of ", nrow(geno), " markers at MAF >= ", threshold)
message("Writing filtered genotype: ", output_file)
fwrite(filtered, output_file, sep = "\t", quote = FALSE, na = "NA")

message("Writing MAF summary: ", summary_file)
fwrite(summary_dt, summary_file)

message("Done.")
