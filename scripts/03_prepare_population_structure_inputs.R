#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

project_dir <- normalizePath(".", mustWork = TRUE)
input_hapmap <- file.path(
  project_dir,
  "data/genotypes/YoGI_genotypes_MAF0.05_303.hmp.txt.gz"
)
maf_summary_file <- file.path(
  project_dir,
  "data/genotypes/YoGI_MAF_summary_179253_SNPs.csv.gz"
)
output_dir <- file.path(
  project_dir,
  "results/population_structure/input"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hapmap_connection <- gzfile(input_hapmap, open = "rt")
hapmap <- read.delim(
  hapmap_connection,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)
close(hapmap_connection)
maf_connection <- gzfile(maf_summary_file, open = "rt")
maf_summary <- read.csv(maf_connection, check.names = FALSE)
close(maf_connection)

metadata_columns <- c(
  "rs", "alleles", "Chromosome", "pos", "strand", "assembly",
  "center", "protLSID", "assayLSID", "panel", "QCcode"
)
stopifnot(identical(names(hapmap)[seq_along(metadata_columns)], metadata_columns))
sample_ids <- names(hapmap)[-(seq_along(metadata_columns))]
stopifnot(length(sample_ids) == 303L, !anyDuplicated(sample_ids))
stopifnot(nrow(hapmap) == 32221L, !anyDuplicated(hapmap$rs))

summary_index <- match(hapmap$rs, maf_summary$SNP)
stopifnot(!anyNA(summary_index))
matched_summary <- maf_summary[summary_index, ]
stopifnot(all(matched_summary$keep_MAF), all(matched_summary$MAF >= 0.05))
stopifnot(all(matched_summary$Chromosome == hapmap$Chromosome))

hapmap$alleles <- paste(
  matched_summary$major_allele,
  matched_summary$minor_allele,
  sep = "/"
)
stopifnot(!anyNA(hapmap$alleles), !any(grepl("NA", hapmap$alleles, fixed = TRUE)))

subgenomes <- list(
  A = list(chromosomes = 1:7, prefix = "YoGI_MAF0.05_A_C17"),
  B = list(chromosomes = 8:14, prefix = "YoGI_MAF0.05_B_C814"),
  D = list(chromosomes = 15:21, prefix = "YoGI_MAF0.05_D_C1521")
)

manifest <- do.call(
  rbind,
  lapply(names(subgenomes), function(subgenome) {
    config <- subgenomes[[subgenome]]
    subset_data <- hapmap[hapmap$Chromosome %in% config$chromosomes, , drop = FALSE]
    output_file <- file.path(output_dir, paste0(config$prefix, ".hmp.txt"))
    write.table(
      subset_data,
      output_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      na = "NA"
    )

    data.frame(
      subgenome = subgenome,
      chromosomes = paste(range(config$chromosomes), collapse = "-"),
      samples = length(sample_ids),
      snps = nrow(subset_data),
      minimum_maf = min(matched_summary$MAF[hapmap$Chromosome %in% config$chromosomes]),
      maximum_maf = max(matched_summary$MAF[hapmap$Chromosome %in% config$chromosomes]),
      hapmap_file = output_file
    )
  })
)

stopifnot(sum(manifest$snps) == 32221L)
write.csv(
  manifest,
  file.path(output_dir, "population_structure_MAF0.05_input_manifest.csv"),
  row.names = FALSE
)
writeLines(sample_ids, file.path(output_dir, "sample_order_303.txt"))

print(manifest[, c("subgenome", "chromosomes", "samples", "snps", "minimum_maf", "maximum_maf")])
