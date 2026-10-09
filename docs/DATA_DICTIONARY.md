# Data dictionary

## Phenotypes

`data/phenotypes/YoGI_traits_303.tsv`

One row per accession. `Taxa` is the YoGI accession identifier. Trait columns use the analysis names: `Plant_height_node`, `Peduncle_length`, `Leaf_length`, `Spike_number`, `Spike_length`, `Spikelet_number`, `Plant_height_spike`, and `Awn`. `Awn_check` is retained for provenance but is not treated as an additional trait in the manuscript.

`data/phenotypes/Figure2_source_data_303.csv`

The eight traits with publication labels, used directly to generate Figure 2.

## Genotypes

`data/genotypes/YoGI_genotypes_MAF0.05_303.hmp.txt.gz`

The exact HapMap-format genotype matrix used for population structure and GWAS: 32,221 SNP rows and 303 accession columns after applying MAF >= 0.05.

`data/genotypes/YoGI_MAF_summary_179253_SNPs.csv.gz`

One row per input SNP. Main fields are SNP identifier, chromosome, position, MAF, call rate, heterozygosity, allele counts, major allele, minor allele, and the Boolean `keep_MAF` decision.

## Population structure

`data/population_structure/Figure1_individual_ancestry_K5.csv`

K = 5 ADMIXTURE ancestry coefficients for each accession and each subgenome, after component relabelling for consistent display.

`Figure1_mean_ancestry_by_continent.csv` and `Figure1_assignment_proportions_Q0.70.csv`

Values used in Figure 1D and 1E. Accessions with maximum Q < 0.70 are classified as admixed.

`ADMIXTURE_CV_K1_K6.csv`

Ten-fold cross-validation errors for K = 1-6 in the A, B, and D subgenomes.

## GWAS

`data/gwas/*_GAPIT_MLM.csv.gz`

GAPIT MLM result tables for the eight manuscript traits. Main fields are SNP identifier, chromosome, position, P value, MAF, number of observations, and estimated effect.

`data/gwas/GWAS_summary_by_trait.csv`

Marker count, strongest SNP, minimum P value, MAF, number of significant SNPs, and genomic inflation factor for each trait.

## Associative transcriptomics

`data/associative_transcriptomics/*_AT_results.tsv.gz`

Regression result tables for 46,248 expression markers for each of the eight traits. `Marker` is the transcript identifier, `log10PA` is the plotted adjusted -log10(P) value, `Graph` encodes chromosome order, and `Sort` is the cumulative plotting position.

## Supplementary tables

`Table_S1_significant_awn_SNPs.csv`

The four genome-wide significant SNPs for awn presence.

`Table_S2_top3_nonsignificant_SNPs.csv`

The three strongest non-significant SNPs for each of the other seven traits.
