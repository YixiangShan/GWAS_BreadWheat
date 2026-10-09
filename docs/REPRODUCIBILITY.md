# Reproducibility notes

## Analysis set

All analyses in the current manuscript use the same 303 accessions. YoGI_342 was excluded because genotype data were unavailable. SNP MAF was calculated within these 303 accessions. Markers with MAF < 0.05 were removed, leaving 32,221 SNPs.

## GWAS

GWAS was run in GAPIT using a mixed linear model. Three principal components recalculated from the MAF-filtered genotype matrix and a kinship matrix were included to control population structure and relatedness. The Bonferroni threshold is 0.05 / 32,221 = 1.55 x 10^-6.

`scripts/02_run_gwas_gapit.R` reruns this analysis. The vendored GAPIT functions used during the project are retained in `scripts/vendor/`.

## Population structure

ADMIXTURE v1.3.0 was run separately for the A, B, and D subgenomes with K = 1-6, 10-fold cross-validation, and random seed 20260901. Figure 1 displays K = 5 as a common resolution. B- and D-subgenome component labels were permuted to maximize sample-wise correlation with the A-subgenome components. This relabelling affects display labels only.

The processed K = 5 coefficients and CV values are provided so Figure 1 can be redrawn without installing ADMIXTURE or PLINK.

## Associative transcriptomics

Each of 46,248 normalized expression markers was tested separately. Q1 and Q2 ancestry coefficients were used as covariates. Supplementary Figure S3 is reproduced from the archived regression result tables. The large RPKM matrix is not duplicated in this repository; the underlying RNA-sequencing reads are available under NCBI BioProject PRJNA912645.

## Expected outputs

- `scripts/04_plot_figure1.R`: Figure 1 and the ADMIXTURE CV plot.
- `scripts/05_plot_figure2.R`: Figure 2 and its source summaries.
- `scripts/06_plot_figures3_5_S1_S2.R`: Figures 3-5, S1, and S2.
- `scripts/07_plot_figureS3_AT.R`: Figure S3.
- `scripts/08_prepare_supplementary_tables.R`: Tables S1 and S2.

Run all scripts from the repository root so relative paths resolve correctly.
