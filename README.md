# GWAS of morphological traits in bread wheat

Data and R scripts supporting GWAS and associative transcriptomics analyses of eight morphological traits in 303 bread wheat accessions.

This repository accompanies the manuscript **"A genome-wide association study (GWAS) to identify genetic factors influencing wheat yield-related morphological traits"**. It contains the source data and plotting code for Figures 1-5, Supplementary Figures S1-S3, and Supplementary Tables S1-S2.

## Study overview

- Panel: 303 accessions from the YoGI bread wheat diversity panel.
- Phenotypes: awn presence, plant height without spike, plant height with spike, peduncle length, leaf length, spikelet number, spike length, and spike number.
- Genotypes: 179,253 transcriptome-derived SNPs before filtering; 32,221 SNPs retained at MAF >= 0.05.
- GWAS: GAPIT mixed linear model with three genotype principal components and a kinship matrix.
- Population structure: ADMIXTURE analyses of the A, B, and D subgenomes; K = 5 is displayed in Figure 1.
- Associative transcriptomics: 46,248 expression markers tested for each trait.

Only awn presence had genome-wide significant SNP associations after MAF filtering. Four significant SNPs mapped to the known *B1 (Tipped 1)* region on chromosome 5A.

## Repository contents

```text
data/
  phenotypes/                    Phenotype table and Figure 2 source data
  genotypes/                     MAF-filtered HapMap matrix and MAF summary
  population_structure/          Figure 1 source tables and CV values
  gwas/                          GAPIT result tables for eight traits
  associative_transcriptomics/   AT result tables for eight traits
  supplementary_tables/          Source and final Tables S1-S2
figures/                          Manuscript figure previews
scripts/                          MAF, GAPIT, table, and plotting scripts
environment/                      R package and session information
docs/                             Data dictionary and reproducibility notes
```

Large text tables are compressed with gzip. The included scripts read them with `R.utils` when available and otherwise use the system `gzip` command.

## Reproduce the figures

Run commands from the repository root:

```bash
Rscript scripts/00_install_packages.R
Rscript scripts/run_all_figures.R
Rscript scripts/08_prepare_supplementary_tables.R
Rscript scripts/09_build_manifest.R
```

Generated files are written to `figures/generated/`. The archived figure previews used in the manuscript are retained in `figures/`.

To rerun the GAPIT analysis from the exact 32,221-SNP matrix:

```bash
Rscript scripts/02_run_gwas_gapit.R
```

GAPIT output is written to `results/gwas/`. This analysis may take substantially longer than redrawing the figures from the archived GWAS result tables.

## Original sequencing data

The RNA-sequencing reads from which SNP genotypes and expression values were derived are available through NCBI BioProject [PRJNA912645](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA912645). The source study is Barratt *et al.* (2023), [doi:10.1111/tpj.16248](https://doi.org/10.1111/tpj.16248).

The large unfiltered genotype matrix is not included. The repository provides the complete marker-level MAF summary for all 179,253 input SNPs and the exact 32,221-SNP matrix used for the MAF-filtered analyses. The processed AT result tables needed to reproduce Supplementary Figure S3 are also included; the large RPKM matrix is not duplicated here.

## Citation

Please cite the associated manuscript and this repository. Citation metadata are provided in `CITATION.cff`.

## Author

Yixiang Shan

## License

The R code is released under the MIT License. See `DATA_LICENSE.md` for the data files and third-party source data.
