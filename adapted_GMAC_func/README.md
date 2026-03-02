# adapted_GMAC_func

Modified GMAC (Genomic Mediation Analysis with Adaptive Confounding) source code with instrumentation for timing confounder selection and inference steps. These timing measurements are used to produce the computation time comparison figures in the MRGN manuscript.

## Files

| File | Description |
|---|---|
| `gmac_get_conf.R` | Implements the GMAC confounder selection procedure. The `gmac_get_conf()` function takes a covariate pool and trio data, uses stratified FDR to identify covariates that are not children/intermediates, and estimates which pool covariates act as confounders for each trio. Supports parallel execution via an optional cluster argument. |
| `gmac_one_trio.R` | Applies the GMAC mediation test to a single trio (SNP, mediator, outcome) with its associated confounders. Computes the indirect effect t-statistic, performs a genotype-stratified permutation test, and returns the p-value and beta-change (proportion of total effect mediated). Defines helper functions `Indirect()`, `nominal.pfun()`, and `get.beta.change()`. |
| `bibrefs.bib` | BibTeX references related to the GMAC method. |

## GMAC_moded/

A lightly modified R package version of GMAC. This is a self-contained R package with the standard structure (`R/`, `man/`, `data/`, `NAMESPACE`, `DESCRIPTION`). The modifications enable timing of internal steps for benchmarking against MRGN.
