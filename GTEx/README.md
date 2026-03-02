# GTEx

Real-data application of MRGN, MRPC, and GMAC to GTEx V8 gene-expression trios. The primary analysis focuses on **WholeBlood** tissue, using trios restricted to protein-coding and lncRNA genes.

## Overview

The analysis pipeline:
1. Filters GTEx trios to retain only protein-coding and lncRNA gene pairs (using BioMart annotations)
2. Selects significant PCs as confounders via correlation testing with FDR correction
3. Runs MRGN, MRPC, and GMAC inference on each trio with its matched confounders
4. Assembles master result tables comparing all three methods

## Subfolders

### `data/`

Input data and intermediate products for the GTEx analysis.

| File | Description |
|---|---|
| `all.data.unqiue.snps.pclrna.only.WholeBlood.RData` | Filtered trio data matrix (n × 3p) containing only protein-coding and lncRNA trios for WholeBlood. |
| `data.snp.cis.trans.final.WholeBlood.V8.unique.snps.RData` | Original n × 3p trio matrix from GTEx V8 with unique SNPs for WholeBlood. |
| `data.with.PCs.WholeBlood.RData` | List of data frames (one per trio), each containing trio data plus its matched PC confounders. |
| `PCs.matrix.WholeBlood.RData` | PC score matrix (n × m) for WholeBlood tissue. |
| `List.significant.asso1.WholeBlood.RData` | List of length m (one per PC) containing column indices of trio variables significantly correlated with each PC. |
| `List.Match.significant.trios.WholeBlood.RData` | List of length p_sub where each element contains indices of PCs correlated with that trio. |
| `kclist_top5_tiss.RData` | Known clinical covariate lists for the top-5 tissues by sample size. |
| `WholeBloodmaster_table.RData` | Master results table combining MRGN, MRPC, and GMAC inference for all WholeBlood trios. |
| `PC_LRNA_PC_Selection_manu.R` | Script that filters trios to protein-coding/lncRNA genes (via BioMart), performs PC-based confounder selection, extracts known covariates (pcr, platform, sex) from GTEx V8 covariate files, and saves the filtered datasets. |

### `scripts/`

| File | Description |
|---|---|
| `make_mrgn_triotables.R` | Builds comprehensive master trio result tables by merging MRGN, MRPC-ADDIS, and GMAC results. Uses BioMart gene annotations to classify genes, identifies M1 subtypes, and combines results from multiple confounder selection strategies. |
| `make_GMAC_tables.R` | Assembles GMAC results tables for top-5 tissues by matching GMAC cis/trans results to MRGN trio tables. Translates significance into mediation categories. |
| `use_liberal_conf_selection_noFDR.R` | Applies a liberal confounder selection procedure (no FDR correction, alpha < 0.01 and alpha < 0.05 cutoffs) across the top-5 tissues. Runs `get.conf.trios()` and `infer.trio()` with permutation testing and saves all outputs. |

### `results/`

Inference outputs organized by method:

| Subfolder | Contents |
|---|---|
| `MRGN/` | MRGN inference results for WholeBlood, including `gtex.analysis.mrgn.R` (analysis script), inferred models, and regression-based confounder selection results. |
| `GMAC/` | GMAC cis/trans inference results for WholeBlood, including `gtex.analysis.gmac.R` (analysis script) and GMAC input/output lists. |
| `MRPC/` | MRPC inference results for WholeBlood. |

### `trios_data_prlnc_liberal_conf_sel/`

MRGN results using the liberal confounder selection (no FDR correction) for top-5 tissues.

| Subfolder | Description |
|---|---|
| `data_with_confs/` | List of data frames for each trio with its selected confounders (`.RData`). |
| `infer_trio_results/` | Full `infer.trio()` output for each trio: a 14 × p matrix saved as `.RData`. |
| `list_output/` | Raw output from `get.conf.trios()` saved as `.RData`. |

### `Prob_trio_plots/`

Diagnostic PDF plots for specific WholeBlood trios, including probability plots and comparisons with/without the permutation test.

### `trio_tables/`

Placeholder for assembled trio result tables (currently empty).

## PC Selection Procedure

1. Compute Pearson correlations between each PC and all trio variables using `psych::corr.test(adjust="none")`
2. Apply q-value FDR correction (10%) to each PC's p-values
3. Match significant PCs to individual trios based on shared variable associations
4. Construct trio-specific data frames incorporating matched PCs as confounders
