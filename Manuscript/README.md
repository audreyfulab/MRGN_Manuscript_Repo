# Manuscript

All materials for the MRGN manuscript: LaTeX source, figure/table generation scripts, final outputs, supplementary materials, and helper functions.

## Subfolders

### `scripts/`

R scripts for generating all manuscript figures, tables, and analyses.

| File | Description |
|---|---|
| `6_6_2023_all_conf_types_simulation_updated_2025.Rmd` | Comprehensive R Markdown simulation report ("MRGN Simulation Set 3"). Documents the simulation design, parameters, confounder types, selection/filtering procedures, and defines all four performance metrics (confounder selection precision/recall, class-based, edge-based, T1–T2 edge-based). |
| `loadedResults.R` | Central data-loading script that reads all simulation results (MRGN, MRPC, GMAC) for both standard and many-confounder scenarios, computation times, confounder lists, and dataset parameters. Used by all figure/table scripts. |
| `create_main_figs.R` | Generates the main manuscript simulation figures. Loads results via `loadedResults.R`, computes T1–T2 edge-based metrics, and produces loess-smoothed performance curves vs. simulation parameters (residual SD, MAF, SNP signal, mediation signal, number of confounders). |
| `create_main_figs_helpers.R` | Utility functions for `create_main_figs.R`: model label conversion, T1–T2 edge extraction from adjacency matrices, precision/recall/F1 computation, edge-based scoring, and the main `plot.sim.metrics()` plotting function. |
| `create_GTEx_figs.R` | Creates GTEx-specific figures: bar charts of inferred model frequencies (M0–M4), cis/trans mediation breakdowns, and permutation test effect comparisons for WholeBlood. Outputs supplementary tables and PDF figures. |
| `create_time_conf_figs.R` | Generates supplementary figures comparing computation times for confounder selection between GMAC and MRGN across varying trio counts and covariate pool sizes (100, 500, 1000). |
| `create_model_misspecification_figs.R` | Creates model misspecification analysis figures. Compares true nominal p-values against median simulated parametric p-values for three misspecification scenarios (STMS, LTMS, SIMS), stratified by MAF. |
| `compute_true_upper_bound.R` | Computes MRGN's "true upper bound" performance by running `infer.trio()` with ground-truth confounders (bypassing confounder selection). Calculates class-based, edge-based, and T1–T2 edge-based metrics. |
| `updated_main_results_table2025.R` | Generates the main manuscript results tables by computing all performance metrics across simulation scenarios. Assembles a master table combining MRGN, MRPC, and GMAC inferences with computation times. |
| `check_mrpc_comp_times.R` | Analyzes MRPC-ADDIS computation times. Plots time vs. number of confounders, identifies problematic trios (>4 hours), and benchmarks against MRGN/GMAC. |
| `gmac_valid_res_6_27_2023.R` | Loads GMAC validation results and creates supplementary figure SF9: GMAC Type I error rate vs. number of confounders at 0.01 and 0.05 significance cutoffs. |
| `Plot_dist_rare_genotype_trios.R` | Generates supplementary figures showing gene expression distributions for GTEx WholeBlood trios with rare MAF (<10%). Includes scatter plots and histograms with skewness/kurtosis statistics. |
| `plot_SF_perm_evidence.R` | Creates figures visualizing the permutation test's effect on MRGN inference. Classifies model changes as "liberal" or "conservative" and shows impact on rare-allele trios. Also exports GTEx result tables. |

**Data files in `scripts/`:**

| File | Description |
|---|---|
| `loadedResults.R` | (see above — also serves as a data pipeline) |
| `mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types.RData` | Combined simulation results for all methods (5K dataset). |
| `mrpc_5k_conf_list_all_confs_filtered_all_mods.RData` | MRPC confounder lists for filtered simulations. |
| `mrpc_v_mrgn_v_gmac_5k_params_all_mods_conf_types.RData` | Simulation parameters for all methods. |
| `all_trios_output_cis.Rdata` | Cis-mediation output for all trios. |

### `figures/`

Final main manuscript figures (PDF):

| File | Description |
|---|---|
| `5 causal models2.pdf` | Diagram of the 5 causal models (M0–M4). |
| `MF1_MRGN.GMAC.MRPC.all.params.t1t2.pdf` | Main Figure 1: T1–T2 edge-based performance across all simulation parameters. |
| `MF3_GTEx.model.and.t1t2.edge.bargraphs.pdf` | Main Figure 3: GTEx model frequency and T1–T2 edge bar graphs. |

### `tables/`

Main manuscript tables:

| File | Description |
|---|---|
| `MT1.xlsx` | Main Table 1. |
| `MT2_Combined_ALL_METICS_15confSIMS.csv` | Main Table 2: All performance metrics for the standard (≤15 confounders) simulations. |
| `MT3_Combined_ALL_METICS_MANYconfSIMS.csv` | Main Table 3: All performance metrics for the many-confounder (15–50) simulations. |

### `supplementary_figures/`

Supplementary figures (PDF):

| File | Description |
|---|---|
| `Confounding Vars2.pdf` | Diagram of confounding variable types. |
| `SF1_comp_times_for_conf_selection_updated.pdf` | SF1: Computation time comparison for confounder selection. |
| `SF3_GTEx_Permutation_and_Model_Changes.pdf` | SF3: GTEx permutation test effect on model assignments. |
| `SF4_Confounder_Selection_Performance.pdf` | SF4: Confounder selection precision/recall. |
| `SF5_Simulation_Permutation_and_Model_Changes.pdf` | SF5: Simulation permutation test and model changes. |
| `SF6_Model_Misspecification.pdf` | SF6: Model misspecification analysis. |
| `SF7_Confounder_Selection_Impact_On_Inference.pdf` | SF7: Impact of confounder selection on inference. |
| `SF8_Trio_Skew_Stat_Distribtuions.pdf` | SF8: Skewness/kurtosis distributions in trio expression data. |
| `SF9_GMAC_Valid_TypeI_error.pdf` | SF9: GMAC Type I error validation. |
| `SF_MRGN.MRPC.all.params.edge.based.pdf` | MRGN vs MRPC edge-based performance comparison. |
| `rare_trio_distributions/` | Expression distribution plots for individual rare-MAF trios. |

### `supplementary_tables/`

Supplementary CSV/Excel tables:

| File | Description |
|---|---|
| `ST_all_results_simulation.csv` | Full simulation results across all methods and metrics. |
| `ST_compute_times.csv` | Computation time comparison table. |
| `ST_GTEx_all_trios_master.csv` | Master table of all GTEx trio results. |
| `ST_sign_rank_tests_skewness_in_trans_genes.csv` | Signed rank test results for skewness in trans genes. |
| `ST10_Model_miss_Type_II_error_rates.csv` | Type II error rates under model misspecification. |
| `S7_GTEx_comp_MRGN_MRPC_GMAC.csv` | GTEx comparison of MRGN, MRPC, and GMAC. |
| `S8_GTEx_comp_MRGN_GMAC_perm_all.csv` | GTEx comparison with permutation test results. |
| `T1.T2.edge.results*.csv` | T1–T2 edge-based metrics for standard and many-confounder simulations. |
| `Results-class-based.csv` / `Results-edge-based.csv` | Class-based and edge-based metric summaries. |

### `supplementary_text/`

| File | Description |
|---|---|
| `Supplementary-Textupdated.docx` | Supplementary text document for the manuscript. |

### `other/`

Helper functions, intermediate data, and auxiliary outputs.

| File | Description |
|---|---|
| `helpers.R` | Reusable functions for generating class-based and edge-based performance metric tables. Defines `generate_class_based_metrics()` and `generate_edge_based_metrics()`. |
| `results_all_gmac.RData` | Aggregated GMAC results. |
| `results_all_MRPC-ADDIS.RData` | Aggregated MRPC-ADDIS results. |
| `results_all_MRPC-LOND.RData` | Aggregated MRPC-LOND results. |
| `results_all_times.RData` | Computation time data for all methods. |
| `results_all_times_mrgn.RData` | MRGN-specific computation time data. |
| `TUB-class-based.csv` | True upper bound — class-based metrics. |
| `TUB-edge-based.csv` | True upper bound — edge-based metrics. |
| `TUB-t1-t2-edge.csv` | True upper bound — T1–T2 edge-based metrics. |
| `TUB-mrgn-all-trios.RData` | True upper bound MRGN results for all trios. |
| `Bad_simulated_trio.csv` | Example of a problematic simulated trio. |
| `tablescraps/` | Working/draft table fragments. |

### `versions/`

LaTeX manuscript source and compiled output. Contains the manuscript `.tex` file, compiled PDF, BibTeX references, and LaTeX style packages (authblk, balance, figcaps, fullpage, sublabel).
