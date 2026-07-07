# Manuscript

All materials for the MRGN manuscript: figure/table generation scripts, final outputs, supplementary materials, and helper functions.

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
| `create_combined_T1T2_edge_table.R` | Standalone script that combines the ≤15 and 15–50 confounder T1–T2 edge results into a single 1800-trio table with confusion matrices, precision, recall, Type I error, and Type II error for each method. Saves to `supplementary_tables/T1.T2.Edge.Results.Combined.csv`. |
| `check_mrpc_comp_times.R` | Analyzes MRPC-ADDIS computation times. Plots time vs. number of confounders, identifies problematic trios (>4 hours), and benchmarks against MRGN/GMAC. |
| `gmac_valid_res_6_27_2023.R` | Loads GMAC validation results and creates supplementary figure SF9: GMAC Type I error rate vs. number of confounders at 0.01 and 0.05 significance cutoffs. |
| `Plot_dist_rare_genotype_trios.R` | Generates supplementary figures showing gene expression distributions for GTEx WholeBlood trios with rare MAF (<10%). Includes scatter plots and histograms with skewness/kurtosis statistics. |
| `plot_SF_perm_evidence.R` | Creates figures visualizing the permutation test's effect on MRGN inference. Classifies model changes as "liberal" or "conservative" and shows impact on rare-allele trios. Also exports GTEx result tables. |
| `MRGN_write_up_helper_functions.R` | Shared helper functions used across multiple manuscript scripts. Provides utilities for loading RData files by name, formatting results, and other common operations. |

**Data files in `scripts/`:**

| File | Description |
|---|---|
| `loadedResults.R` | (see above — also serves as a data pipeline) |
| `mrgn_v_gmac_v_mrpc_5k_datasets_all_mods_conf_types.RData` | Combined simulation results for all methods (5K dataset). |
| `mrpc_5k_conf_list_all_confs_filtered_all_mods.RData` | MRPC confounder lists for filtered simulations. |
| `mrpc_v_mrgn_v_gmac_5k_params_all_mods_conf_types.RData` | Simulation parameters for all methods. |
| `all_trios_output_cis.Rdata` | Cis-mediation output for all trios. |

### `figures/`

Final main manuscript figures (PDF).

### `supplementary_figures/`

Supplementary figures (PDF).

### `supplementary_tables/`

Supplementary CSV/Excel tables.

### `supplementary_text/`

Supplementary text for simulation study on model misspecification.

### `other/`

Helper functions, intermediate data, and auxiliary outputs.

| File | Description |
|---|---|
| `helpers.R` | Reusable functions for generating class-based and edge-based performance metric tables. Defines `generate_class_based_metrics()`, `generate_edge_based_metrics()`, and `generate_t1_t2_results()` (which produces the combined T1–T2 edge table with Type I/II error). |
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
| `mrpc_ADDIS_time_to_compute_trios.pdf` | Plot of MRPC-ADDIS computation time per trio. |
| `mrpc_LOND_time_to_compute_trios.pdf` | Plot of MRPC-LOND computation time per trio. |
| `mrpc_time_to_compute_trios.pdf` | Plot of overall MRPC computation time per trio. |
| `tablescraps/` | Working/draft table fragments. |
