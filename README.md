# MRGN Manuscript Repository

This repository contains all code, data, and manuscript materials for the **MRGN (Mendelian Randomization with Graph-based Nuisance-variable selection)** method paper. It includes simulation studies comparing MRGN to GMAC and MRPC, real-data application to GTEx gene-expression trios, permutation test analyses, and all scripts for generating manuscript figures and tables.

## Repository Structure

| Folder | Description |
|---|---|
| [`adapted_GMAC_func/`](adapted_GMAC_func/) | Modified GMAC source code with timing instrumentation for confounder selection and inference steps |
| [`Simulation/`](Simulation/) | Simulation study: data generation, inference scripts (MRGN, MRPC, GMAC), and results for all confounder types |
| [`GTEx/`](GTEx/) | GTEx real-data analysis: trio data, PC selection, confounder selection, and inference results for WholeBlood tissue |
| [`Permutation_test_analysis/`](Permutation_test_analysis/) | Analysis of the MRGN permutation test on GTEx WholeBlood trios under different confounder selection strategies |
| [`Manuscript/`](Manuscript/) | Manuscript LaTeX source, figure/table generation scripts, supplementary materials, and helper functions |

## Methods Compared

- **MRGN** — Mendelian Randomization with Graph-based Nuisance-variable selection (the proposed method)
- **GMAC** — Genomic Mediation Analysis with Adaptive Confounding adjustment
- **MRPC** — Mendelian Randomization using the PC algorithm (ADDIS and LOND variants)

## Key Analyses

1. **Simulation Study** — 1,500 trios (5 causal models × 300 each) with 1–15 confounders, plus 300 trios with 15–50 confounders. Performance evaluated via class-based, edge-based, and T1–T2 edge-based metrics.
2. **GTEx Application** — Protein-coding and lncRNA trios from GTEx V8 WholeBlood tissue, with PC-based and regression-based confounder selection.
3. **Permutation Test Analysis** — Comparison of MRGN inference with and without the permutation test under standard and liberal confounder selection.
4. **GMAC Validation** — Type I error assessment of GMAC across varying confounder counts.
5. **Model Misspecification** — Analysis of MRGN robustness under STMS, LTMS, and SIMS misspecification scenarios.

## Requirements

The analyses are implemented in **R**. Key packages used include:
- `MRGN` — the proposed method
- `GMAC` — adapted version included in `adapted_GMAC_func/GMAC_moded/`
- `MRPC` — Mendelian Randomization PC algorithm
- `psych`, `qvalue` — for correlation testing and FDR correction
- `ggplot2`, `gridExtra` — for figure generation

See individual subfolder READMEs for detailed descriptions of each file.
