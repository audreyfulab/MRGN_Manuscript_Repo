# MRGN Manuscript Repository

This repository contains all code, data, and manuscript materials for the **MRGN (Mendelian Randomization for Genomic Network)** method paper. It includes simulation studies comparing MRGN to GMAC and MRPC, real-data application to GTEx gene-expression trios, permutation test analyses, and all scripts for generating manuscript figures and tables.

## Repository Structure

| Folder | Description |
|---|---|
| [`adapted_GMAC_func/`](adapted_GMAC_func/) | Modified GMAC source code with timing instrumentation for confounder selection and inference steps |
| [`Simulation/`](Simulation/) | Simulation study: data generation, inference scripts (MRGN, MRPC, GMAC), and results for all confounder types |
| [`GTEx/`](GTEx/) | GTEx real-data analysis: trio data, PC selection, confounder selection, and inference results for WholeBlood tissue |
| [`Manuscript/`](Manuscript/) | Figure/table generation scripts, supplementary materials, and helper functions |

## Methods Compared

- **MRGN** — Mendelian Randomization for Genomic Networks
- **GMAC** — Genomic Mediation Analysis with Adaptive Confounding adjustment
- **MRPC** — Mendelian Randomization using the PC algorithm (with ADDIS and LOND for online FDR control)

## Key Analyses

1. **Simulation Study** — 1,500 trios (5 causal models × 300 each) with 1–15 confounders, plus 300 trios with 15–50 confounders. Performance evaluated via class-based, edge-based, and T1–T2 edge-based metrics.
2. **GTEx Application** — Protein-coding and lncRNA trios from GTEx V8 WholeBlood tissue, with PC-based and regression-based confounder selection.
3. **GMAC Validation** — Type I error assessment of GMAC across varying confounder counts.
4. **Model Misspecification** — Analysis of MRGN robustness under STMS, LTMS, and SIMS misspecification scenarios.

## Requirements

The analyses are implemented in **R**. Key packages used include:
- `MRGN` — the proposed method
- `GMAC` — adapted version included in `adapted_GMAC_func/GMAC_moded/`
- `MRPC` — Mendelian Randomization PC algorithm
- `psych`, `qvalue` — for correlation testing and FDR correction
- `ggplot2`, `gridExtra` — for figure generation

See individual subfolder READMEs for detailed descriptions of each file.
