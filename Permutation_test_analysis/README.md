# Permutation Test Analysis

Analysis of the effect of the MRGN permutation test on GTEx WholeBlood trio inference. Compares MRGN results with and without the permutation test under different confounder selection strategies.

## Overview

The permutation test is an optional step in MRGN's `infer.trio()` function (`use.perm = TRUE/FALSE`) that provides additional validation of mediation signals. This folder investigates how enabling the permutation test changes model assignments, particularly for trios with rare minor allele frequencies.

Two confounder selection strategies are compared:
1. **Standard (manuscript settings)** — FDR-corrected confounder selection
2. **Liberal (no FDR)** — Alpha < 0.05 cutoff without FDR correction

## Subfolders

### `scripts/`

| File | Description |
|---|---|
| `investigate.perm.R` | Runs `infer.trio()` with and without the permutation test on GTEx WholeBlood trios under both the standard and liberal confounder selection strategies. Saves paired results for comparison. |
| `investigate.perm.liberal.confsel.noFDR.R` | Same analysis as above but exclusively uses the liberal confounder selection (no FDR, alpha < 0.01). |

### `results/`

Output files from the permutation test analysis:

| File | Description |
|---|---|
| `no.perm.all.trios.WB.RData` | MRGN results **without** permutation test (standard confounder selection). |
| `perm.all.trios.WB.RData` | MRGN results **with** permutation test (standard confounder selection). |
| `no.perm.all.trios.WB.liberal.confs.alpha05.RData` | MRGN results **without** permutation test (liberal confounder selection, alpha < 0.05). |
| `perm.all.trios.WB.liberal.confs.alpha05.RData` | MRGN results **with** permutation test (liberal confounder selection, alpha < 0.05). |
| `GTEx_results_*.txt` | Human-readable summary tables of inference results for each combination of permutation/no-permutation and confounder selection strategy. |
| `GTEx_Analysis_All_Results_Table_WholeBlood*.txt` | Comprehensive results tables for all WholeBlood trios, including a rare-trios-only variant. |
