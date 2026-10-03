# `bioinfo_revision/gtex_rerun/` — GTEx inference without the known covariates

The GTEx v8 whole-blood expression was PEER-normalized with the GTEx-reported covariates
already regressed out. The original analysis then added three of them — `pcr`, `platform`
and `sex` (columns 4–6 of every trio data frame, from `GTEx/data/kclist_top5_tiss.RData`) —
again as known confounders, so they were adjusted for twice. This folder reruns the
**inference step only** with those three covariates removed, and rebuilds manuscript
Figure 4 from the result.

Confounder selection is **not** rerun: MRGN uses the same CS-q and CS-α PC sets as before,
and GMAC's internal PC selection never reads the known covariates, so it selects the same
PCs as before (the script asserts this).

## What was rerun

| method | input (covariate columns 4–6 dropped) | settings | rerun? |
| --- | --- | --- | --- |
| MRGN, CS-q | `GTEx/data/data.with.PCs.WholeBlood.RData` | `infer.trio(use.perm = FALSE)` | yes |
| MRGN, CS-α | `GTEx/trios_data_prlnc_liberal_conf_sel/data_with_confs/trios.with.confs.WholeBlood.RData` | `infer.trio(use.perm = FALSE)` | yes |
| GMAC | `GTEx/results/GMAC/GMAC_input_list_for_WholeBlood.RData`, empty known-confounder matrix | as `GTEx/results/GMAC/gtex.analysis.gmac.R`: nperm = 1000, nominal p, fdr = 0.05, fdr_filter = 0.1; cis and trans mediator runs | yes |
| MRPC-ADDIS | — | — | **no** |

The MRGN settings are the ones that produced the existing `MRGN.Inferred.Model.CSq` /
`.Csalpha` columns. Before each arm runs, the script reruns 100 random trios **with** the
covariates kept and stops unless all 100 reproduce the stored calls.

MRPC was not rerun. On a test of five trios at n = 670 with 18–38 CS-q confounders, three
did not finish within 10 minutes, so a full run was not practical. The original
`MRPC.Addis.Inferred.Model` column is unchanged and is what the figure plots for MRPC.

## Scripts

Run both from the **repository root**, in this order.

| script | does | writes |
| --- | --- | --- |
| `rerun_gtex_no_known_covs.R` | MRGN (both arms) and GMAC (cis + trans) without the covariates; aligns GMAC to the master table by trio IDs; appends the new columns to the master table; old-vs-new summaries | `updated_results/*`, new columns in `TableS3_GTEx_all_trios_master_updatedOct1.csv` |
| `make_gtex_fig_no_known_covs.R` | Rebuilds the GTEx figure (style from `Manuscript/scripts/create_GTEx_figs.R`, panel order of the v10 manuscript figure) for one MRGN confounder set, chosen with `--arm=CSq` (default) or `--arm=CSalpha`, and prints the counts quoted in the text | `updated_figures/fig_gtex_noKC_MRGN_<arm>.{pdf,png}` |

```
Rscript bioinfo_revision/gtex_rerun/rerun_gtex_no_known_covs.R --cores=4
Rscript bioinfo_revision/gtex_rerun/make_gtex_fig_no_known_covs.R --arm=CSq
Rscript bioinfo_revision/gtex_rerun/make_gtex_fig_no_known_covs.R --arm=CSalpha
```

`--cores` sets the size of the parallel cluster (default 4). The trio lists take several GB
in memory, so the script loads, runs and releases one arm at a time. Each stage
checkpoints to `updated_results/` and is skipped on a rerun if its checkpoint exists. To
force a stage to rerun, delete its `.RData`.

Requirements: R 4.4.x with `MRGN` (≥ 1.3.9), `qvalue`, `parallel`, and for the figure
`ggpubr`, `ggthemes`, `gridExtra`. GMAC comes from `adapted_GMAC_func/GMAC_moded/R/GMAC.R`.
The two MRGN trio lists are gitignored GTEx inputs and must be present locally.

## Files

### `TableS3_GTEx_all_trios_master_updatedOct1.csv`

The supplementary table of all 3,248 whole-blood trios. The original columns are kept
unchanged, and the rerun appends these at the end:

| column | meaning |
| --- | --- |
| `MRGN.Inferred.Model.CSq.noKC` | MRGN model, CS-q confounders, no known covariates |
| `MRGN.Inferred.Model.Csalpha.noKC` | MRGN model, CS-α confounders, no known covariates |
| `GMAC.Inferred.Mediation.Type.Qval.noKC` | GMAC mediation call, q-value < 0.1 |
| `GMAC.Inferred.Mediation.Type.alpha01.noKC` | GMAC mediation call, p < 0.01 (used in Figure 4C) |
| `GMAC.Inferred.Mediation.Type.alpha05.noKC` | GMAC mediation call, p < 0.05 |
| `GMAC.cis.pval.noKC`, `GMAC.trans.pval.noKC` | GMAC mediation p-values (known + selected PCs), cis and trans gene as mediator |

GMAC calls are one of `Cis Mediated`, `Trans Mediated`, `Undirected` (both) or `No
Mediation`, built as in `GTEx/scripts/make_mrgn_triotables.R`.

### `updated_results/`

| file | contents |
| --- | --- |
| `mrgn_noKC_CSq.RData`, `mrgn_noKC_CSalpha.RData` | `res`: 18 × 3248 `infer.trio` output, one column per trio in master-table order. Rows 1–6 are the indicators b11, b12, b21, b22, V1:T1, V1:T2; 7–12 their p-values; 13 minor allele frequency; 14–17 coefficients (new in MRGN 1.3.9); 18 `Inferred.Model` |
| `gmac_noKC_cis.RData`, `gmac_noKC_trans.RData` | `out`: raw `gmac()` output (`pvals`, `beta.change`, `sel.conf.ind`, `comp.time`). Rows are in GMAC input order, which matches the master table |
| `gmac_noKC_table.RData` | `gmac.tab`: post-processed GMAC p-values, effect changes, mediation calls and PC counts, aligned to the master table |
| `TableS3_GTEx_all_trios_master_updatedOct1.csv` | copy of the updated master table |
| `summary_old_vs_new.csv` | count of each call per method, with vs. without the known covariates |
| `agreement_old_vs_new.csv` | fraction of trios whose call is unchanged, per method |
| `rerun.log` | run log. The bare `[1] "..."` lines are GMAC's own progress prints |

### `updated_figures/`

| file | contents |
| --- | --- |
| `fig_gtex_noKC_MRGN_CSq.pdf` / `.png` | Manuscript Figure 4 rebuilt, MRGN on CS-q confounders |
| `fig_gtex_noKC_MRGN_CSalpha.pdf` / `.png` | The same figure with MRGN on CS-α confounders (cf. Supplementary Figure 9) |

Both figures have the same panels: (A) selected PCs per trio for that MRGN confounder set; (B)
inferred models, MRGN vs MRPC; (C) T1–T2 edge, all three methods; (D) cis- vs trans-gene
mediation, MRGN vs MRPC. MRGN and GMAC come from the rerun without covariates, and MRPC
from the original run, so the GMAC and MRPC bars are identical in the two figures. Panel
C's MRGN indicators come from the same run as its model labels.
