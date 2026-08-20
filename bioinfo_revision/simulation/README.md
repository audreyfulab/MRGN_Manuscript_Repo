# `simulation/` — data generation for the revised study

Generates the 3,750 simulated trios the revision's inference runs on, and checks that
the confounding they carry matches real GTEx Whole Blood.

The methodology and every measured calibration number live in
[`../METHODS.md`](../METHODS.md); this file describes the files and how to run them.

## Contents

| file | what it is |
| --- | --- |
| `updated_data_simulation.R` | **driver** — builds the scenario grid and generates every dataset |
| `simulation_utils.R` | helper functions sourced by the driver; no top-level side effects, safe to source anywhere |
| `verify_simulation.R` | calibration check — simulated vs real confounding, plus a diagnosis of why the confounder filter struggles at small `n` |
| [`simulated_data/`](simulated_data/) | the generated trios (untracked, ~368 MB per file) |

## Running

From the **repository root**:

```r
setwd("path/to/MRGN_Manuscript_Repo")
source("bioinfo_revision/simulation/updated_data_simulation.R")
```

then the calibration check, either in R or from a shell:

```sh
Rscript bioinfo_revision/simulation/verify_simulation.R
Rscript bioinfo_revision/simulation/verify_simulation.R path/to/some_other_trios.RData
```

`verify_simulation.R` takes one optional argument: the simulated-data file to check.
It defaults to `bioinfo_revision/simulation/simulated_data/simulated_trios.RData`.
(The usage lines in the script's own header comment still show the pre-reorganization
path `bioinfo_revision/verify_simulation.R`; the file now lives in this folder.)

Requires `MRGN`, `mvtnorm`, `gridExtra` for the driver and `MRGN`, `ggplot2` for the
check. The driver seeds with `set.seed(234)` and reads one GTEx input,
`GTEx/data/kclist_top5_tiss.RData` — the clinical covariates (`pcr`, `platform`, `sex`)
for the 670 Whole Blood donors, used as the `K` block at n = 670 only.

`verify_simulation.R` additionally reads
`pc_distribution_invest/data/real_pc_effect_pools.RData`, so run that stage first.

## `updated_data_simulation.R`

Builds `scenarios`, one row per dataset:

```
5 models x 5 sample sizes (50, 150, 300, 670, 1000) x 3 effect sizes x 50 replicates
  = 3,750 datasets
```

then draws the remaining parameters per row and hands the table to
`simulate.all.datasets()`. Parameters drawn in the driver:

| column | draw |
| --- | --- |
| `minor.freq` | `Uniform{0.01, 0.02, ..., 0.50}` |
| `b.snp`, `b.med` | by `draw.effect.sizes()`, see below |
| `SD` | fixed at 1 (Yang et al. 2017) |
| `W_n`, `Z_n` | 1 each |
| `K_n` | 3 at n = 670, else 0 |
| `U_n` | `Uniform{1..50}`, capped at `sample.size - W_n - Z_n - 4` |

`effect_sizes` is a two-level list, one tertile split per parameter:

| stratum | `b.snp` | `b.med` |
| --- | --- | --- |
| small | 0.05–0.50 | 0.05–0.35 |
| medium | 0.55–1.00 | 0.40–0.70 |
| large | 1.05–1.50 | 0.75–1.00 |

The SNP gets the wider range because it drives everything downstream — `cor(V1, T1)`
caps `cor(V1, W)` and `cor(V1, Z)`. Neither range may start at 0: `b.snp = 0` removes
the `V1 → T1` edge so a "model0" trio would carry no edges and its `M0.1` truth label
would be wrong, and `b.med = 0` does the same to model1, model2 and model4.
`draw.effect.sizes()` errors on a zero lower bound.

Two traps the code guards against explicitly, worth knowing before editing:

- `scenarios` is built with `stringsAsFactors = FALSE`. `gen.graph.skel()` dispatches
  on the model with `switch(model, model0 = ..., ...)`, and a factor is silently
  treated as its integer code.
- `effect_sizes` is a `list()`, not `c()`. `c(small = c(0.1, 0.3), ...)` flattens and
  renames the elements `small1`, `small2`, so `effect_sizes[["small"]]` would be a
  subscript error.

The driver `setwd()`s into `simulated_data/` to write the output and restores the root
afterwards.

## `simulation_utils.R`

| function | what it does |
| --- | --- |
| `draw.effect.sizes(scenarios, effect_sizes, step = 0.05)` | fills the `b.snp` and `b.med` columns per effect-size stratum, drawing the two independently within a stratum; errors if a range starts at 0 |
| `simulate.dataset(settings, clinical.covs, verbose)` | one dataset: draws the `U` block from `rmvnorm(mean = 0, sigma = I)`, prepends the `K` block when `K_n > 0`, and calls `MRGN::simData.from.graph()` in a loop until all three genotype classes appear; returns the data and the resample count |
| `name.trio.columns(data, index, K_n, kc.names)` | gives the `K` columns their clinical names and suffixes every simulated confounder with `.<dataset index>`, which is what keeps names unique once the confounders of all datasets are pooled |
| `conf.r.squared(data, K_n, U_n)` | **unadjusted** `R²` of T1 and T2 on their own `U` block; `NA` for a saturated fit |
| `write.rdata(data, filename)` | `save()` wrapper, defaults to `simulated_trios.RData` in the working directory |
| `simulate.all.datasets(...)` | loops over `scenarios`, assembles the three-element record per dataset, optionally saves |

### Generating models and confounder blocks

Topology comes from `MRGN::gen.graph.skel()`, data from `MRGN::simData.from.graph()`.
The cis gene is always `T1`, so the truth label is always the `.1` variant.

| model | trio structure | truth label |
| --- | --- | --- |
| model0 | `V1 → T1`, `T2` independent | M0.1 |
| model1 | `V1 → T1 → T2` | M1.1 |
| model2 | `V1 → T1 ← T2` | M2.1 |
| model3 | `V1 → T1`, `V1 → T2` | M3 |
| model4 | `V1 → T1 → T2`, `V1 → T2` | M4 |

| block | role | count | `conf.coef.ranges` |
| --- | --- | --- | --- |
| `K` | known clinical covariates | 3 at n = 670, else 0 | `c(0, 0)` — deliberately null |
| `U` | unobserved confounders, `U → T1` and `U → T2` | `Uniform{1..50}` | `c(0, 0.3)` |
| `W` | intermediate, `T1 → W → T2` | 1 | `c(0.05, 0.5)` |
| `Z` | common child, `T1 → Z ← T2` | 1 | `c(1, 1.5)` |

All four entries must stay in the `conf.coef.ranges` list even when a block is empty:
`gen.graph.skel()` indexes it positionally (K, U, W, Z). `conf.mat` must carry column
names matching the node names (`U1`..`Ud`, `K1`..`K3`) — an unnamed matrix silently
produces all-`NA` confounder columns that propagate into T1/T2 as `NaN`, and an unnamed
data frame errors with "undefined columns selected".

> **Note on the `U` bound.** The code currently uses `U = c(0, 0.3)`; `METHODS.md` §3
> describes the change as `c(0, 0.2)`, the central 95% of the real per-PC effects. The
> data in `simulated_data/simulated_trios.RData` was generated with 0.3, which lands the
> realized adjusted `R²` at ~0.37 against the real 0.41. Treat the script as
> authoritative and the METHODS figure as the interval it was derived from.

## `verify_simulation.R`

Read-only. Nothing it computes feeds back into the simulation; it reports measured
against real and writes one figure,
`../simulation_results/simulated_vs_real_conf_effects.png`.

Five sections, all printed to the console:

1. **Per-confounder standardized effect**, `cor(U_i, gene)`, simulated vs the real
   pool. Restricted to n = 670 datasets — a sample correlation carries about `1/√n` of
   noise, so only the sample size the real pool was measured at is a like-for-like
   comparison.
2. **Adjusted `R²`** of each gene on its own `U` block, by sample size, against the real
   medians of 0.412 (cis) and 0.307 (trans). Adjusted, not raw: at n = 50 with `U_n = 26`
   the raw value is inflated by about 0.35 whatever the truth is, which is the trap
   `conf.r.squared()` falls into.
3. **`R²` against confounder count**, `n ≥ 300`, binned by `U_n` against the real curve.
4. **SNP effect** — `median |cor(V1, T1)|` overall and by effect-size stratum.
5. **Intermediate / common-child detectability** — whether `get.conf.trios()` can see
   `W` and `Z`, separating the three obstacles (effect size, sample size, multiplicity)
   by comparing the pass rate at the pooled q-value threshold against the rate at an
   uncorrected α = 0.05, and de-noising the observed correlations with
   `√(mean(r²) − 1/n)`. This is the section to rerun after any change to the effect
   ranges — it is the quickest check on whether the confounder filter has anything to
   work with.

Reference values are hard-coded at the top in `REAL`, sourced from
`pc_distribution_invest/`.
