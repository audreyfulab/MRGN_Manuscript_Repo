# `simulation/simulated_data/`

Output of `../updated_data_simulation.R`. **Not tracked in git** — the files run to hundreds of MB,
well past GitHub's limits. Regenerate by rerunning the driver from the repository root;
it seeds with `set.seed(234)`, so a rerun on the same R version reproduces the same
trios.

| file | size | what it is |
| --- | --- | --- |
| `simulated_trios.RData` | 146 MB | **current** — the 1,500 trios the inference stage consumes. Regenerated 2026-08-22 under raised `b.snp`/`b.med` floors and narrowed `U`/`W`/`Z` coefficient ranges; the run it superseded is described in [`../../simulation_results/legacy/first_pass/README.md`](../../simulation_results/legacy/first_pass/README.md) |
| `simulated_trios_precalibration.RData` | 367 MB | superseded — 3,750 trios at a single sample size, from before the effect-size and confounder recalibration; kept for comparison |
| `simulated_trios_u_only.RData` | 43 MB | 300 trios, n = 670, confounders only (`W_n = 0, Z_n = 0`). From `../confounder_structure_simulation.R`, seed 2341 |
| `simulated_trios_u_w.RData` | 42 MB | 300 trios, n = 670, confounders + 1 intermediate (`W_n = 1, Z_n = 0`). Seed 2342 |
| `simulated_trios_u_z.RData` | 43 MB | 300 trios, n = 670, confounders + 1 common child (`W_n = 0, Z_n = 1`). Seed 2343 |

The three `_u_*` files are the confounder-structure cases, MRGN only; the fourth structure
(confounders + intermediate + common child) is the n = 670 group of `simulated_trios.RData`
and is not duplicated. They share `simulated_trios.RData`'s structure exactly — the only
difference is which of the `W`/`Z` columns exist.

## Structure

Both files load a single object named `data`: a list with one element per dataset — 1,500
in `simulated_trios.RData`, 3,750 in the precalibration file — each element
a three-element list.

```r
sim <- MRGN::loadRData("bioinfo_revision/simulation/simulated_data/simulated_trios.RData")
length(sim)            # 1500
names(sim[[1]])        # "data" "params" "conf.effects"
dim(sim[[1]]$data)     # 50 24
```

### `$data` — the dataset itself

A `sample.size x (3 + K_n + U_n + 2)` data frame:

| columns | contents |
| --- | --- |
| 1 | `V1` — genotype |
| 2 | `T1` — cis gene |
| 3 | `T2` — trans gene |
| next `K_n` | clinical covariates `pcr`, `platform`, `sex` (n = 670 datasets only) |
| next `U_n` | `U1.<i>` … `U<U_n>.<i>` — unobserved confounders |
| then | `W1.<i>` — intermediate |
| last | `Z1.<i>` — common child |

Every simulated confounder is suffixed with `.<dataset index>` by
`name.trio.columns()`. That suffix is what lets the inference stage pool all the
confounders of a sample-size group into one matrix and still attribute a selected
column back to its trio and its block — see `own.block.names()` in
`../../simulation_results/inference_utils.R`.

### `$params` — a one-row data frame

`dataset`, `model`, `sample.size`, `effect_size`, `replicate`, `minor.freq`, `b.snp`,
`b.med`, `SD`, `W_n`, `Z_n`, `K_n`, `U_n`, `n.resamples`, `R2.T1.U`, `R2.T2.U`.

`n.resamples` is how many draws it took before all three genotype classes appeared —
up to 1,577 at θ = 0.01, n = 50, which inflates the realized MAF roughly 1.7× for
nominal MAF ≤ 0.05 (see [`../../METHODS.md`](../../METHODS.md) §2).

`R2.T1.U` / `R2.T2.U` are **unadjusted** `R²`, inflated whenever `U_n` is an
appreciable fraction of `n`. Use `../verify_simulation.R`'s adjusted figures for any
comparison against the real 0.41 / 0.31.

### `$conf.effects`

`list(cis = ..., trans = ...)`, each the `summary(lm(...))$coefficients` matrix from
regressing that gene on *all* covariates. These condition on the collider `Z` and the
mediator `W` and omit `V1`, so they are **not** comparable to the pools in
`pc_distribution_invest/data/`, which regress on PCs only.

## Grid coverage

`simulated_trios.RData` is balanced by construction:

| factor | levels | datasets each |
| --- | --- | --- |
| `model` | model0 … model4 | 750 |
| `sample.size` | 50, 150, 300, 670, 1000 | 750 |
| `effect_size` | small, medium, large | 1,250 |

## Why the precalibration file is kept

It was generated with the shared effect-size range and the pre-revision confounder
bound, and is the "before" side of the recalibration described in `METHODS.md` §3:

| | `b.snp` range | `b.med` range | median unadjusted `R2.T1.U` |
| --- | --- | --- | --- |
| `simulated_trios_precalibration.RData` | 0.10 – 1.00 | 0.10 – 1.00 | 0.725 |
| `simulated_trios.RData` | 0.05 – 1.50 | 0.05 – 1.00 | 0.456 |

Nothing in the pipeline reads it. It can be passed to `verify_simulation.R` as its
optional argument to reproduce the before-and-after comparison, and can be deleted once
the revision is settled.
