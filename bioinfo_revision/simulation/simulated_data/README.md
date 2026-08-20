# `simulation/simulated_data/`

Output of `../updated_data_simulation.R`. **Not tracked in git** — each file is ~368 MB,
well past GitHub's limits. Regenerate by rerunning the driver from the repository root;
it seeds with `set.seed(234)`, so a rerun on the same R version reproduces the same
trios.

| file | size | what it is |
| --- | --- | --- |
| `simulated_trios.RData` | 368 MB | **current** — the 3,750 trios the inference stage consumes |
| `simulated_trios_precalibration.RData` | 367 MB | superseded — the run before the effect-size and confounder recalibration, kept for comparison |

## Structure

Both files load a single object named `data`: a length-3,750 list, one element per
dataset, each a three-element list.

```r
sim <- MRGN::loadRData("bioinfo_revision/simulation/simulated_data/simulated_trios.RData")
length(sim)            # 3750
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
`../../simulation_results/updated_simulation_inference.R`.

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
