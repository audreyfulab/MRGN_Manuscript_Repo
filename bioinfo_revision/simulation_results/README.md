# `simulation_results/` — inference over the simulated trios

Runs MRGN, MRPC and GMAC on every simulated trio and consolidates the six inferences,
the generating parameters, and a scoring of each confounder selection into one wide
data frame with one row per trio.

Rationale and measured results are in [`../METHODS.md`](../METHODS.md) §4; this file
describes the files and how to run them.

## Contents

| file | what it is |
| --- | --- |
| `updated_simulation_inference.R` | **driver** — configuration, per-group orchestration, checkpointing, and the final combine |
| `inference_utils.R` | helpers: confounder selection, and one `apply.*()` per method |
| `simulated_vs_real_conf_effects.png` | figure — written by `../simulation/verify_simulation.R`, not by anything here |

Running the driver additionally writes, into this folder:

| file | when |
| --- | --- |
| `inference_group_n50.RData` … `inference_group_n1000.RData` | one checkpoint per sample-size group, saved as that group finishes |
| `inference_results.RData` | all groups combined, object name `inference.results` |
| `inference_results.csv` | the same table as CSV |

None of these are in the repository yet — they appear on the first full run, and are
large enough that they should stay untracked.

## Running

From the **repository root**:

```r
setwd("path/to/MRGN_Manuscript_Repo")
source("bioinfo_revision/simulation_results/updated_simulation_inference.R")
```

Sourcing the file runs the whole thing: the last two lines call `run.all.inference()`
and then `combine.groups()`. To drive it manually, comment those two lines out and call
the functions yourself — `run.all.inference(sample.sizes = 670)` for one group,
`combine.groups(sample.sizes = c(50, 150))` to stitch a subset of the checkpoints.

**Smoke test first.** Set `max.per.group` (top of the file) to something small, e.g.
`max.per.group <- 5`, and optionally `n.bootstrap <- 50`, to run the pipeline end to
end in minutes. A full run is 3,750 trios × (3 MRGN fits, each with 1,000 bootstrap
replicates, + 2 MRPC fits + a per-trio GMAC timing call) plus two `get.conf.trios()`
selections and two group-wide `gmac()` calls per sample size — plan on a long run, and
note that a group whose checkpoint already exists is *not* skipped on a rerun.

### Configuration block

| setting | default | meaning |
| --- | --- | --- |
| `sim.data.file` | `bioinfo_revision/simulation/simulated_data/simulated_trios.RData` | input |
| `clinical.file` | `./GTEx/data/kclist_top5_tiss.RData` | clinical covariates |
| `out.dir` | `bioinfo_revision/simulation_results` | where checkpoints and results go |
| `n.bootstrap` | 1000 | bootstrap replicates per MRGN fit |
| `mrpc.timeout` | 120 s | MRPC has taken hours on trios with many confounders |
| `gmac.nperm` | 1000 | GMAC permutations |
| `selection.alpha` | 0.05 | cutoff for the GMAC mediation calls |
| `n.cores` | `detectCores() - 2` | one cluster shared by GMAC and the MRGN bootstrap |
| `max.per.group` | `NULL` | cap the trios per group, for smoke tests |

Requires `MRGN`, `MRPC`, `R.utils`, `parallel`, and sources
`../simulation/simulation_utils.R` plus the adapted GMAC code in
`../../adapted_GMAC_func/`. Cluster workers `setwd()` to the recorded root and source
the GMAC files themselves — `gmac()` dispatches `child.p`/`conf.fdr`/`getp.func`
through `parLapply`, and those helpers live in sourced files rather than a package. The
stream is seeded with `clusterSetRNGStream(cl, 234)`.

## Why it processes one sample size at a time

Both `get.conf.trios()` and `gmac()` index a covariate pool **by row**, so every trio in
a call must have the same number of observations. The driver therefore groups the
datasets by sample size, builds one pooled covariate matrix per group from every trio's
own `U`/`W`/`Z` columns (~21,000 columns over 750 trios), and checkpoints each group as
it finishes. The `K` block is passed separately as `known.conf` rather than pooled,
since all trios in a group share the same three clinical columns.

Two consequences of the pooled design, recorded rather than fixed (`METHODS.md` §4):
the q-value filter tests ~15.8M correlations at once, so a covariate must clear
`|r| ≈ 0.55` at n = 50 to survive; and CS-α at α = 0.01 over 21,000 columns implies
~210 false-positive confounders per trio. Every method call is therefore wrapped in the
local `safely()`, which converts an error into a row of `NA`s plus an error message
instead of taking down the run.

## The six inferences per trio

| column prefix | method | confounders |
| --- | --- | --- |
| `mrgn.truth.` | MRGN | the true ones: trio + `K` + that trio's own `U` block |
| `mrgn.CSq.` | MRGN | CS-q selected |
| `mrgn.CSa.` | MRGN | CS-α selected |
| `mrpc.CSq.` | MRPC | CS-q selected |
| `mrpc.CSa.` | MRPC | CS-α selected |
| `gmac.` | GMAC | whatever GMAC selects for itself |

`ground.truth.input()` deliberately excludes `W` and `Z` — the `U` variables are the
true confounders; `W` is an intermediate and `Z` a common child, and conditioning on
either is a mistake, not a baseline.

Truth labels come from an explicit `TRUTH.LABEL` lookup rather than
`MRGN::convert.truth()`, which maps by sorted position and silently mislabels when the
input does not contain all five models. The cis gene is always `T1`, so the truth is
always the `.1` variant.

## Output columns

One row per trio. Beyond `dataset` and every column of `$params` (see
[`../simulation/simulated_data/`](../simulation/simulated_data/)):

| group | columns |
| --- | --- |
| truth | `truth.model`, `n.true.confs`, `true.confs` (`;`-separated) |
| selection scoring, one block each for `CSq`, `CSa`, `gmac` | `<p>.n.selected`, `<p>.n.tp`, `<p>.n.fp`, `<p>.n.fn`, `<p>.n.fp.other.trio`, `<p>.has.common.child`, `<p>.has.intermediate`, `<p>.selected` |
| filter fallback | `CSq.filter_int_child`, `CSa.filter_int_child` — `FALSE` marks a group where `get.conf.trios()` found no intermediate or common child and the selection fell back to `filter_int_child = FALSE` |
| MRGN, per setting | `<p>.model`, `<p>.correct`, `<p>.correct.coarse`, `<p>.time.seconds`, `<p>.boot.model`, `<p>.boot.min.edge.prob`, `<p>.boot.p.V1T1`, `<p>.boot.p.T1T2`, `<p>.boot.p.V1T2`, `<p>.boot.p.T2T1`, `<p>.bootstrap.time.seconds`, `<p>.error` |
| MRPC, per setting | `<p>.model`, `<p>.correct`, `<p>.correct.coarse`, `<p>.time.seconds`, `<p>.timed.out`, `<p>.error` |
| GMAC | `gmac.model`, `gmac.cispval`, `gmac.transpval`, `gmac.ciseffect`, `gmac.transeffect`, `gmac.time.seconds`, `gmac.error` |

`n.fp.other.trio` counts false positives borrowed from a *different* trio's confounder
block, as opposed to this trio's own `W`/`Z` — the suffix in the column name is what
makes that distinguishable.

`correct.coarse` compares only the model family (`M1.1` and `M1.2` both count as `M1`).
GMAC gets **no** correctness flag: it reports a mediation call
(`"Cis Mediated"` / `"Trans Mediated"` / `"Undirected"` / `"No Mediation"`), not a model
label, and cannot separate M0 from M2/M3, so cross-tabbing it against the truth is a
decision for the analysis stage.

Each group's results carry `selection.time.seconds`, `gmac.time.seconds` and
`group.time.seconds` as attributes; note that `attr()`s do not survive the `rbind()` in
`combine.groups()`, so read them from the per-group checkpoints.

## `inference_utils.R`

Every `apply.*()` runs one method on a **single** trio and attaches the time the
inference itself took, in seconds, so the three methods are directly comparable.

| function | what it does |
| --- | --- |
| `select.confounders(trios, cov.pool, known.conf, ...)` | calls `MRGN::get.conf.trios()` **once** per group with `adjust_by = "all"`, giving **CS-q** (FDR 5%); **CS-α** (per-test α = 0.01) is derived by thresholding the returned `reg.pvalues`, which is exactly what `adjust_by = "none"` does internally. Returns, per setting, one `trio + known + selected` data frame per trio, the raw selection, the timing, and whether the intermediate/common-child filter actually applied |
| `apply.mrgn(data, bootstrap, number_of_samples, cl, ...)` | `MRGN::infer.trio()` plus optional bootstrap; `time.seconds` covers the fit only, the bootstrap is timed separately |
| `boostrap_edge_probabilities(trio, number_of_samples, cl)` | resamples rows, keeps the first six elements of each `infer.trio()` result, and returns per-indicator means, the majority-vote model, and the weakest supported edge probability |
| `apply.mrpc(data, timeout, ...)` | `MRPC()` capped with `R.utils::withTimeout()`, then classified against the eight trio adjacencies |
| `classify.mrpc.adj(adj)` / `mrpc.truth.adjacencies()` | matches an inferred 3×3 adjacency against `M0.1`, `M0.2`, `M1.1`, `M1.2`, `M2.1`, `M2.2`, `M3`, `M4`, else `"Other"`; the reference matrices are built once and cached |
| `apply.gmac(trio, confounders, known.conf, ...)` | per-trio `gmacOneTrio()` in both directions — this exists to **time** GMAC on one trio; the reported results come from `run.gmac.all()` |
| `run.gmac.all(trios, cov.pool, known.conf, ...)` | GMAC over a whole group, which is where GMAC does its own confounder selection (stratified FDR across the trios), so it cannot be done one trio at a time |
| `gmac.model.call(p.cis, p.trans, alpha)` | maps the two mediation p-values to GMAC's four-way call |

Three behaviours worth knowing before editing:

- The bootstrap keeps `b11, b12, b21, b22, V1:T1, V1:T2` rather than just the
  adjacency. The two interaction terms are what `MRGN::class.vec()` uses to separate M2
  from M4, and they are needed exactly when all four edge indicators are 1 — an averaged
  adjacency alone cannot name the supported model.
- `apply.mrgn()` reads `result$Inferred.Model` **by name**. `infer.trio()` returns it as
  element 18; the `result[[14]]` used by the original scripts now picks up `coef11`.
- `apply.mrpc()` handles its own timeout and returns a complete record
  (`model = NA`, `timed.out = TRUE`) rather than letting a driver loop lose the trio.
  `withTimeout()` interrupts between R evaluation steps, which is where MRPC spends its
  time; it cannot break into long-running compiled code.

`select.confounders()` also guards three failure modes with explicit errors, a
normalization, or a documented fallback:

- `blocksize` must not exceed the number of trios — `propagate::bigcor()` mis-blocks the
  correlation otherwise and fails later with an opaque `'x' is empty`.
- `get.conf.trios()` stops outright with "No common child or intermediate variables
  detected" when nothing in the pool associates with any variant. Caught, retried with
  `filter_int_child = FALSE`, and recorded in the `*.filter_int_child` output columns.
  Because one call now serves both settings, the two columns always agree.
- `sig.asso.covs` comes back from `apply(reg.sigmat, 1, which)`, which simplifies to a
  **matrix** if every trio happens to select the same number of covariates. Everything
  downstream indexes it with `[[i]]`, so `as.cov.list()` normalizes it to a list of
  integer vectors first.

## `simulated_vs_real_conf_effects.png`

Written by `../simulation/verify_simulation.R`, and lives here because it is a result
rather than an input. Two facets, cis and trans gene, each overlaying the density of
the simulated per-confounder standardized effects (n = 670 datasets) on the real GTEx
Whole Blood pool, with each source's `sd` and `n` printed into the panel.

The two distributions overlay closely (simulated sd 0.123 against real 0.117 cis /
0.104 trans). Two visible differences, both expected: the real densities dip at zero,
because those PCs were *selected* at FDR 0.05 so near-null effects are filtered out by
construction, and the real tails reach ±0.66 where the simulated ones stop near ±0.35.

![simulated vs real confounder effects](simulated_vs_real_conf_effects.png)
