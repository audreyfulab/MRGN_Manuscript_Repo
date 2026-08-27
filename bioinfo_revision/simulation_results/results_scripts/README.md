# `results_scripts/` — summarising the simulation inference

Turns the wide per-trio tables written by the inference stage
(`../data/inference_mrgn.RData`, `../data/inference_gmac.RData`) into confusion matrices, organised
by sample size. Output goes to [`../tables/`](../tables); nothing here touches the
inference results themselves.

The inference stage records, per trio, *whether* a call was correct
(`mrgn.*.correct.coarse`) but not *what it was confused with*. These tables answer the
second question — when MRGN misses `M1`, does it land on `M0`, on `M2`, or bail to
`Other`? — and score GMAC on the one question it can answer, whether a T1–T2 edge is
present.

## Contents

| file | what it is |
| --- | --- |
| `confusion_utils.R` | label sets, the results loader, the `confusion()` builder, `scored.table()`, the Markdown writer |
| `confusion_mrgn.R` | MRGN matrices for all three confounder arms |
| `confusion_mrpc.R` | MRPC matrices, `truth`/`CSq` arms, with the timeout column scored. CS-α is excluded outright — see `inference_config.R` |
| `confusion_gmac.R` | GMAC T1–T2 edge matrices, selected and truth arms |
| `confusion_mrggi.R` | MR-GGI T1–T2 edge matrices, four arms |
| `make_all_tables.R` | **driver** — sources all four, then writes the three output files |
| `confusion_structures.R` | **standalone** — MRGN across the four confounder structures at n = 670; writes its own two files |
| `selection_metrics.R` | confounder-selection scoring: per-trio precision/recall for CS-q and CS-α, pool width, the 2×2 selection table. Writes nothing |
| `make_figures.R` | **standalone** — Figures 3 and 4 reproduced on the revised simulation; writes PNGs to `../../reports/figures/` |
| `make_selection_report.R` | **driver** — writes `../../reports/CONFOUNDER_SELECTION.md`. Run `make_figures.R` first |
| `make_inference_report.R` | **driver** — writes `../../reports/INFERENCE_PERFORMANCE.md`, every matrix in both views. Run `make_all_tables.R` first |
| `make_structure_report.R` | **driver** — writes `../../reports/CONFOUNDER_STRUCTURE.md`: which of `W`/`Z` damages MRGN, and which generating models it hits. Needs the three `data_structures/` runs |
| `compute_time.R` | **standalone** — median and IQR of per-trio compute time for all four methods. Writes `../tables/compute_time{,_long}.csv`, three PNGs, and `../../reports/COMPUTE_TIME.md`. Needs ggplot2 |
| `make_edge_pr_figure.R` | **standalone** — edge precision/recall for all four methods in a sample-size × effect-treatment grid, bootstrapped. Writes `../../reports/figures/fig_edge_pr_grid.png` (embedded in `INFERENCE_PERFORMANCE.md` §2) and `../tables/edge_pr_grid.csv`. Run `make_all_tables.R` first; needs ggplot2 |

Each method script runs standalone and leaves its counts in a
`<method>.confusion.long` global, but **writes nothing**. Only the drivers write.
`make_all_tables.R` writes three files; `confusion_structures.R` is separate — it reads
from `../data_structures/`, is not part of the main report, and writes
`../tables/confusion_structures.md` and `../tables/structure_comparison.csv`.

Two things worth knowing before reading the output:

- **MRPC's `Failed` column is a result, not missing data.** Each fit is capped at
  `mrpc.timeout`; a trio that expires has no model. At the previous 120 s cap that was 182
  of 300 trios at n = 670. Timed-out trios count against accuracy, the same way MRGN's
  `Other` and MR-GGI's `Weak instrument` do.
- **MR-GGI's four arms give one edge table, not four.** Its raw-p call is arm-invariant by
  construction, and `confusion_mrggi.R` asserts it rather than assuming it — if the
  assertion ever fires, `X` has stopped lining up with the columns of `y` in
  `mrggi.one.trio()`. The arms differ only at the `edge.fdr` level.

**No package dependencies.** Pure base R — no MRGN, no ggplot2, no `knitr`/`xtable`. In
particular `confusion_utils.R` deliberately does *not* source `../inference_utils.R`,
which pulls in `adapted_GMAC_func/*` at its top level; the single helper needed from it
(`coarse.model()`, `inference_utils.R:32`) is re-declared instead. The tables can therefore
be regenerated without a working MRGN install.

## Running

From the **repository root**:

```r
setwd("path/to/MRGN_Manuscript_Repo")
source("bioinfo_revision/simulation_results/results_scripts/make_all_tables.R")
```

or

```
Rscript bioinfo_revision/simulation_results/results_scripts/make_all_tables.R
```

Takes a couple of seconds. It reads the two `.RData` files, so run the inference stage
first (see [`../README.md`](../README.md)). Output files are overwritten in place.

## What gets written

```
../tables/
├── confusion_matrices.md      the pooled matrices, rendered for reading
├── confusion_counts_long.csv  every cell of every matrix, tidy
└── edge_comparison.csv        MRGN vs GMAC vs MR-GGI on the T1–T2 edge
```

**Three files, and `confusion_counts_long.csv` is the one to parse.** Its columns are
`method`, `arm`, `level`, `sample_size`, `effect_size`, `truth`, `predicted`, `n`, with
`effect_size = "all"` marking the pooled rows. It holds *every* cell the stage computes —
both levels, all arms, all five sample sizes, pooled and split by effect size — so any
matrix can be rebuilt from it. `matrix.from.long()` in `make_all_tables.R` does exactly
that, and is what the Markdown report is rendered from.

> An earlier layout also wrote 32 per-arm CSVs under `mrgn/`, `gmac/`, `mrggi/` and
> `by_effect_size/`, each stacking the five sample-size tables under `n = <size>` caption
> rows. They were dropped: a stacked file is not one rectangle, so `read.csv` could not
> load them, and every number in them was already a row of `confusion_counts_long.csv`.
> The archived set in [`../legacy/first_pass/tables/`](../legacy/first_pass/tables/) still
> has them, from the superseded first run.

**Two levels.** `level = "model"` is the six-way M0–M4/`Other` call; `level = "edge"` is
the same trios re-scored on whether a T1–T2 edge was found. MRGN has both, GMAC only the
second. Filter on `level` before aggregating — a `(method, arm, sample_size, effect_size)`
group holds both cross-tabs and summing over it double-counts every trio.

## Reading the tables

**Columns are the generating model, rows the inferred label** — the layout of the
manuscript tables in `Manuscript/other/tablescraps/MRGN.GMAC.class.inference.50conf`. Each
table carries four margins:

| margin | what it is |
| --- | --- |
| `Total` column | trios given this inferred label |
| `Precision` column | of those, the share whose generating model maps to this label |
| `Total` row | trios generated under this model |
| `Recall` row | of those, the share given the label that model maps to |
| bottom-right cell | overall accuracy — every correct cell over every trio |

Row percentages are gone; precision and recall replace them. An inferred label that is
never the right answer for any model — MRGN's `Other` — has no precision and is left
**blank rather than 0**, because 0 would read as "always wrong" when the truth is that the
question does not apply. A label nothing was assigned to is blank for the same reason
(0/0).

The design is fully crossed — 5 models × 5 sample sizes × 3 effect sizes × 20 replicates =
1,500 trios per method — so every pooled column totals 60, every pooled table 300, and
every per-effect-size column 20.

Truth labels come from `TRUTH.LABEL` (`inference_utils.R:29`); topologies are
`METHODS.md` Table 2:

| row | `model` | trio structure |
| --- | --- | --- |
| `M0` | `model0` | `V1 → T1`, `T2` independent |
| `M1` | `model1` | `V1 → T1 → T2` — the cis gene mediates |
| `M2` | `model2` | `V1 → T1 ← T2` — the trans gene mediates |
| `M3` | `model3` | `V1 → T1`, `V1 → T2` |
| `M4` | `model4` | `V1 → T1 → T2`, `V1 → T2` |

**Labels are coarse.** `M0.1` and `M0.2` are both reported as `M0`, which is the collapse
behind the `mrgn.*.correct.coarse` column. Nothing is lost on the truth side: the cis gene
is always `T1`, so a truth label is never a `.2` variant in the first place.

One row is not a model:

- `Other` — the fit succeeded but matched none of the eight M-topologies.

There is no `Failed` row. A trio whose fit errored would arrive as `NA`, and `confusion()`
maps `NA` to a `Failed` level that is deliberately absent from `MRGN.LEVELS` /
`GMAC.EDGE.LEVELS`, so the run stops with a named error instead of quietly reporting a
missing call as a category. There are none at present; the 114 that used to be `NA` were
recovered by `backfill_mrgn_models.R`.

### MRGN: three arms

MRGN was fitted three times per trio against three confounder sets, and each gets its own
matrices. The diagonal is the coarse-correct count.

| arm | confounders |
| --- | --- |
| `truth` | the true confounders — an **oracle**, i.e. the ceiling MRGN could reach if selection were perfect. Not an attainable result. |
| `CSq` | the CS-q selection |
| `CSa` | the CS-alpha selection |

The gap between `truth` and the other two is the cost of confounder *selection* rather
than of MRGN, which is what `METHODS.md` §5 is about. Coarse accuracy as currently
generated:

| arm | n=50 | n=150 | n=300 | n=670 | n=1000 |
| --- | --- | --- | --- | --- | --- |
| `truth` | 14.3% | 46.3% | 63.3% | 75.3% | 79.3% |
| `CSq` | 14.7% | 22.7% | 34.3% | 51.3% | 55.3% |
| `CSa` | 0.0% | 6.7% | 36.3% | 54.0% | 55.7% |

For MRGN the `truth` arm really is a ceiling. **That is not true of MRPC**, whose oracle
arm is beaten by CS-q at every sample size — see `../../reports/INFERENCE_PERFORMANCE.md`
§3.1 and `METHODS.md` Table 12c.

### GMAC: scored on the T1–T2 edge

GMAC answers a different question from MRGN. `gmac.model.call()`
(`inference_utils.R:885`) returns one of four labels by thresholding the two mediation
p-values at `selection.alpha` (0.05):

| call | cis p | trans p |
| --- | --- | --- |
| `Cis Mediated` | significant | not |
| `Trans Mediated` | not | significant |
| `No Mediation` | not | not |
| `Undirected` | significant | significant |

**The tables collapse those four to two.** GMAC's statistic is the Wald test on the
mediator coefficient in `outcome ~ mediator + treatment + confounders`, and that regression
is symmetric in T1 and T2: under `V1 → T1 → T2` the reverse coefficient is nonzero too, so
the trans test fires whenever the cis one does. GMAC has no mechanism for orienting the
edge, and the run bears it out — across all 1,500 trios `Cis Mediated` and `Trans Mediated`
hold 2 and 3 trios, and 1,356 are `Undirected`. Splitting those cells reports noise as a
direction call.

So `gmac.edge()` pools the three edge-present calls and the table asks what GMAC can
answer: is there a T1–T2 edge? That is the framing the manuscript already used —
`Manuscript/other/tablescraps/MRGN.GMAC.class.inference.50conf` reports GMAC as
`T1 - T2 Edge Absent` / `T1 - T2 Edge Present` against the same five truth models. The
four-way call is untouched in `../data/inference_gmac.csv` if it is wanted.

Precision and recall are well defined for the two-row table because each generating model
has an edge status (`GMAC.CORRECT`):

| edge | models |
| --- | --- |
| absent | `M0`, `M3` |
| present | `M1`, `M2`, `M4` |

Recall on `M1`/`M2`/`M4` is ~1.00 at every sample size above 50 — GMAC finds a real edge
essentially always. The accuracy is spent on `M0` and `M3`, where recall runs 0.05–0.43:

| | n=50 | n=150 | n=300 | n=670 | n=1000 |
| --- | --- | --- | --- | --- | --- |
| edge accuracy | 56.0% | 66.0% | 67.7% | 73.3% | 70.0% |
| `M0` recall | 0.000 | 0.117 | 0.183 | 0.433 | 0.300 |
| `M3` recall | 0.050 | 0.200 | 0.217 | 0.267 | 0.217 |

**That is a failure of confounder selection, not of the test.** GMAC's `fdr_filter` step is
meant to strip common children and intermediates from the candidate pool, but it only tests
covariate∼SNP association and q-values across the whole trio × pool p-value matrix, so a
covariate needs p ≈ 1e-7 to be filtered. A common child of T1 and T2 carries only a
doubly-attenuated SNP signal and survives; the selection step then scores on joint
association with T1 *and* T2, which is exactly what a collider maximises. GMAC ends up
adjusting for the trio's own common child in 98% of trios at n=50 and 50% at n=670, and
conditioning on a collider manufactures the edge — at n=670, of the `M0`/`M3` trios where
it selected the common child, 65 of 65 were called edge-present; where it did not, 42 of 55
were called edge-absent. Given the true confounders instead, the same test calls `M0` and
`M3` edge-absent.

This reproduces the original simulation rather than departing from it: the old
300-dataset run selected the own common child in 51.3% of trios and called 66 of 66 such
`model0`/`model3` trios `Both`.

### MRGN vs GMAC on the T1–T2 edge

The model-level tables cannot be compared across the two methods — GMAC never names a
model. So MRGN's call is collapsed onto GMAC's two rows using the same `EDGE.CORRECT`
mapping (`mrgn.edge()`), giving the `level = "edge"` rows: same columns, same right
answers, same margins as GMAC's own table, so the two are read off each other directly.
`edge_comparison.csv` and the **T1-T2 edge** section of `confusion_matrices.md` put the
rates side by side.

MRGN's edge table has a **third row**. `Other` means `infer.trio()` fitted the trio but its
four edge indicators matched none of the eight topologies — not a wrong edge call, *no*
edge call. It is not folded into either row, because that would invent output MRGN did not
produce and there is far too much of it to hide (300 of 300 trios in the CS-alpha arm at
n=50). It counts against accuracy, and it is the structural difference from GMAC, which
always answers.

| edge accuracy | n=50 | n=150 | n=300 | n=670 | n=1000 |
| --- | --- | --- | --- | --- | --- |
| MRGN `truth` | 0.233 | 0.617 | 0.710 | 0.747 | 0.730 |
| MRGN `CSq` | 0.077 | 0.263 | 0.443 | 0.613 | 0.610 |
| MRGN `CSa` | 0.000 | 0.080 | 0.383 | 0.613 | 0.593 |
| GMAC | 0.560 | 0.660 | 0.677 | 0.733 | 0.700 |

Read that alongside the two recalls before concluding anything from it. GMAC's edge-present
recall is ~0.99 at every n above 50 against MRGN `CSq`'s 0.32–0.80, but its edge-absent
recall is *lower* than MRGN `CSq`'s at every sample size (0.350 vs 0.392 at n=670) — GMAC
scores well on accuracy largely by calling almost everything an edge, which is right 60% of
the time by construction (`M1`/`M2`/`M4` are three of the five models). MRGN's accuracy is
held down by `Other` rather than by wrong calls: its edge-absent *precision* is 0.93–0.96
from n=150 on.

### A missing arm

GMAC selects its own confounders and MRGN's `CSq`/`CSa` arms select theirs, so those three
are comparable. There is no GMAC equivalent of MRGN's `truth` arm — no run of the GMAC test
against the true confounder set — which means the table cannot separate "GMAC's test" from
"GMAC's confounder selection" the way it can for MRGN, even though selection is where the
GMAC losses demonstrably come from. Adding it is a call to `run.gmac.all()`'s per-trio path
with the trio's own `U` block; see `apply.gmac()` in `../inference_utils.R`.

## Correctness check

`confusion_mrgn.R` asserts, for each of the 15 (arm × sample size) cells, that the matrix
diagonal equals `sum(mrgn.<arm>.correct.coarse)` on the same subset. Both are the same
collapse applied to the same two columns, so this is an identity, not a tolerance. If it
ever fires, the truth/inferred alignment is wrong and the script stops rather than writing
misleading tables.

## Not covered yet

- **MRPC.** It runs to completion — the first pass produced all five group checkpoints and
  a 1,500-row `inference_mrpc.*` — but nothing here tabulates it. A `confusion_mrpc.R` is a
  copy of `confusion_mrgn.R` against the `mrpc.CSq.*` columns, with two things to handle.
  Timed-out fits carry `model = NA`, which `confusion()` maps to `"Failed"` and then stops
  on, so they need a level of their own rather than an error — in the first pass that was
  224 of 1,500 trios, not a rounding error. And `mrpc.arms` is `c("CSq")` only
  (`../inference_config.R`), so it is one arm, not three.
- **Bootstrap calls.** `mrgn.*.boot.model` is a second full label set per arm and would
  double the table count. The loop in `confusion_mrgn.R` is parameterised by column name,
  so adding it is a small change.
