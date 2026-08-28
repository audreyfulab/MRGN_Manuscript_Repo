# Revised simulation study — methods

This is the methods record for the revision: what is generated, how confounders are
selected, which methods are fitted, and how their output is scored. Every setting below is
read from the scripts rather than from prose, and each is cited to the line that sets it.
It is deliberately **methods only** — no accuracy, precision or recall figures appear here.
Those live in [`../simulation_results/tables/`](../simulation_results/tables/).

**At a glance.** The study simulates 1,500 trios — a genotype `V1`, a cis gene `T1` and a
trans gene `T2` — under a fully crossed design of 5 generating topologies × 5 sample sizes
(50, 150, 300, 670, 1000) × 3 effect-size strata × 20 replicates, from a single seed
(`set.seed(234)`). Every trio additionally carries 1–50 mutually orthogonal unobserved
confounders acting on both genes, one intermediate sitting between them, one common child
below them, and — at n = 670 only — three real GTEx clinical covariates with no effect.
The confounder coefficients are calibrated so that the variance each confounder block
explains in a gene matches what selected principal components explain in real GTEx Whole
Blood, which makes the selection problem as hard here as it is there rather than easier.

**What is done to them.** Four methods are fitted per trio — MRGN, MRPC, GMAC and MR-GGI —
each under several *arms* that differ only in which covariates the method is handed: an
oracle arm carrying the trio's true confounders, two arms carrying the output of a
selection rule (CS-q, an FDR-controlled q-value screen; CS-α, an uncorrected per-test
threshold), and, for GMAC, the set GMAC selects internally. Holding the method fixed and
varying only the covariate set is what separates the cost of **confounder selection** from
the cost of the **inference method**. Scoring happens at two levels: a five-way model-label
call, and a two-row call on whether a `T1`–`T2` edge is present, the latter being the only
question all four methods can be asked in common.

## Contents

- [1. Pipeline](#1-pipeline)
- [2. Data generation](#2-data-generation)
  - [2.1 Generating models](#21-generating-models)
  - [2.2 Covariate blocks](#22-covariate-blocks)
  - [2.3 Scenario grid and parameter settings](#23-scenario-grid-and-parameter-settings)
  - [2.4 Confounder-structure variants](#24-confounder-structure-variants)
- [3. Effect sizes on a common scale](#3-effect-sizes-on-a-common-scale)
- [4. Confounder selection](#4-confounder-selection)
  - [4.1 CS-q and CS-alpha](#41-cs-q-and-cs-alpha)
  - [4.2 GMAC's internal selection](#42-gmacs-internal-selection)
  - [4.3 The oracle arm](#43-the-oracle-arm)
- [5. Methods and their settings](#5-methods-and-their-settings)
- [6. Evaluation](#6-evaluation)
- [7. Design rationale and known limitations](#7-design-rationale-and-known-limitations)
- [8. Reproducing](#8-reproducing)

---

## 1. Pipeline

| # | stage | script | writes |
| --- | --- | --- | --- |
| 1 | real confounder (PC) effects | `pc_distribution_invest/compute_pc_dist_bounds.R` | `real_pc_effect_pools.RData`, 2 PNGs |
| 2 | real SNP effects | `pc_distribution_invest/compute_effects_snp_on_gene.R` | 2 PNGs |
| 3 | simulate 1,500 trios | `simulation/updated_data_simulation.R` | `simulation/simulated_data/simulated_trios.RData` |
| 3b | simulate the structure variants | `simulation/confounder_structure_simulation.R` | `simulated_trios_u_{only,w,z}.RData` |
| 4 | calibration check | `simulation/verify_simulation.R` | console report + 2 PNGs |
| 5 | confounder selection | `simulation_results/run_confounder_selection.R` | `data/selection_group_n*.RData`, `selection_results.*` |
| 5b | add CS-i to the cached selections | `simulation_results/backfill_csi.R` | rewrites `data/selection_group_n*.RData` |
| 6 | inference | `simulation_results/run_all_inference.R` | `data/inference_{method}.*`, `inference_results.*` |
| 6b | the CS-i arm | `simulation_results/apply_{method}_csi.R`, then `merge_csi.R` | `data/inference_{method}_csi.*`, merged into `inference_{method}.*` |
| 7 | scoring | `results_scripts/make_all_tables.R` | `tables/` (3 files) |

**Table 1. The pipeline in dependency order.** Stages 1–2 measure the real GTEx data that
sets the simulation's targets; stage 3 generates the trios; stage 4 checks what stage 3
realized against what stage 1 measured; stages 5–6 select confounders and fit the methods;
stage 7 scores. Stages 1 and 2 are independent of each other, as are 4 and 5.

Every script runs from the **repository root**, not from the folder it lives in. Stages 5
and 6 are split because selection is the expensive part and does not depend on which method
is fitted afterwards; stage 5 caches per sample-size group, and stage 6 loads the cache,
validating it against the current request so a cache built from different simulated data is
detected rather than silently reused. Within stage 6 each method checkpoints separately per
group, so re-running one method leaves the others untouched.

---

## 2. Data generation

Topology comes from `MRGN::gen.graph.skel()` and the data from
`MRGN::simData.from.graph()`, both driven by `simulate.dataset()` in
`simulation/simulation_utils.R`.

### 2.1 Generating models

| model | trio structure | truth label | `T1`–`T2` edge |
| --- | --- | --- | --- |
| model0 | `V1 → T1`, `T2` independent | M0.1 | absent |
| model1 | `V1 → T1 → T2` (cis mediates) | M1.1 | present |
| model2 | `V1 → T1 ← T2` (trans mediates) | M2.1 | present |
| model3 | `V1 → T1`, `V1 → T2` | M3 | absent |
| model4 | `V1 → T1 → T2`, `V1 → T2` | M4 | present |

**Table 2. The five generating models, the label an inference must return to be scored
correct, and the edge status used at the edge level (§6).** The simulation always puts the
cis gene in `T1`, so the truth is always the `.1` variant; M3 is symmetric in the two genes
and M4 takes a single label. The lookup is explicit as `TRUTH.LABEL`
(`inference_utils.R:29`) rather than `MRGN::convert.truth()`, which maps by sorted position
and mislabels when the input does not contain all five models.

### 2.2 Covariate blocks

| block | role | count per trio | coefficient range |
| --- | --- | --- | --- |
| `K` | known clinical covariates (`pcr`, `platform`, `sex`) | 3 at n = 670, else 0 | `c(0, 0)` |
| `U` | unobserved confounders, `U → T1` and `U → T2` | 1–50 | `c(0.05, 0.3)` |
| `W` | intermediate, `T1 → W → T2` (reversed for model2) | 1 | `c(0.05, 0.3)` |
| `Z` | common child, `T1 → Z ← T2` | 1 | `c(0.3, 0.5)` |

**Table 3. The four covariate blocks attached to every trio**, from `conf.coef.ranges` at
`simulation_utils.R:133-148`. `gen.conf.coefs()` draws each coefficient's magnitude
uniformly from the interval and flips its sign at `neg.freq = 0.5`. The ranges are **raw
slopes** — see [§3](#3-effect-sizes-on-a-common-scale).

Only `U` is a confounder in the causal sense: a common parent of both genes, and the block
a selection method is scored on recovering. `W` and `Z` are covariates a method must learn
to **reject** — conditioning on `W` blocks the mediated path, conditioning on `Z` opens a
collider — which is why the oracle arm withholds them ([§4.3](#43-the-oracle-arm)). All
four entries must stay in the list even when a block is empty, since `gen.graph.skel()`
indexes it positionally.

`U` is drawn `rmvnorm(mean = 0, sigma = I)`, so the confounders are mutually orthogonal, as
principal components are in the real data. The `K` block is deliberately null: the clinical
covariates are real observed data but carry no effect on either gene, so the n = 670 arm
differs from the others only by three extra columns handed to every method as known
confounders. They are available at n = 670 only, since `pcr`/`platform`/`sex` are observed
for exactly the 670 Whole Blood donors.

### 2.3 Scenario grid and parameter settings

| parameter | setting | how it is set |
| --- | --- | --- |
| seed | `234` | `updated_data_simulation.R:13` |
| models | model0 … model4 | grid |
| sample sizes | 50, 150, 300, 670, 1000 | grid |
| effect strata | small / medium / large | grid |
| replicates | 20 | `updated_data_simulation.R:30` |
| **total datasets** | **1,500** — 300 per sample-size group | `5 × 5 × 3 × 20` |
| `b.snp` (`V1 → T1`) | small `[0.1, 0.5]`, medium `[0.5, 1.0]`, large `[1.0, 1.5]` | `runif` within the stratum |
| `b.med` (`T1 → T2`) | small `[0.1, 0.3]`, medium `[0.3, 0.5]`, large `[0.5, 1.0]` | `runif`, drawn independently of `b.snp` |
| `minor.freq` (θ) | `seq(0.01, 0.50, 0.01)`, sampled uniformly | resampled until all 3 genotype classes appear |
| `SD` (residual σ) | 1 | fixed, following Yang et al. 2017 |
| `b0.1` (intercept) | 0 | fixed, `simulation_utils.R:121` |
| `U_n` | `sample(1:50)`, capped at `n − W_n − Z_n − 4` | drawn per dataset |
| `W_n`, `Z_n` | 1 each | fixed |
| `K_n` | 3 at n = 670, else 0 | `updated_data_simulation.R:109-110` |

**Table 4. Every generating parameter, as currently coded** (`updated_data_simulation.R:13,
30, 58-64, 70-74, 88-115`). The first five rows are the factorial grid; the rest are drawn
independently for each of its 1,500 rows, so every cell of the design spans the full range
of MAF and confounder count rather than holding them fixed. `b.snp` and `b.med` are **raw
slopes**, like the coefficient ranges of Table 3 — see
[§3](#3-effect-sizes-on-a-common-scale).

Four properties of this grid are load-bearing:

- **The effect strata are contiguous and the draws continuous.** Each stratum's upper bound
  is the next one's lower bound, and values come from `runif()` over the interval rather
  than from a coarse grid, so no value in the range is unreachable. `draw.effect.sizes()`
  errors if the strata do not join up (`simulation_utils.R:56-62`).
- **Neither `b.snp` nor `b.med` may reach 0.** A `b.snp` near zero removes the `V1 → T1`
  edge, so a model0 trio would carry no edges while still bearing the `M0.1` truth label;
  `b.med` near zero does the same to model1, model2 and model4. Both small strata are
  floored at 0.1.
- **`b.snp` and `b.med` are drawn independently within a stratum**, so a trio can pair a
  strong SNP with weak mediation or the reverse — except in the large stratum, where the
  ranges do not overlap and `b.snp` always exceeds `b.med`.
- **`SD` is not scaled by `b.med`.** Because `b.med` is drawn from the same stratum as
  `b.snp`, scaling the noise by it would hold `b.snp/SD` roughly constant across strata and
  flatten the effect-size factor entirely.

The `U_n` cap keeps the confounder count below the residual degrees of freedom. It binds
only at n = 50, where it truncates the draw at 44 — which is why adjusted `R²` is
unreliable in that group. Genotypes are drawn under Hardy-Weinberg and resampled until all
three genotype classes appear; this inflates the realized MAF at low θ and small `n`, and
the count is recorded per dataset as `n.resamples`.

### 2.4 Confounder-structure variants

The main simulation gives every trio both hazards at once — one intermediate and one common
child — so it cannot say which one drives a failure. Three further simulations isolate
them, at **n = 670 only** and **MRGN only**:

| case | structure | `W_n` | `Z_n` | seed | `filter_int_child` |
| --- | --- | --- | --- | --- | --- |
| `u_only` | confounders only | 0 | 0 | 2341 | `FALSE` |
| `u_w` | + 1 intermediate | 1 | 0 | 2342 | `TRUE` |
| `u_z` | + 1 common child | 0 | 1 | 2343 | `TRUE` |
| `u_w_z` | + both — **the main simulation** | 1 | 1 | 234 | `TRUE` |

**Table 5. The four covariate structures** (`confounder_structure_simulation.R:48-60`,
`run_structure_sims.R:45-51`). Each of the three variants generates its own 300 trios
(5 models × 3 effect sizes × 20 replicates) with effect strata, MAF, `U_n` range, residual
SD and coefficient ranges identical to Table 4; only `W_n` and `Z_n` differ, so a
difference between the resulting tables has one possible cause. Each case draws under its
own seed, so the three are independent rather than one set of trios with columns deleted.
The fourth row is not regenerated — it is the n = 670 group of the main simulation.

`filter_int_child` is off for `u_only` because there is nothing to filter: no trio in that
group contributes a `W` or `Z`, and `MRGN::get.conf.trios()` stops rather than no-ops in
that case. `select.confounders()` already catches that and falls back, so setting it
explicitly changes no result; it makes the intent visible rather than leaving the right
answer to an error handler.

---

## 3. Effect sizes on a common scale

Three scales appear in this study and they are not interchangeable. Stating the conversion
once, here, is what the rest of the document relies on.

- **Raw slope.** `b.snp`, `b.med` and every entry of `conf.coef.ranges` are regression
  coefficients passed to `simData.from.graph()`. This is the scale the settings are *set*
  on (Tables 3 and 4).
- **Correlation.** A standardized effect `b · sd(X) / sd(Y)`. For the `U` block this is
  *exactly* `cor(U_i, gene)`, because a standardized coefficient equals the marginal
  correlation when the predictors are mutually uncorrelated, and `U` is drawn orthogonal by
  construction. The same holds for real PCs.
- **`R²`, variance explained.** The comparison scale used throughout. Over an orthogonal
  block it is additive, `R² = Σ_i cor(U_i, gene)²`, so `R² ≈ U_n × E[r²]`: matching the
  per-confounder effect *and* the confounder count is the same as matching the total.

**A raw slope is not a correlation, and the difference is not a constant.** They differ by
`sd(T1)`, which exceeds 1 and grows with `U_n`, so a nominal `U` upper bound of 0.3
realizes a standardized effect well below 0.3 and one that drifts with the confounder
count. Two consequences follow, and both matter more than the nominal number does:

1. The coefficient ranges of Table 3 are **not** directly comparable to the real bounds in
   Table 6. Only the *realized* correlation and `R²` are, which is what
   `verify_simulation.R` measures.
2. The calibration is **not portable**. Change `U_n`'s range, `SD`, or `b.snp`, and the
   same interval realizes a different correlation. Rerun `verify_simulation.R` after any
   such change rather than assuming it still holds.

### Calibration targets

| quantity | target | measured on |
| --- | --- | --- |
| per-confounder effect on the cis gene, sd | 0.117 | 95,564 PC–gene correlations |
| per-confounder effect on the trans gene, sd | 0.1035 | the same |
| `R²` of the cis gene on its own confounder block | 0.412 | 3,248 trios, median |
| `R²` of the trans gene on its own confounder block | 0.307 | the same |
| SNP partial correlation with a gene | modal ±0.13, central 99% inside ±0.45 | one estimate per trio |

**Table 6. The real GTEx Whole Blood values the simulation is calibrated against**,
hard-coded as `REAL` at `verify_simulation.R:26-29` and measured by
[`../pc_distribution_invest/`](../pc_distribution_invest/). The confounder rows pool every
selected PC of every one of the 3,248 trios; the central 95% of that distribution runs to
±0.20 and the central 99% to ±0.27. The SNP row is a **partial** correlation — each gene
regressed on the genotype adjusted for that trio's PCs — so the comparable simulated
quantity is `cor(V1, T1 | U)`, not the marginal `cor(V1, T1)`.

Two properties of these targets constrain how they can be read. They are **conditional on
selection**: these are the PCs that an FDR-0.05 screen already retained, so the
distribution describes confounders strong enough to have been found, and near-null
confounders are absent from it by construction. And only the **small** `b.snp` stratum is
calibrated to GTEx; medium and large deliberately overshoot the real range so that
confounder selection has enough signal to be testable at all. Results pooled across strata
are therefore an upper bound on what the same methods would achieve on real trios.

`verify_simulation.R` is a read-only check — nothing it computes feeds back into
generation — and it reports the realized distributions against Table 6 for both genes,
together with adjusted `R²` binned by confounder count. Use its **adjusted** figures. The
`R2.T1.U` / `R2.T2.U` columns recorded per dataset are **unadjusted** and are inflated
whenever `U_n` is an appreciable fraction of `n`, so they are not comparable to the 0.412 /
0.307 above.

---

## 4. Confounder selection

Selection is a separate pipeline stage from inference, and it operates on a **pooled
covariate matrix**: every trio in a sample-size group contributes its own `U`/`W`/`Z`
columns to one shared pool, against which every trio in that group is scored
(`group.cov.pool()`, `inference_utils.R:83-88`). The `K` block is not pooled — it is passed
separately as `known.conf` and included with every trio. Because selection indexes the pool
by row, datasets are processed in groups sharing a sample size, and each group is
checkpointed. The consequences of this design are in
[§7](#7-design-rationale-and-known-limitations).

### 4.1 CS-q and CS-alpha

| setting | value |
| --- | --- |
| `selection_fdr` | 0.05 |
| `filter_fdr` | 0.1 |
| `alpha` (CS-α per-test) | 0.01 |
| `filter_int_child` | `TRUE` |
| `adjust_by` | `"all"` |
| `blocksize` | `min(500, n trios)` |

**Table 7. Confounder-selection settings** (`inference_config.R:122-126`,
`inference_utils.R:213-237`). These produce **three** confounder sets per trio:

- **CS-q** — the q-value FDR screen at 5% under `adjust_by = "all"`, taken directly from
  `MRGN::get.conf.trios()`.
- **CS-α** — no multiplicity correction, per-test α = 0.01, derived by thresholding the
  `reg.pvalues` matrix the same call already returns.
- **CS-i** — the q-value FDR screen at 5% under `adjust_by = "individual"`, likewise
  derived from `reg.pvalues`. This is the setting the published GTEx analysis actually ran
  (`GTEx/data/PC_LRNA_PC_Selection_manu.R:127`) and, as §4.2 shows, it reproduces GMAC's
  internal selection exactly.

`get.conf.trios()` is called **once per group** and serves all three settings. Everything
expensive inside it — the correlation matrix over the covariate pool, and the per-trio
per-covariate regressions that produce `reg.pvalues` — is common to all three and does not
depend on `adjust_by`, so calling it once per setting would triple the cost for nothing.

### The three settings differ in one thing: the multiplicity family

All three threshold the **same** `reg.pvalues` matrix — a 2-df F test of
`lm(covariate ~ T1 + T2)`, one cell per (trio, covariate). They differ only in the family
the correction is computed over:

| setting | `adjust_by` | family | size at n = 50 |
| --- | --- | --- | ---: |
| CS-q | `"all"` | every cell of the trio × pool matrix | ~2,394,300 |
| CS-i | `"individual"` | one family per **covariate**, across trios | 300 |
| CS-α | `"none"` | no family | 1 |

**Table 7b. The multiplicity family is the whole difference between the settings.**
`reg.pvalues` is trios × covariates, so `get.q.sig()`'s `"individual"` branch —
`apply(pvalues, 2, adjust.q, ...)` — corrects down each **column**, i.e. per covariate
across trios.

The three miss in different directions, and CS-i is not between the other two by
interpolation — it is a different question. Measured on the n = 50 group, all three
thresholding the one cached matrix:

| | CS-q | CS-i | CS-α | true confounders |
| --- | ---: | ---: | ---: | ---: |
| selected/trio (mean) | 0.49 | 2.60 | 82.27 | ~25 |
| realised p threshold | 3.0 × 10⁻⁶ | 4.8 × 10⁻⁴ | 10⁻² | — |
| rejections in the group | 146 | 720 | 24,673 | — |

**Table 7c. What each family costs, n = 50.** CS-q's correction across 2.4M tests is severe
enough that it returns far fewer covariates than a trio actually has. CS-α's uncorrected
threshold returns roughly `α × pool width` false positives by arithmetic alone
(`7,981 × 0.01 ≈ 80`, against 82.3 observed). At the smaller sample sizes CS-α therefore
hands a method more covariates than there are observations, making every such fit
rank-deficient — which is why each fit is wrapped in `safely()`. **The n = 50 CS-α arm
should be read as uninformative rather than as a measurement of how a method behaves.**

Note that **equal signal density does not imply equal threshold**, and an earlier version
of these methods drew that inference. The density genuinely is the same either way — 25 of
7,981 pooled (0.0031), and 1 of 300 per covariate column (0.0033), since the private-pool
design (§7) gives each covariate exactly one trio it truly confounds. But a BH/Storey
cutoff is `q·k/m`, and `m` differs by a factor of ~8,000, so the realised threshold differs
by 159× as tabulated above. Density says how many signals exist; the family size says which
of them clear.

If `get.conf.trios()` raises "No common child or intermediate variables detected",
`select.confounders()` retries with `filter_int_child = FALSE` and records that it did so
in the `CSq.filter_int_child` / `CSa.filter_int_child` / `CSi.filter_int_child` columns, so
the fallback is visible in the results rather than silent. All three settings come from the
one call, so the three columns always agree.

Selection results are cached as `selection_group_n<size>.RData` and validated on sample
size, dataset indices, covariate names and the settings list, so a cache built against
different simulated data is recomputed rather than reused.

### 4.2 GMAC's internal selection

GMAC never reads the CS cache (`apply_gmac.R:43` sets `needs.selection = FALSE`); it
selects inside `gmac()` in two stages:

1. **Filter, `fdr_filter = 0.1`.** For each pool covariate, the p-value of
   `lm(covariate ~ SNP)`, q-valued across the whole trio × pool matrix. A covariate stays a
   candidate only if it is *not* associated with the SNP — the step meant to strip common
   children and intermediates.
2. **Select, `fdr = 0.05`.** For each surviving covariate, an F-test of `lm(cov ~ T1 + T2)`
   across trios, with stratified q-values.

The reported test uses the `Known_sel_pool` p-values — the known confounders plus the
selected pool covariates.

#### GMAC's selection is CS-i

Step 2 is `conf.fdr()` (`adapted_GMAC_func/gmac_get_conf.R:73-101`), and it is dispatched
as `lapply(1:num_pool, conf.fdr, ...)` — a loop over **covariates**, calling `qvalue()` on
that covariate's vector of per-trio p-values. That is a per-covariate family across trios,
which is exactly what `adjust_by = "individual"` applies (Table 7b). The test is identical
too: GMAC's `summary(lm(cov ~ T1 + T2))$fstatistic` through `pf()` is the same statistic
`MRGN:::p.from.reg()` computes. Same test, same FDR 0.05, same 0.1 pre-filter, same family.

This is not an analogy. Comparing the CS-i sets against GMAC's own `gmac.selected` column,
**trio by trio across all 1,500 trios**:

| n | 50 | 150 | 300 | 670 | 1000 |
| --- | ---: | ---: | ---: | ---: | ---: |
| CS-i selected/trio | 2.60 | 4.57 | 9.55 | 19.23 | 24.01 |
| GMAC selected/trio | 2.60 | 4.57 | 9.55 | 19.23 | 24.01 |
| Jaccard(CS-i, GMAC) | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

**Table 8b. CS-i and GMAC select identical confounder sets.** Not merely the same count —
the same covariates, in every trio, at every sample size. Two implementation differences
were expected to leave a residual and do not: GMAC restricts a covariate's family to the
trios that passed its child filter and estimates π₀ on that subset, where
`get.conf.trios()` drops filtered covariates from the trio's design instead; and GMAC's
vendored `qvalue()` uses a fixed `lambda = seq(0.05, 0.95, 0.05)` where MRGN's `adjust.q()`
uses `seq(0.05, max(p), 0.05)`. On this design neither changes a single selection. The
agreement is a measurement, not a proof, and should be re-checked if the pool structure
changes.

Two consequences follow. First, **the CS-i arm puts GMAC's exact covariates in front of
another method** — currently MRGN (§5), with the machinery in place for MRPC and MR-GGI.
That is what lets the tables begin to separate GMAC's *test* from GMAC's *selection*, which
no previous arm could do, since none of them handed another method GMAC's covariates.
Second, the earlier framing of GMAC as "selecting its own confounders by a different
procedure" was wrong: it uses the same procedure MRGN ships, under a setting the simulation
had not been running.

What this does **not** do is close the gap in MRGN's accuracy — see §7, Table 15. CS-i
recovers more true confounders than CS-q at every sample size and still moves MRGN's coarse
accuracy by at most 1.4 points, and moves it *down* by 5 points at n = 50. The selection
difference is real and large; its downstream effect on this design is small.

### 4.3 The oracle arm

`ground.truth.input()` (`inference_utils.R:43-49`) builds the trio's three columns, the `K`
block, and **that trio's own `U` block only**. `W` and `Z` are deliberately excluded: an
intermediate and a collider are covariates a method should reject (§2.2), so including them
would not be a baseline but a mistake.

This arm is an **oracle** — the ceiling a method could reach if selection were perfect, not
an attainable result. The gap between it and the CS-q/CS-α arms is the cost of confounder
selection rather than of the method.

---

## 5. Methods and their settings

| method | settings | source |
| --- | --- | --- |
| MRGN | `MRGN::infer.trio()`, `use.perm = FALSE`; model read as `$Inferred.Model`. Bootstrap: 1,000 resamples with replacement, non-finite replicates dropped, edge indicators averaged and called at ≥ 0.5, relabelled by `MRGN::class.vec()` | `inference_config.R:41`, `inference_utils.R:694-761` |
| MRPC | `GV = 1`, `FDR = 0.05`, `alpha = 0.01`, `indepTest = "gaussCItest"`, `FDRcontrol = "ADDIS"`; classical correlation sufficient statistic; 180 s cap per fit via `withTimeout()`; scored on the 3×3 trio block of the fitted graph | `inference_config.R:42`, `inference_utils.R:883-927` |
| GMAC | 1,000 permutations, `nominal.p = TRUE`, `fdr = 0.05`, `fdr_filter = 0.1`; mediation called at α = 0.05; run twice per group, cis-mediator and trans-mediator | `inference_config.R:43-44`, `inference_utils.R:1015-1061` |
| MR-GGI | α = 0.05, `cor.thr = 0`, first-stage `F ≥ 10` gate on `V1 → T1`; only `T1` is instrumented, `T2` and every covariate receive a zero column | `inference_config.R:192-244`, `inference_utils.R:1339-1398` |

**Table 8. Per-method tuning parameters.** All are set in `inference_config.R` so the four
methods stay comparable; changing one in an `apply_*.R` script instead breaks that.
`run_all_inference.R` launches one process per method and splits cores between them, giving
MRPC exactly one (it is single-threaded) and dividing the rest 0.45 / 0.35 / 0.20 between
MRGN, MR-GGI and GMAC.

Three settings need their rationale recorded, because each looks like a choice and is not:

- **MR-GGI's multiplicity correction is always holm.** `mrggi.p.adjust` is set to
  `"bonferroni"`, but the package body calls `p.adjust(pval.idx, method = p.adjust.methods)`
  — with a trailing `s`, base R's vector of *all* method names — so `match.arg()` silently
  takes the first element. Read `mrggi.<arm>.FDR.T1T2` as holm-adjusted whatever is passed,
  and recompute from `mrggi.<arm>.p.T1T2` if a different correction is wanted.
- **`cor.thr = 0` is the package default and is left there deliberately.** It screens gene
  pairs by correlation before testing; with a single pair per trio it only decides whether
  that trio is tested at all, so a nonzero value would discard trios untested and cap
  recall before any inference happened.
- **`mrggi.min.F = 10` is a weak-instrument gate, not a tuning knob.** With a single
  instrument MR-GGI's p-value reduces to the instrument→outcome t-statistic, so the
  exposure's first stage cancels out of the test and nothing otherwise stops it reporting a
  near-zero p for a ratio with a near-zero denominator. F > 10 is the conventional
  Staiger–Stock rule.

### The arms

| method | arms with results | covariates each arm receives |
| --- | --- | --- |
| MRGN | `truth`, `CSq`, `CSa`, **`CSi`** | oracle set / CS-q / CS-α / CS-i |
| MRPC | `truth`, `CSq` | as above; `truth` attempted only at n ≤ 300, `CSa` disabled |
| GMAC | `gmac` only | GMAC's own selection. A `gmac.truth` oracle arm is coded but has no results — see below |
| MR-GGI | `none`, `truth`, `CSq`, `CSa`, **`CSi`** | bare trio / oracle set / CS-q / CS-α / CS-i |

**Table 9. Which covariate set each method sees.** The rows within a method differ *only*
in the covariates handed to it, which is what makes the comparison isolate selection from
inference. Every arm writes its own prefixed block of columns (`mrgn.truth.*`,
`mrpc.CSq.*`, `mrggi.CSa.*`, …), and a disabled or skipped arm still emits
its columns with the reason in `<method>.<arm>.error` — so *not attempted* stays distinct
from *attempted and did not finish*, and the results schema does not change with the
configuration.

Two arm-level caveats:

- **MRPC's `truth` arm is capped at n ≤ 300** (`mrpc.truth.max.n`) and its CS-α arm is off
  (`mrpc.arms <- c("truth", "CSq", "CSi")`). Both are budget controls, not claims that the
  arms are uninformative: MRPC's cost is bimodal in the confounder count, and above roughly
  20 confounders essentially every fit reaches the cap. The recipe for re-measuring and
  raising either threshold is documented in `inference_config.R:46-120`.
- **GMAC's oracle arm is coded but has no results.** `run.gmac.group()` builds a
  `gmac.truth.*` block by calling `apply.gmac()` against that trio's true `U` columns
  (`inference_utils.R`, the `gmac.truth` block), so the arm exists in the pipeline. It has
  never been run: the code was added on 2026-08-26, while every `gmac_group_n*.RData`
  checkpoint dates from 2026-08-23 and `inference_gmac.RData` is only a re-combine of them.
  There is therefore **no `gmac.truth.*` column in the results**, and GMAC appears in the
  tables with its own selection alone. Re-running `apply_gmac.R` with
  `--rerun-inference 1` would produce it; nothing else needs to change.
- **CS-i is scored for MRGN and MR-GGI.** The selection itself exists for every method — it
  is in all five `selection_group_n*.RData` caches. MR-GGI picked it up when the stage was
  re-run on 2026-08-27, so `mrggi.CSi.*` is fully populated across all 1,500 trios. `MRPC`
  lists `CSi` in `mrpc.arms` but has not been re-run since, so `mrpc.CSi.*` and
  `gmac.CSi.*` do not exist in the results. Every scoring
  script filters arms on column presence, so they report what is there rather than failing
  or inventing an empty arm. `apply_mrpc_csi.R`, `apply_mrggi_csi.R` and `apply_gmac_csi.R`
  are written and partially checkpointed (n ≤ 300 for MRPC and MR-GGI, n ≤ 150 for GMAC)
  if the other three are wanted later; each resumes from its own checkpoints.
- **The CS-i arm was run as a separate pass** (`apply_mrgn_csi.R`), not by re-running
  `apply_mrgn.R` with `CSi` added to the arm list. The three arms already on disk are
  unchanged by it, and re-running them would have repeated the bootstrap for all three.
  `merge_csi.R` joins the `CSi.*` and `mrgn.CSi.*` columns into the master on `dataset`,
  and asserts that every column both passes actually computed — the id columns and the
  CS-q/CS-α selection scores — agrees before it writes. The CS-i selection itself did not
  need the ~40 min per group `get.conf.trios()` call at all: `backfill_csi.R` derives it
  from the `reg.pvalues` already in the cache in ~20 s per group.
- **The MRGN CS-i arm carries no bootstrap.** It was run with `--bootstrap 0`, so
  `mrgn.CSi.boot.*` is NA where the other three arms have values. This changes no model
  call and no column the tables or figures read — the bootstrap is a second, separate label
  set that nothing in the scoring stage currently consumes — and it is the difference
  between about a minute and about six hours for the arm.
- **MR-GGI's five arms are not a confounder adjustment OF THE EDGE, and must not be read as
  one.** `MRggi()` has no covariate argument; the arms differ in which covariates ride along
  as extra columns of `y`, and the estimator is strictly pairwise, so `B.T1T2` and `p.T1T2`
  are *identical* in all five. What the covariates change is the multiplicity correction
  across each gene's pairs. MR-GGI therefore writes two edge calls: `edge`, from the raw p
  — arm-invariant, and the column comparable with MRGN and GMAC — and `edge.fdr`, from the
  adjusted p, which is the only edge column that varies by arm. `confusion_mrggi.R` asserts
  the invariance rather than assuming it.

  **Its MODEL call is a different matter and IS adjusted.** `mrggi.<arm>.model` feeds
  `MRGN::class.vec()` two pairwise TSLS indicators, which no covariate set can move, plus
  two conditional regressions that ARE adjusted for that arm's covariates — exactly as
  MRGN's are. So the model call varies by arm while the raw-p edge call does not. See
  [`../MRGGI_METHODS.md`](../MRGGI_METHODS.md) for the full derivation, the trio adaptation
  and its limitations.
- **MR-GGI's correlation screen is set to `cor.thr = 0.1`, and it removes a quarter of the
  trios — unevenly.** `MRggi()` estimates a gene pair only when `|cor(g1, g2)| > cor.thr`.
  It is not confounder selection and never enters `.TSLS()`, so for a surviving pair `B` and
  `p` are unchanged; what it changes is the multiplicity family behind `edge.fdr` and, for
  the trio, **whether the trio is analysed at all**. If `|cor(T1, T2)|` falls below the
  threshold there is no estimate, and the trio is a no-call. Measured over all 1,500 trios:

  | generating model | true T1–T2 edge | screened out |
  | --- | --- | --- |
  | model0 | **none** | **65.7%** |
  | model3 | **none** | **43.0%** |
  | model1 | yes | 8.3% |
  | model2 | yes | 5.0% |
  | model4 | yes | 6.3% |

  **Table 9b. The screen falls almost entirely on the two models with no T1–T2 edge**, which
  are the trios whose T1–T2 correlation is near zero — exactly what the screen exists to
  remove. ~26% overall, and stable across sample sizes (21.7 / 24.0 / 29.3 / 27.7 / 25.7% at
  n = 50 → 1000), so it is a property of the models rather than of power. The consequence is
  structural and must be quoted with any MR-GGI edge figure: **the method is excused from
  answering on most of its true negatives, which raises its edge precision and lowers its
  edge-absent recall by construction.** The screen is kept because it is the package working
  as designed — MR on an uncorrelated pair is not meaningful — and its cost is reported
  rather than hidden by disabling it. It was `0` in the runs before 2026-08-27, which was the
  package default and harmless while MR-GGI saw only the bare trio.

---

## 6. Evaluation

### What is recorded per trio

The inference stage writes one row per trio: the generating parameters, the truth label,
the selection scoring, and one prefixed block per method arm. Selection is scored by
`score.selection()` (`inference_utils.R:52-65`) against that trio's own `U` columns, giving
`n.selected`, `n.tp`, `n.fp`, `n.fn`, `n.fp.other.trio` (false positives borrowed from a
*different* trio's block), `has.common.child` and `has.intermediate` (whether the trio's own
`Z` or `W` was selected).

MRGN and MRPC additionally carry `correct` (exact label match) and `correct.coarse`
(`M0.1` and `M0.2` both collapse to `M0`). GMAC and MR-GGI carry **no** correctness flag —
GMAC reports a mediation call, which is not a model label, so scoring it is a cross-tab
decision left to the scoring stage. MR-GGI reports an edge call and, since 2026-08-27, a
model label as well — see Table 10.

> **"Recall" means two different things in this study.** *Selection recall* is the fraction
> of a trio's true confounders that a selection rule returns (§4). *Model recall* and *edge
> recall* are per-class rates in the confusion tables (below). They are unrelated
> quantities; the surrounding section always disambiguates which is meant.

### The two scoring levels

| level | question | rows |
| --- | --- | --- |
| `model` | which of the five topologies is this? | `M0`–`M4`, plus `Other` (and `Failed` for MRPC) |
| `edge` | is there a `T1`–`T2` edge? | `Edge Absent`, `Edge Present`, plus a no-call row |

**Table 10. The two levels every method is scored at.** MRGN, MRPC and MR-GGI have both;
**GMAC has only the edge level**, because a mediation call is not a model label. The edge
level exists so all four methods face **identical rows, columns and right answers** and can
be read off each other directly. Columns are always the five generating models
(`TRUTH.LEVELS`); only the row set differs by method.

**MR-GGI's model level is reported because it fails, and the failure is the result.**
`mrggi.<arm>.model` comes from `MRGN::class.vec()`, fed MR-GGI's two pairwise TSLS
estimates plus the instrument-gene and marginal tests. With a **single** instrument the
Wald-ratio p-value reduces to the instrument→outcome t-statistic, so those two "causal"
indicators are not independent of the two marginals: `b12` is really the `V1 → T2` test and
`b22` the `V1 → T1` test, measured identical to them on **100%** of trios. The six-vector
therefore carries four distinct tests, not six, and three of the five models become
unreachable:

| model | needs | why it cannot be formed |
| --- | --- | --- |
| **M0** | `b22 = 0` | `b22` proxies `V1 → T1`, which exists in every model; it fires on 92% of all trios |
| **M3** | `b22 = 0` | same |
| **M2** | `b12 ≠ 0` | `b12` proxies `V1 → T2`, and M2 is `V1 → T1 ← T2`, where the two are marginally independent |

Measured recall at n = 1000, CS-q arm: **M0 0.017, M2 0.000, M3 0.000** against M1 0.667 and
M4 0.750. M1 and M4 are recovered because their signatures coincide with what the degenerate
tests emit, not because their structure is identified. Read the MR-GGI model rows as a
measurement of an identification limit -- one instrument identifies one edge -- rather than
as a competitive score. Full derivation in [`../MRGGI_METHODS.md`](../MRGGI_METHODS.md)
§5.2.

The right answer at the edge level is a property of the simulation, not of a method
(`EDGE.CORRECT`, `confusion_utils.R:195-196`): `M0` and `M3` have no `T1`–`T2` edge; `M1`,
`M2` and `M4` do. GMAC's four-way mediation call collapses to these two rows by mapping
`No Mediation` to absent and everything else to present; MRGN's and MRPC's model labels
collapse by the same `EDGE.CORRECT` mapping.

### Margins

Each table carries four margins (`scored.table()`, `confusion_utils.R:300-343`):

| margin | definition |
| --- | --- |
| `Total` column | trios given this inferred label |
| `Precision` column | of those, the share whose generating model maps to this label |
| `Total` row | trios generated under this model |
| `Recall` row | of those, the share given the label that model maps to |
| bottom-right | overall accuracy — every correct cell over every trio |

**Table 11. The margins reported on every confusion table.** An inferred label that is
never the right answer for any model is left **blank rather than 0**, because 0 would read
as "always wrong" when the truth is that the question does not apply. A label nothing was
assigned to is blank for the same reason.

### No-call rows

Four rows mean "no call was produced", and they are **not** the same failure:

- **`Other`** (MRGN, MRPC) — the fit succeeded but matched none of the eight M-topologies.
- **`Failed`** (MRPC only) — the fit did not finish within `mrpc.timeout`. This is a result,
  not missing data: it is a property of the method on this design.
- **`Weak instrument`** (MR-GGI) — the first-stage F fell below `mrggi.min.F`. MR-GGI tried
  to test the edge and could not.
- **`Screened out`** (MR-GGI) — `|cor(T1, T2)|` fell at or below `mrggi.cor.thr`, so the pair
  never entered `MRggi()`'s `calc.idx`. MR-GGI never looked at this trio. Distinct from
  `Weak instrument` — one is a failed test, the other an absent one — and by far the larger
  of the two: ~26% of trios at `cor.thr = 0.1`, against a weak-instrument rate that falls
  from 70% at n = 50 to 15% at n = 1000.

All four **count against accuracy**, and none is folded into an edge row, because doing so
would invent output the method did not produce. That matters most for `Screened out`: the
screen removes 65.7% of model0 and 43.0% of model3 — the two models whose right answer *is*
`Edge Absent` — so folding it in would count a refusal to answer as a correct rejection and
flatter MR-GGI badly. GMAC has no such row: it always answers. Read the no-call rate
alongside accuracy rather than as one rate.

An `NA` model is fatal rather than a category — `confusion()` maps it to a `Failed` level
that is deliberately absent from the MRGN, GMAC and MR-GGI level sets, so the run stops
with a named error instead of quietly reporting a missing call as a class. MRPC is the one
exception, where `Failed` is a legitimate reported level.

### Outputs

`make_all_tables.R` writes three files to `tables/`:

| file | contents |
| --- | --- |
| `confusion_counts_long.csv` | every cell of every matrix — `method, arm, level, sample_size, effect_size, truth, predicted, n`. The parseable artifact; any matrix can be rebuilt from it |
| `confusion_matrices.md` | the pooled matrices rendered for reading, one per method × arm × sample size |
| `edge_comparison.csv` | all methods side by side on the edge level: accuracy, and precision/recall for each edge status |

**Table 12. What the scoring stage writes.** `effect_size = "all"` marks the pooled rows;
every matrix is also built once per effect-size stratum. Filter on `level` before
aggregating — a `(method, arm, sample_size, effect_size)` group holds both cross-tabs, and
summing over it double-counts every trio. `confusion_structures.R` is standalone and writes
two further files covering MRGN across the four structures of Table 5.

Two assertions stop the run rather than producing a misleading table: the MRGN matrix
diagonal must equal `sum(mrgn.<arm>.correct.coarse)` on the same subset for each of the 15
arm × sample-size cells, and MR-GGI's raw-p edge call must be identical across all four
arms.

---

## 7. Design rationale and known limitations

### Confounder strength is calibrated down to GTEx

The pre-revision simulation used a wider `U` interval, which realized roughly **twice the
confounding** that real Whole Blood shows on both the per-confounder correlation and the
aggregate `R²`. The current interval sits essentially on the real values (Table 6).

This matters for how selection performance should be read, and the mechanism is worth
stating precisely: **confounder strength is nearly irrelevant when the confounders are
actually adjusted for** — regressing them out removes them whatever their size — and costs
performance only in proportion to what is left *unadjusted*. So the recalibration does not
make inference easier; it reduces the bias from confounders that selection fails to find,
which is the regime the CS-q arm operates in. The published >90% selection recall is
contingent on a confounding level roughly twice the real one, and should not be quoted
alongside GTEx-calibrated numbers without that qualification.

One caveat attaches to the benchmark itself: it is **mildly circular** — the real PC effect
pool is conditional on selection, so generating from that distribution and then asking a
selection rule to recover it should in principle be easy.

A second caveat, that the simulation did not run the procedure the published GTEx analysis
applied (`adjust_by = 'individual'` against `select.confounders()`'s `'all'`), **has now
been measured rather than merely recorded**: that setting is the CS-i arm (§4.1), it is
scored for MRGN alongside CS-q and CS-α, and it turned out to be GMAC's selection rule as
well (§4.2).

**As a selection the gap is large; as a change in MRGN's inference it is about one point.**
Both halves of that need stating, because either alone misleads:

| MRGN, coarse model accuracy | n = 50 | n = 150 | n = 300 | n = 670 | n = 1000 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `truth` (oracle) | 14.3% | 46.3% | 63.3% | 75.3% | 79.3% |
| `CSq` | 14.7% | 22.7% | 34.3% | 51.3% | 55.3% |
| **`CSi`** | **9.7%** | **23.0%** | **35.7%** | **52.7%** | **56.0%** |
| `CSa` | 0.0% | 6.7% | 36.3% | 54.0% | 55.7% |
| true confounders found, CS-q | 0.0 | 0.3 | 3.5 | 13.9 | 19.0 |
| true confounders found, CS-i | 0.1 | 1.3 | 5.8 | 14.9 | 19.4 |

**Table 15. CS-i against CS-q, MRGN, all 1,500 trios.** CS-i selects roughly five times as
many covariates at n = 50 (2.60 against 0.49) and recovers more true confounders at every
sample size, but converts that into only **+0.3 to +1.4 accuracy points** at n ≥ 150 — and
it is **5 points worse at n = 50**, where the extra covariates push more fits into
rank-deficiency and the `Other` rate rises from 65.7% to 73.3%.

So the earlier estimate that the setting was "worth about 2 points at n = 670" was
approximately right about *inference*, and the mismatch is not a large correction to the
reported accuracies. What was wrong was the reasoning offered for it — that the two schemes
are equivalent because they see the same signal density (§4.1) — and the conclusion drawn
from it, that the setting could be left unreconciled. It could not: CS-i is a materially
different *selection*, it is the one the deployed analysis used, and it is the rule GMAC
applies internally, which is what makes the GMAC comparison in §4.2 possible at all.

### The covariate pool is private, not shared

Every trio contributes its own `U`/`W`/`Z` columns to the group pool, and **no covariate is
ever a candidate confounder for more than the one trio that produced it**. The real GTEx
analysis does the opposite: one modest pool of PCs that every trio draws on. Two
consequences, recorded rather than fixed:

- **Signal density is ~15× lower.** A trio's true confounders as a fraction of the
  candidates it is scored against is ~0.3% here against ~4.5% in the real analysis, and
  that fraction — not the raw test count — is what a Benjamini-Hochberg threshold responds
  to. This depresses selection recall independently of confounder strength, and the two
  interact: the private pool costs little when confounders are strong enough to clear any
  threshold, and much more at realistic strength. That is why the design looked harmless in
  the published simulation.
- **Every covariate holds exactly one role, permanently.** In real data a single PC can be
  a confounder for one trio and a common child for another, which is the entire premise of
  `filter_int_child` — it flags covariates correlated with *many* variants across the
  group. Here each `W` and `Z` correlates with precisely one variant out of 300, so the
  filter is not being tested on the problem it was designed for. This is arguably the more
  serious of the two mismatches.

Moving to a shared pool is not implemented. It would require the pool to be built per
sample-size group (natural — selection already runs per group), `W` and `Z` to come from
the shared pool too (otherwise no covariate holds different roles for different trios, and
the second point above is unaddressed), and would make trios within a group non-independent,
which affects the standard errors on per-cell rates.

### Why 20 replicates

Because the pool grows with the trio count, selection costs `O(trios × pool)` and is
therefore **quadratic in the replicate count** — halving the count quarters the bill. Going
from 50 replicates to 20 cuts a full run from roughly 20 hours to roughly 3. It is nearly
free statistically: the pooled q-value threshold scales as `1/pool_size`, so a 2.5× smaller
pool barely moves the cutoff. What shrinks is per-cell precision, not sensitivity — with 75
cells and 20 replicates the standard error is near 0.11 on a per-cell rate of 0.5, so
**read the marginals rather than individual cells**. Note this cost is a property of the
private-pool design, not of the methods: a shared pool of fixed width would make selection
linear in replicates.

### Open items

- The **cis/trans asymmetry is not reproduced.** `gen.graph.skel()` draws `T1`'s and `T2`'s
  weights from the same `conf.coef.ranges$U`, so both genes share one effect distribution,
  where real data has `R²` 0.412 cis against 0.307 trans. Fixing it means leaving
  `simData.from.graph()`.
- **The calibration is not portable** (§3): the coefficient intervals are calibrated only
  for the current `U_n`, `SD` and `b.snp`.
- **`verify_simulation.R` has not been rerun since the trios were regenerated
  (2026-08-22)**, so realized calibration figures quoted anywhere in this repository predate
  the current coefficient ranges. Rerun it before quoting any of them.
- **`conf.r.squared()` reports unadjusted `R²`** (§3), inflated at small `n`.
- **The low-MAF resampling distortion is unaddressed** (§2.3): forcing all three genotype
  classes to appear inflates the realized MAF at low θ and small `n`.
- **`filter_int_child` is being evaluated off-design** — see the private-pool discussion
  above. The n = 50 group in particular should be caveated in any analysis.
- **The `W` interval may be a copy rather than a considered choice.** It matches `U`'s
  range; raising it is the remaining lever if the intermediate needs to be more detectable,
  since `W` has one parent and so no structural ceiling on `cor(T1, W)`, whereas `Z` with
  two parents approaches `1/√2`.

---

## 8. Reproducing

Run everything from the **repository root**:

```r
setwd("path/to/MRGN_Manuscript_Repo")
source("bioinfo_revision/pc_distribution_invest/compute_pc_dist_bounds.R")   # stage 1
source("bioinfo_revision/simulation/updated_data_simulation.R")              # stage 3
source("bioinfo_revision/simulation/verify_simulation.R")                    # stage 4
```

```
Rscript bioinfo_revision/simulation_results/run_confounder_selection.R         # stage 5
Rscript bioinfo_revision/simulation_results/backfill_csi.R                     # stage 5b
Rscript bioinfo_revision/simulation_results/run_all_inference.R                # stage 6
Rscript bioinfo_revision/simulation_results/results_scripts/make_all_tables.R  # stage 7
```

**Stage 5b and 6b exist because the CS-i arm was added after the rest had run.** From a
clean start they are unnecessary: `select.confounders()` now emits CS-i alongside CS-q and
CS-α, and `mrgn.arms` / `mrpc.arms` / `mrggi.arms` all list `CSi`, so a fresh stage 5 + 6
produces the CS-i columns directly. On an existing run, 5b adds CS-i to the cached
selections from the `reg.pvalues` already in them (~20 s per group, against ~40 min to
recompute), and 6b fits the one new arm without disturbing the arms already on disk:

```
Rscript bioinfo_revision/simulation_results/apply_mrgn_csi.R  --cores 2      # stage 6b
Rscript bioinfo_revision/simulation_results/apply_mrpc_csi.R  --cores 2
Rscript bioinfo_revision/simulation_results/apply_mrggi_csi.R --cores 2
Rscript bioinfo_revision/simulation_results/apply_gmac_csi.R  --cores 2
Rscript bioinfo_revision/simulation_results/merge_csi.R
```

The four `apply_*_csi.R` runs are independent and can run concurrently; `merge_csi.R` needs
all four finished. It joins on `dataset` and stops rather than writing if any column both
passes computed disagrees.

The structure variants of Table 5 are `simulation/confounder_structure_simulation.R`
followed by `simulation_results/run_structure_sims.R`, then
`results_scripts/confusion_structures.R`.

**Requirements.** R 4.4.x with `MRGN`, `MRPC`, `MRggi` (from
`remotes::install_github("hiows/MRggi")`, which depends on `susieR`), `mvtnorm`, `ggplot2`,
`patchwork`, `gridExtra`, `R.utils`, `parallel`, and — via `MRGN` — `propagate` and
`qvalue`. Stage 6 also sources the adapted GMAC code in
[`../../adapted_GMAC_func/`](../../adapted_GMAC_func/). The scoring stage is pure base R and
needs none of these.

**Data files are not in git.** `pc_distribution_invest/data/`,
`simulation/simulated_data/` and `simulation_results/data/` hold generated `.RData` and are
untracked, as is the GTEx input `GTEx/data/data.with.PCs.WholeBlood.RData` — so stages 1
and 2 cannot be rerun from a fresh clone without obtaining that file separately.

**Useful command-line overrides** (stage 6): `--sizes 50,150` restricts to sample sizes,
`--max-per-group N` runs a smoke test, `--cores N` sets the cluster size, and
`--sim-file` / `--out-dir` / `--filter-int-child` point the whole stage at a different
simulation without editing `inference_config.R`. Give a different simulation its own
`--out-dir`: the selection cache is keyed only on sample size, so two runs over different
data sharing one directory would fight over one cache file.

---

## Related documents

| document | covers |
| --- | --- |
| [`../METHODS.md`](../METHODS.md) | the full working record, including the history of how each parameter got its value and the measurements behind each decision |
| [`../MRGGI_METHODS.md`](../MRGGI_METHODS.md) | MR-GGI in depth: the trio adaptation, what it can and cannot answer, and its known limitations |
| [`../simulation_results/legacy/first_pass/README.md`](../simulation_results/legacy/first_pass/README.md) | the superseded first inference pass, and the authoritative record of what changed between the two runs |
| [`../simulation_results/tables/`](../simulation_results/tables/) | the results — confusion matrices and the edge comparison |
