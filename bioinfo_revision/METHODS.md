# Revised simulation study — methodology

Working record of the simulation and inference methodology for the revision. Kept
current as the work proceeds. Numbers here are measured, not intended.

Run everything from the **repository root**.

| stage | script | output |
| --- | --- | --- |
| real-data confounder effects | `pc_distribution_invest/compute_pc_dist_bounds.R` | `pc_distribution_invest/data/real_pc_effect_pools.RData`, 2 PNGs |
| real-data SNP effects | `pc_distribution_invest/compute_effects_snp_on_gene.R` | 2 PNGs |
| simulation | `simulation/updated_data_simulation.R` | `simulation/simulated_data/simulated_trios.RData` |
| calibration check | `simulation/verify_simulation.R` | console report + `simulation_results/simulated_vs_real_conf_effects.png` |
| inference | `simulation_results/updated_simulation_inference.R` | `simulation_results/inference_results.RData` |

Shared helpers live in `simulation/simulation_utils.R` and
`simulation_results/inference_utils.R`, both free of top-level side effects so
they can be sourced anywhere.

---

## 1. Generating models

Each dataset is one trio — genotype `V1`, cis gene `T1`, trans gene `T2` — plus
four covariate blocks. Topology comes from `MRGN::gen.graph.skel()` and the data
from `MRGN::simData.from.graph()`.

| model | trio structure | inference label |
| --- | --- | --- |
| model0 | `V1 → T1`, `T2` independent | M0.1 |
| model1 | `V1 → T1 → T2` (cis mediates) | M1.1 |
| model2 | `V1 → T1 ← T2` (trans mediates) | M2.1 |
| model3 | `V1 → T1`, `V1 → T2` | M3 |
| model4 | `V1 → T1 → T2`, `V1 → T2` | M4 |

The cis gene is always `T1`, so the truth is always the `.1` variant. The label
lookup is explicit in `updated_simulation_inference.R` rather than
`MRGN::convert.truth()`, which maps by sorted position and mislabels when the
input does not contain all five models.

### Covariate blocks

| block | role | count | `conf.coef.ranges` |
| --- | --- | --- | --- |
| `K` | known clinical covariates (`pcr`, `platform`, `sex`) | 3 at n = 670, else 0 | `c(0, 0)` |
| `U` | unobserved confounders, `U → T1` and `U → T2` | `Uniform{1..50}` | `c(0, 0.2)` |
| `W` | intermediate, `T1 → W → T2` (reversed for model2) | 1 | `c(0.05, 0.5)` |
| `Z` | common child, `T1 → Z ← T2` | 1 | `c(1, 1.5)` |

`U` is drawn `rmvnorm(mean = 0, sigma = I)`, so the confounders are mutually
orthogonal, as principal components are in the real data.

- **The `K` block is deliberately null**, carried over verbatim from
  `Simulation/sim_data.R`. The clinical covariates are real observed data but
  carry no effect on `T1` or `T2`, so the n = 670 arm differs from the others
  only by three extra columns handed to every method as known confounders.
  Verified in the generated data: `K` coefficients ≈ 0, p > 0.09.
- **`K` is only available at n = 670**, since `pcr`/`platform`/`sex` are observed
  for exactly the 670 Whole Blood donors.

---

## 2. Scenario grid

`5 models × 5 sample sizes (50, 150, 300, 670, 1000) × 3 effect sizes × 50
replicates = 3,750 datasets`, one row per dataset in `scenarios`.

| parameter | draw |
| --- | --- |
| `minor.freq` (θ) | `Uniform{0.01, 0.02, …, 0.50}` |
| `b.snp` | stratum tertile of `(0, 1.5]`: small `[0.05, 0.50]`, medium `[0.55, 1.00]`, large `[1.05, 1.50]`, step 0.05 |
| `b.med` | stratum tertile of `(0, 1.0]`: small `[0.05, 0.35]`, medium `[0.40, 0.70]`, large `[0.75, 1.00]`, step 0.05 |
| `SD` (residual σ) | 1, following Yang et al. 2017 |
| `U_n` | `Uniform{1..50}`, capped at `sample.size − 6` |
| `W_n`, `Z_n` | 1 each |
| `K_n` | 3 at n = 670, else 0 |

Genotypes are drawn under Hardy-Weinberg and **resampled until all three
genotype classes appear**. At θ = 0.01 and n = 50 this has taken up to 1,577
resamples and inflates the realized MAF roughly 1.7× for nominal MAF ≤ 0.05, so
the lowest-MAF scenarios do not simulate quite what their θ says. Known and
accepted.

**`b.snp` and `b.med` are drawn from separate ranges**, both indexed by the same
`effect_size` stratum, and independently within a stratum. The SNP gets the wider
range because it drives everything downstream: `cor(V1, T1)` caps `cor(V1, W)`
and `cor(V1, Z)`, since `W` and `Z` are children of the genes and correlation
multiplies along a path. Measured on the previous run, common-child detection
correlates 0.52 with `b.snp` but only 0.04 with `b.snp/b.med` — absolute SNP
size matters, its size relative to the mediation effect does not. The range
`(0, 1.5]` also matches Yang et al. 2017, who drew `β₁c ~ U(0.5, 1.5)`.

**Neither range may include 0.** `b.snp = 0` removes the `V1 → T1` edge, so a
"model0" trio would carry no edges and its `M0.1` truth label would be wrong;
`b.med = 0` does the same to model1, model2 and model4. `draw.effect.sizes()`
errors on a zero lower bound. One consequence of the asymmetric ranges: in the
large stratum `b.snp` always exceeds `b.med` (1.05–1.50 vs 0.75–1.00), so that
cell cannot produce strong mediation with a weak SNP.

**`SD` must not be scaled by `b.med`.** `Simulation/sim_data.R` used
`SD = sample(seq(0.3, 1.5, 0.1)) * b.med`, but `b.med` here is drawn from the same
stratum as `b.snp`, so scaling the noise by it holds `b.snp/SD` roughly constant
across strata (measured 1.06 / 1.38 / 1.22) and the effect-size factor goes flat
and non-monotone — partial `|cor(V1,T1)|` of 0.55 / 0.63 / 0.48 for small /
medium / large. The original script had no effect-size strata, so its formula
does not transfer.

---

## 3. Confounder effect size

### The change

`conf.coef.ranges$U` was `c(0.05, 0.5)`; it is now `c(0, 0.2)`.

`gen.conf.coefs()` draws `|a| ~ Uniform` over this interval and flips the sign at
`neg.freq = 0.5`. The bound of 0.2 is the central 95% of the real confounder
effects on cis and trans genes in Whole Blood, measured by
`pc_distribution_invest/compute_pc_dist_bounds.R` over 3,248 trios:

| per-PC effect on the cis gene, 95,564 values | 2.5% | 50% | 97.5% | sd | max |
| --- | --- | --- | --- | --- | --- |
| standardized `b·sd(PC)/sd(Y)` | −0.203 | 0.000 | 0.200 | 0.117 | 0.661 |

The previous upper bound of 0.5 sat well above the real 97.5% point and produced
`R²(T1 | U) = 0.65` against 0.41 in real data — roughly 1.9× too much
confounding, since `R²` is additive over an orthogonal block
(`R² = Σᵢ cor(Uᵢ, T)²`) and each trio carries ~26 confounders.

### One caveat on scale

`conf.coef.ranges` is consumed as a **raw regression slope**, while the 0.2 above
is a **correlation**. The two differ by `sd(T1)`, which is not 1 and which grows
with `U_n`, so the realized standardized effect comes out somewhat below 0.2 and
drifts with the confounder count. The interval is therefore a defensible bound
rather than an exact match; `verify_simulation.R` reports what is realized.

### Measured result

Per-confounder standardized effect, n = 670 datasets against the real pool.
`verify_simulation.R` also writes `simulated_vs_real_conf_effects.png`, which
overlays the two densities for the cis and trans gene.

| | n | sd | median \|r\| | max \|r\| |
| --- | --- | --- | --- | --- |
| simulated cis | 19,591 | 0.123 | 0.098 | 0.336 |
| real cis | 95,564 | 0.117 | 0.101 | 0.661 |
| simulated trans | 19,591 | 0.123 | 0.097 | 0.362 |
| real trans | 95,564 | 0.104 | 0.077 | 0.653 |

The distributions overlay closely. Two differences are visible in the plot: the
real densities dip at zero, because those PCs were *selected* at FDR 0.05 and
near-null effects are filtered out by construction, and the real tails reach
±0.66 where the simulated ones stop near ±0.35.

Adjusted `R²` of each gene on its own `U` block:

| `U_n` | simulated `R²` | real |
| --- | --- | --- |
| ≤ 10 | 0.107 | 0.084 (at 9.5 PCs) |
| 10–20 | 0.267 | 0.279 (at 18) |
| 20–30 | 0.382 | 0.378 (at 26) |
| 30–40 | 0.460 | 0.451 (at 34) |
| 40–50 | 0.523 | 0.521 (at 43) |

The curve tracks the real one closely across the whole confounder-count range.
Overall median **0.373 cis against the real 0.412** — the small gap is because
`U_n` has median 26 while real trios carry a median of 30 selected PCs, and
`R² = U_n × E[r²]`. Before the change the median was 0.65.

By sample size the median adjusted `R²` is 0.321 / 0.379 / 0.378 / 0.374 / 0.388
at n = 50 / 150 / 300 / 670 / 1000 — flat, as it should be, since the confounding
is a property of the generating model rather than of `n`.

### The SNP effect now exceeds real Whole Blood

Median `|cor(V1, T1)| = 0.286`, separating across the strata at **0.121 / 0.308 /
0.473**. Real GTEx cis-eQTL partial correlations peak at ±0.13 with the central
99% inside ±0.45, so the medium and large strata sit above what is observed, and
the large stratum's median exceeds the real 99th percentile.

This is a deliberate trade, made to give the confounder-selection filter enough
signal to be testable (§4), and it has precedent — Yang et al. 2017 simulate
`β₁c ~ U(0.5, 1.5)` with no confounding on the cis gene at all, which yields
`cor(L, C) ≈ 0.39`. It should nonetheless be stated plainly in the manuscript
that the simulated SNP effects are stronger than GTEx's, so method performance
here is an upper bound on what the same methods would achieve on real trios.

### Known limitation

`gen.graph.skel()` draws `T1`'s and `T2`'s weights from the same
`conf.coef.ranges$U`, so both genes share one effect distribution. Real data has
`R²` 0.412 cis versus 0.307 trans; the simulation cannot reproduce that asymmetry
without leaving `simData.from.graph()`.

---

## 4. Confounder selection and inference

`updated_simulation_inference.R` runs six inferences per trio:

| method | confounders |
| --- | --- |
| MRGN | the true ones (trio + `K` + that trio's own `U`) |
| MRGN | CS-q selected |
| MRGN | CS-α selected |
| MRPC | CS-q selected |
| MRPC | CS-α selected |
| GMAC | whatever GMAC selects for itself |

`MRGN::get.conf.trios()` runs in two settings via `select.confounders()`:

- **CS-q** — q-value FDR at 5%, `adjust_by = "all"`
- **CS-α** — no multiplicity correction, per-test α = 0.01, `adjust_by = "none"`

Selection and GMAC both index a covariate pool by row, so datasets are processed
in groups sharing a sample size, and each group is checkpointed.

### The pooled covariate design

Each trio contributes its own private `U`/`W`/`Z` columns to a shared pool, so a
group's pool is ~21,000 columns over 750 trios. This differs from the real
analysis, where every trio is scored against one common ~670-PC matrix. Two
consequences, recorded rather than fixed:

- `filter_int_child = TRUE` q-values 750 × 21,000 ≈ 15.8M correlations together,
  so a covariate must clear `|r| ≈ 0.55` at n = 50 to survive.
- CS-α at α = 0.01 over 21,000 columns implies ~210 false-positive confounders
  per trio, which is why every fit in `run.group()` is wrapped in `safely()`.

### Confounder selection has no power at n = 50

`get.conf.trios()` finds no intermediate or common child variables in the n = 50
group and raises "No common child or intermediate variables detected";
`select.confounders()` catches this and falls back to `filter_int_child = FALSE`,
recording it in the `CSq.filter_int_child` / `CSa.filter_int_child` columns.

**This is not a calibration artefact — it persists after the confounder effects
are corrected.** There are three distinct obstacles, and which one binds depends
on `n`.

`verify_simulation.R` separates them. The observed correlations *fall* with `n`
because a sample correlation carries about `1/√n` of noise; de-noising with
`√(mean(r²) − 1/n)` recovers a true effect that is constant across `n`, which is
the point:

| n | med \|rVW\| | true rW | med \|rVZ\| | true rZ | W pooled | Z pooled | W α=.05 | Z α=.05 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 50 | 0.135 | 0.128 | 0.207 | 0.269 | 0.000 | 0.063 | 0.147 | 0.380 |
| 150 | 0.098 | 0.134 | 0.183 | 0.270 | 0.028 | 0.240 | 0.288 | 0.537 |
| 300 | 0.081 | 0.119 | 0.167 | 0.268 | 0.073 | 0.331 | 0.352 | 0.637 |
| 670 | 0.067 | 0.124 | 0.164 | 0.276 | 0.164 | 0.491 | 0.465 | 0.717 |
| 1000 | 0.061 | 0.121 | 0.161 | 0.268 | 0.233 | 0.557 | 0.497 | 0.752 |

1. **Effect size.** The true `|cor(V1, W)| ≈ 0.121` and `|cor(V1, Z)| ≈ 0.268`.
   From `n = (z/r)²`, `W` needs n ≈ 262 uncorrected or ≈ 1,337 at the pooled
   threshold; `Z` needs ≈ 54 and ≈ 273. `Z` is roughly 2.2× better coupled
   because it hangs off the trio by two edges with coefficients `U(1, 1.5)`,
   while `W` hangs off one with `U(0.05, 0.5)`.
2. **Sample size.** At n = 50 the uncorrected pass rate for `W` is 0.147 against
   a 0.05 null rate — power is the binding constraint there.
3. **Multiplicity.** At n = 1000 the uncorrected rate for `W` is 0.497 but the
   pooled rate is 0.233; for `Z`, 0.752 against 0.557. At large `n` the pooled
   threshold, not power, is what binds.

Widening `b.snp` to `(0, 1.5]` roughly doubled both effects (`rW` 0.068 → 0.121,
`rZ` 0.163 → 0.268) and the required sample sizes fell accordingly — `Z` now
clears the pooled threshold from about n = 273, so the n = 300, 670 and 1000
groups detect common children in 33%, 49% and 56% of trios.

**The n = 50 group no longer returns nothing.** 6.3% of its trios now have a
detectable common child, so `get.conf.trios()` should no longer raise "No common
child or intermediate variables detected" and `select.confounders()` should no
longer fall back to `filter_int_child = FALSE`. Confirm this on the next
inference run via the `CSq.filter_int_child` / `CSa.filter_int_child` columns.
Detection at n = 50 is still poor in absolute terms, so the group's confounder
sets will contain colliders for most trios and should be caveated.

Both blocks are capped above by `|cor(V1, T1)|`, since `W` and `Z` are downstream
of `T1` and correlation multiplies along a path. Measured at n = 1000:

| link | median \|r\| |
| --- | --- |
| `cor(V1, T1)` variant → cis gene | 0.280 |
| `cor(T1, W)` | 0.325 |
| `cor(T1, Z)` | 0.664 |

`0.280 × 0.664 = 0.186` against an observed `cor(V1, Z)` of 0.161. `Z`'s effect on
the *genes* is large, but the filter tests it against the *variant*, and the first
link is the bottleneck — which is why widening `b.snp` was the effective lever and
recalibrating the confounders was not.

`Z` is already at its structural ceiling — with two roughly equal parents
`cor(T1, Z) → 1/√2 = 0.707` as the coefficients grow, and `U(1, 1.5)` already
reaches 0.664. `W` is not: a single parent means `cor(T1, W) → 1`, and
`U(0.05, 0.5)` reaches only 0.325. Moving `W` to `Z`'s range would roughly triple
`cor(V1, W)` and is the remaining lever if the intermediate needs to be more
detectable. `W`'s range is the same `c(0.05, 0.5)` the `U` block used before it
was recalibrated, which may be a copy rather than a considered choice.

### The failure is conditional, not universal

Detection depends on `b.snp` and `n`, not on the model or on `b.snp` relative to
`b.med`. The ratio is irrelevant — `cor(Z pass, b.snp/b.med) = 0.036`, against
`cor(Z pass, b.snp) = 0.523` — which follows from `b.med` being the `T1 → T2`
edge, not on the path from `V1` to `W` or `Z`.

So the filter works where the SNP effect is strong and `n` is adequate, and n = 50
remains the weakest case by a wide margin. Rerun the stratum breakdown after any
change to the effect ranges — it is the quickest check on whether the filter has
anything to work with.

Two per-model details, both correct behaviour rather than defects: in **model2**
the intermediate is upstream (`T2 → W → T1`), so `W` and `V1` are both parents of
`T1` and marginally independent — `cor(V1, W) = 0.022`, and the filter can never
find it. In **model3** `V1` reaches `Z` through both genes with independently
signed coefficients, so the two paths partly cancel and `cor(V1, Z) = 0.046`,
the lowest of the five models.

There is also a design mismatch worth recording. `get.conf.trios()` identifies an
intermediate or common child by its correlation with variants **across the whole
group**: in real GTEx every trio is scored against the same PC matrix, so a PC
that is a common child for many trios correlates with many variants and stands
out. Here each trio owns a private `W` and `Z` correlated with exactly one
variant out of 750, which is not the setting the filter was designed for.

The consequence is material: with the filter off, the common child `Z` is
selected as a confounder for essentially every n = 50 trio (measured on the
pre-revision data: median p = 8e−10 against `T1` or `T2`; `Z` beat every true `U`
in 96.5% of trios). MRGN and MRPC then condition on a collider. **The n = 50
group should be caveated or excluded in the analysis.**

---

## 5. Diagnostics recorded per dataset

`params` carries every generating parameter plus `n.resamples` and
`R2.T1.U` / `R2.T2.U`.

> `R2.T1.U` / `R2.T2.U` are **unadjusted** `R²`, inflated whenever `U_n` is an
> appreciable fraction of `n` — at n = 50 with `U_n = 26` the inflation is about
> 0.35 regardless of the true value. Compare `verify_simulation.R`'s adjusted
> figures against the real 0.41 / 0.31 instead.

`conf.effects` holds the cis and trans regressions of each gene on all
covariates. These condition on the collider `Z` and the mediator `W` and omit
`V1`, so they are not comparable to `compute_pc_dist_bounds.R`, which regresses
on PCs only.

---

## 6. Open items

- The real effect pools are **conditional on selection** (those PCs were retained
  at FDR 0.05), so they describe confounders strong enough to have been found.
  Recovery rates should not be read as covering near-null confounders.
- The cis/trans asymmetry (§3) is not reproduced.
- `conf.coef.ranges` is a raw slope while the bound is a correlation (§3).
- `conf.r.squared()` reports unadjusted `R²` (§5).
- The low-MAF resampling distortion (§2) is unaddressed.
