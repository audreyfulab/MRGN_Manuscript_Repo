# Revised simulation study — methodology

Working record of the simulation and inference methodology for the revision. Kept
current as the work proceeds. Numbers here are measured, not intended.

Run everything from the **repository root**.

| stage | script | output |
| --- | --- | --- |
| real-data confounder effects | `pc_distribution_invest/compute_pc_dist_bounds.R` | `real_pc_effect_pools.RData`, 2 PNGs |
| real-data SNP effects | `pc_distribution_invest/compute_effects_snp_on_gene.R` | 2 PNGs |
| simulation | `updated_data_simulation.R` | `simulated_data/simulated_trios.RData` |
| calibration check | `verify_simulation.R` | console report |
| inference | `updated_simulation_inference.R` | `simulated_data/inference_results.RData` |

Shared helpers live in `simulation_utils.R` and `inference_utils.R`, both free of
top-level side effects so they can be sourced anywhere.

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
| `b.snp` | `Uniform` over the stratum: small `[0.1, 0.3]`, medium `[0.3, 0.5]`, large `[0.5, 1.0]`, step 0.05 |
| `b.med` | same stratum as `b.snp`, drawn independently |
| `SD` (residual σ) | 1, following Yang et al. 2017 |
| `U_n` | `Uniform{1..50}`, capped at `sample.size − 6` |
| `W_n`, `Z_n` | 1 each |
| `K_n` | 3 at n = 670, else 0 |

Genotypes are drawn under Hardy-Weinberg and **resampled until all three
genotype classes appear**. At θ = 0.01 and n = 50 this has taken up to 1,577
resamples and inflates the realized MAF roughly 1.7× for nominal MAF ≤ 0.05, so
the lowest-MAF scenarios do not simulate quite what their θ says. Known and
accepted.

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
| simulated cis | 19,003 | 0.101 | 0.081 | 0.295 |
| real cis | 95,564 | 0.117 | 0.101 | 0.661 |
| simulated trans | 19,003 | 0.101 | 0.080 | 0.292 |
| real trans | 95,564 | 0.104 | 0.077 | 0.653 |

The distributions overlay closely. Two differences are visible in the plot: the
real densities dip at zero, because those PCs were *selected* at FDR 0.05 and
near-null effects are filtered out by construction, and the real tails reach
±0.66 where the simulated ones stop at ±0.30.

Adjusted `R²` of each gene on its own `U` block:

| `U_n` | simulated `R²` | real |
| --- | --- | --- |
| ≤ 10 | 0.056 | 0.084 (at 9.5 PCs) |
| 10–20 | 0.161 | 0.279 (at 18) |
| 20–30 | 0.236 | 0.378 (at 26) |
| 30–40 | 0.303 | 0.451 (at 34) |
| 40–50 | 0.358 | 0.521 (at 43) |

Overall median 0.226 cis and 0.221 trans, against 0.412 and 0.307 in real data.
**The simulation is now under-confounded by roughly a factor of two**, where
before it was over-confounded by 1.6× (`R²` 0.65). Two causes, both arithmetic:
the realized effect sd is 0.101 against 0.117, and `U_n` has median 26 against 30
selected PCs in real trios, and `R² = U_n × E[r²]`. Raising the interval to about
`c(0, 0.25)` would land the median `R²` on the real 0.41; `c(0, 0.2)` keeps the
effect *range* inside the real central 95% and errs conservative.

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

**This is a power limit, not a calibration artefact.** `W` and `Z` are downstream
of `T1`, so `|cor(V1, W)| ≤ |cor(V1, T1)|`, and the real cis-eQTL partial
correlation is ~0.13. At n = 50 the standard error of a correlation is 0.14, so
an `r` of 0.15 is undetectable as a *single* test (p ≈ 0.30), before any
multiplicity correction.

It persists after the confounder effects are reduced — measured on the current
simulation, the fraction of each block clearing the pooled threshold:

| n | median \|cor(V1,W)\| | median \|cor(V1,Z)\| | W pass | Z pass |
| --- | --- | --- | --- | --- |
| 50 | 0.111 | 0.163 | 0.000 | 0.013 |
| 150 | 0.073 | 0.111 | 0.003 | 0.087 |
| 300 | 0.059 | 0.112 | 0.007 | 0.184 |
| 670 | 0.046 | 0.101 | 0.060 | 0.297 |
| 1000 | 0.041 | 0.098 | 0.083 | 0.367 |

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

- **Confounding is about half the real level** (`R²` 0.23 vs 0.41). `c(0, 0.25)`
  would match it; `c(0, 0.2)` matches the effect *range* instead and errs
  conservative. One number, §3.
- The real effect pools are **conditional on selection** (those PCs were retained
  at FDR 0.05), so they describe confounders strong enough to have been found.
  Recovery rates should not be read as covering near-null confounders.
- The cis/trans asymmetry (§3) is not reproduced.
- `conf.coef.ranges` is a raw slope while the bound is a correlation (§3).
- `conf.r.squared()` reports unadjusted `R²` (§5).
- The low-MAF resampling distortion (§2) is unaddressed.
