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
| inference | `simulation_results/updated_simulation_inference.R` | `simulation_results/inference_group_n*.RData`, then `inference_results.RData` + `.csv` |

Shared helpers live in `simulation/simulation_utils.R` and
`simulation_results/inference_utils.R`, both free of top-level side effects so
they can be sourced anywhere. Each folder carries a `README.md` describing its
scripts, data files and figures; this document covers the methodology only.

Every measured number below was last reproduced on 2026-08-19 against the current
`simulation/simulated_data/simulated_trios.RData` (3,750 datasets, `set.seed(234)`).
The §3 and §4 tables come straight from `verify_simulation.R`; the per-model and
resampling figures in §1, §2 and §4 are separate one-off diagnostics on the same
file.

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
| `U` | unobserved confounders, `U → T1` and `U → T2` | `Uniform{1..50}` | `c(0, 0.3)` |
| `W` | intermediate, `T1 → W → T2` (reversed for model2) | 1 | `c(0.05, 0.5)` |
| `Z` | common child, `T1 → Z ← T2` | 1 | `c(1, 1.5)` |

`U` is drawn `rmvnorm(mean = 0, sigma = I)`, so the confounders are mutually
orthogonal, as principal components are in the real data.

- **The `K` block is deliberately null**, carried over verbatim from
  `Simulation/sim_data.R`. The clinical covariates are real observed data but
  carry no effect on `T1` or `T2`, so the n = 670 arm differs from the others
  only by three extra columns handed to every method as known confounders.
  Verified in the generated data, over the 750 n = 670 datasets: the `K`
  coefficients on `T1` are centred at zero (medians −0.0019 / −0.0033 / +0.0023)
  and their p-values are uniform — median 0.48 and a rejection rate at α = 0.05
  of 0.047 / 0.057 / 0.044, i.e. the nominal null rate. Individual trios do reach
  small p-values (min 1e−04 over 2,250 tests), exactly as a null should.
- **`K` is only available at n = 670**, since `pcr`/`platform`/`sex` are observed
  for exactly the 670 Whole Blood donors.

---

## 2. Scenario grid

`5 models × 5 sample sizes (50, 150, 300, 670, 1000) × 3 effect sizes × 20
replicates = 1,500 datasets`, 300 per sample-size group, one row per dataset in
`scenarios`.

The replicate count is set by the cost of **confounder selection**, not by the
cost of simulating — see §4. Cutting it from 50 to 20 leaves detection power
essentially unchanged and costs only per-cell precision: with 75 cells, 20
replicates gives a standard error near 0.11 on a per-cell accuracy of 0.5, so
read the marginals rather than individual cells.

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
genotype classes appear**. This bites only at low θ and small `n`: 430 of the
3,750 datasets needed more than one draw, the median is 1, and the worst case
took 1,156 resamples at θ = 0.01, n = 50. Where it does bite it inflates the
realized MAF — 1.29× on average for nominal θ ≤ 0.05 pooled over sample sizes,
and 1.92× for those same θ at n = 50 alone. Above θ = 0.10 the distortion is
gone (ratio 1.007 for θ ∈ (0.10, 0.25], 0.997 above that). So the lowest-MAF,
smallest-`n` scenarios do not simulate quite what their θ says. Known and
accepted.

**`b.snp` and `b.med` are drawn from separate ranges**, both indexed by the same
`effect_size` stratum, and independently within a stratum. The SNP gets the wider
range because it drives everything downstream: `cor(V1, T1)` caps `cor(V1, W)`
and `cor(V1, Z)`, since `W` and `Z` are children of the genes and correlation
multiplies along a path. Measured on the current run, common-child detection
correlates 0.46 with `b.snp` and 0.42 with `b.med`, but −0.04 with `b.snp/b.med`
— the two absolute effect sizes both matter, their ratio does not (§4). The range
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

`conf.coef.ranges$U` was `c(0.05, 0.5)`; it is now `c(0, 0.3)`.

`gen.conf.coefs()` draws `|a| ~ Uniform` over this interval and flips the sign at
`neg.freq = 0.5`. The target is the distribution of real confounder effects on
cis and trans genes in Whole Blood, measured by
`pc_distribution_invest/compute_pc_dist_bounds.R` over 3,248 trios:

| per-PC effect on the cis gene, 95,564 values | 2.5% | 50% | 97.5% | sd | max |
| --- | --- | --- | --- | --- | --- |
| standardized `b·sd(PC)/sd(Y)` | −0.203 | 0.000 | 0.200 | 0.117 | 0.661 |

That is ±0.20 at the central 95% and ±0.27 at the central 99%. The bound is set
at 0.3 rather than 0.2 because **the interval is a raw slope while the target is
a correlation** — the two differ by `sd(T1)`, which exceeds 1 and grows with
`U_n`, so a nominal 0.3 realizes a standardized effect well below 0.3 (measured
median `|r|` 0.098, sd 0.123, max 0.336 against the real 0.101 / 0.117 / 0.661).
The bound was chosen on the realized scale, which is the one that has to match;
see the caveat below.

The previous upper bound of 0.5 sat well above the real distribution on either
scale and produced an adjusted `R²(T1 | U)` of 0.676 against 0.412 in real data
— roughly 1.6× too much confounding, since `R²` is additive over an orthogonal
block (`R² = Σᵢ cor(Uᵢ, T)²`) and each trio carries ~26 confounders.

### One caveat on scale

`conf.coef.ranges` is consumed as a **raw regression slope**, while the ±0.20 /
±0.27 real bounds above are **correlations**. The two differ by `sd(T1)`, which
is not 1 and which grows with `U_n`, so the realized standardized effect comes
out well below the nominal 0.3 and drifts with the confounder count. The
nominal interval is therefore not directly comparable to the real quantiles;
what matters is the realized distribution, which `verify_simulation.R` reports
and which is tabulated next.

This also means the nominal bound is not portable: change `U_n`'s range, `SD`,
or `b.snp`, and the same 0.3 realizes a different correlation. Rerun
`verify_simulation.R` after any such change rather than assuming the calibration
holds.

### Measured result

Per-confounder standardized effect, n = 670 datasets against the real pool.
`verify_simulation.R` also writes `simulated_vs_real_conf_effects.png`, which
overlays the two densities for the cis and trans gene.

| | n | sd | median \|r\| | max \|r\| |
| --- | --- | --- | --- | --- |
| simulated cis | 8,021 | 0.122 | 0.098 | 0.320 |
| real cis | 95,564 | 0.117 | 0.101 | 0.661 |
| simulated trans | 8,021 | 0.123 | 0.098 | 0.360 |
| real trans | 95,564 | 0.104 | 0.077 | 0.653 |

![confounder effects, simulated vs real](simulation_results/simulated_vs_real_conf_effects.png)

The distributions overlay closely. Two differences are visible in the plot: the
real densities dip at zero, because those PCs were *selected* at FDR 0.05 and
near-null effects are filtered out by construction, and the real tails reach
±0.66 where the simulated ones stop near ±0.35.

Adjusted `R²` of each gene on its own `U` block:

| `U_n` | simulated `R²` | real |
| --- | --- | --- |
| ≤ 10 | 0.109 | 0.084 (at 9.5 PCs) |
| 10–20 | 0.260 | 0.279 (at 18) |
| 20–30 | 0.388 | 0.378 (at 26) |
| 30–40 | 0.468 | 0.451 (at 34) |
| 40–50 | 0.521 | 0.521 (at 43) |

![R2 vs confounder count](simulation_results/r2_vs_confounder_count.png)

The figure plots this per trio rather than binned: adjusted `R²` against the
number of confounders acting on the trio, simulated (`U_n`) against real
(selected PCs), with a median line each. The two medians sit almost on top of one
another across the whole range, which is the relationship the calibration is
really trying to reproduce — a single median `R²` would hide whether the slope is
right. Simulated points are restricted to n ≥ 300, since adjusted `R²` is noisy
once `U_n` is an appreciable fraction of `n`.

Overall median **0.367 cis against the real 0.412** — the small gap is because
`U_n` has median 26 while real trios carry a median of 30 selected PCs, and
`R² = U_n × E[r²]`. Before the change the median was 0.676.

One caveat on that before-and-after: `simulated_trios_precalibration.RData` was
generated with both the old `U` bound *and* the old shared `b.snp`/`b.med` range
of `(0.1, 1.0]`, so the drop from 0.676 to 0.373 is not attributable to the `U`
bound alone. A larger `b.snp` puts more variance into `T1` that the `U` block
does not explain, which lowers `R²` on its own. Both changes push the same way.

By sample size the median adjusted `R²` is 0.267 / 0.368 / 0.365 / 0.392 / 0.391
at n = 50 / 150 / 300 / 670 / 1000 — flat apart from n = 50, where adjusted `R²`
is unreliable because `U_n` reaches 44 against 50 observations.

### The SNP effect brackets real Whole Blood

The real reference is a **partial** correlation: `compute_effects_snp_on_gene.R`
regresses each gene on the genotype adjusted for that trio's PCs. The comparable
simulated quantity is therefore `cor(V1, T1 | U)`, not the marginal
`cor(V1, T1)` — conditioning on `U` removes ~37% of `T1`'s variance, so the two
differ substantially (median 0.364 partial against 0.276 marginal). Comparing the
marginal to the real value understates the simulated effect, which an earlier
draft of this document did.

| stratum | 2.5% | median | 97.5% |
| --- | --- | --- | --- |
| small | 0.015 | **0.146** | 0.388 |
| medium | 0.109 | 0.398 | 0.669 |
| large | 0.185 | 0.572 | 0.780 |

Real GTEx cis-eQTL partial correlations peak near ±0.13 with the central 99%
inside ±0.45. **The small stratum is the GTEx-like case** — its median of 0.146
sits essentially on the real modal value — while medium and large deliberately
extend past the real range to give the confounder-selection filter something to
detect (§4). Overall, 62.9% of simulated trios fall inside the real central 99%.

The extension has precedent: Yang et al. 2017 simulate `β₁c ~ U(0.5, 1.5)` with
no confounding on the cis gene at all, yielding `cor(L, C) ≈ 0.39`. The
manuscript should still say plainly that the medium and large strata are stronger
than GTEx, so results pooled across strata are an upper bound on what the same
methods would achieve on real trios; the small stratum is the honest estimate.

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

`select.confounders()` produces two confounder sets per trio:

- **CS-q** — q-value FDR at 5%, `adjust_by = "all"`
- **CS-α** — no multiplicity correction, per-test α = 0.01

`MRGN::get.conf.trios()` is called **once per group**, with `adjust_by = "all"`.
CS-q is its output directly; CS-α is derived by thresholding the returned
`reg.pvalues` at α, which is exactly what `adjust_by = "none"` does internally —
the package builds `sigmat = reg.pvalues < alpha` and then
`apply(sigmat, 1, which)`, so the two routes are identical.

Everything expensive inside `get.conf.trios()` — `propagate::bigcor()` over the
variant/covariate correlations and the per-trio-per-covariate regressions in
`p.from.reg()` — is common to both settings and does not depend on `adjust_by`,
so calling it twice doubled the cost for nothing. Measured at ~0.955 ms per
(trio × covariate) test, a 300-trio group with an 8,340-column pool takes ~40 min
per call.

Two consequences for the recorded results:

- `CSq.filter_int_child` and `CSa.filter_int_child` now **always agree**, since
  one call serves both settings and the `filter_int_child = FALSE` fallback is
  therefore shared. Both columns are still written, so the results schema is
  unchanged. (They agreed under the two-call version as well — the fallback is
  deterministic — but only by construction rather than by design.)
- `sel$time.seconds` is now the wall clock of the single call, not the sum of two.

Selection and GMAC both index a covariate pool by row, so datasets are processed
in groups sharing a sample size, and each group is checkpointed.

### The pooled covariate design

Each trio contributes its own private `U`/`W`/`Z` columns to a shared pool, so a
group's pool is ~8,340 columns over 300 trios (mean 27.8 confounder columns per
trio). This differs from the real analysis, where every trio is scored against
one common ~670-PC matrix, and from Yang et al. 2017, whose pool is 350 variables
for 1,000 trios. Two consequences, recorded rather than fixed:

- `filter_int_child = TRUE` q-values 300 × 8,340 ≈ 2.5M correlations together, so
  a covariate must clear `|r| ≈ 0.60` at n = 50 to survive.
- CS-α at α = 0.01 over 8,340 columns implies ~83 false-positive confounders per
  trio, which is why every fit in `run.group()` is wrapped in `safely()`.

**This is what sets the replicate count.** The pool grows with the trio count, so
selection costs `O(trios × pool) = O(replicates²)`:

| replicates | trios/group | pool | min/group | ×5 groups |
| --- | --- | --- | --- | --- |
| 50 | 750 | 20,850 | 249 | 20.7 h |
| **20 (chosen)** | **300** | **8,340** | **40** | **3.3 h** |
| 10 | 150 | 4,170 | 10 | 0.8 h |

Cutting replicates is nearly free statistically. The pooled q-value threshold is
roughly `0.2 / pool_size`, so a 2.5× smaller pool moves the z cutoff only from
4.42 to 4.22 — measured, the common-child pass rate at n = 1000 went 0.557 (50
replicates) to 0.547 (20). What shrinks is per-cell precision, not sensitivity.

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
| 50 | 0.124 | 0.114 | 0.212 | 0.273 | 0.000 | 0.073 | 0.130 | 0.400 |
| 150 | 0.094 | 0.130 | 0.178 | 0.284 | 0.040 | 0.260 | 0.273 | 0.530 |
| 300 | 0.078 | 0.118 | 0.161 | 0.269 | 0.093 | 0.360 | 0.323 | 0.617 |
| 670 | 0.059 | 0.116 | 0.173 | 0.282 | 0.177 | 0.513 | 0.423 | 0.743 |
| 1000 | 0.070 | 0.116 | 0.166 | 0.262 | 0.243 | 0.547 | 0.557 | 0.753 |

1. **Effect size.** The true `|cor(V1, W)| ≈ 0.116` and `|cor(V1, Z)| ≈ 0.262`.
   From `n = (z/r)²`, `W` needs n ≈ 286 uncorrected or ≈ 1,332 at the pooled
   threshold; `Z` needs ≈ 56 and ≈ 260. `Z` is roughly 2.3× better coupled
   because it hangs off the trio by two edges with coefficients `U(1, 1.5)`,
   while `W` hangs off one with `U(0.05, 0.5)`.
2. **Sample size.** At n = 50 the uncorrected pass rate for `W` is 0.130 against
   a 0.05 null rate — power is the binding constraint there.
3. **Multiplicity.** At n = 1000 the uncorrected rate for `W` is 0.557 but the
   pooled rate is 0.243; for `Z`, 0.753 against 0.547. At large `n` the pooled
   threshold, not power, is what binds.

The recalibration roughly doubled both effects — de-noised at n = 1000, `rW`
went 0.063 → 0.121 and `rZ` 0.130 → 0.268 against the pre-recalibration run —
and the required sample sizes fell accordingly. Widening `b.snp` is the main
lever, but the lower `U` bound pushes the same way (less `U` variance in `T1`
leaves the SNP a larger share), and the two cannot be separated from the two
saved runs. `Z` now
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
link is the bottleneck — which is why the levers that worked were the ones acting
on `cor(V1, T1)`: a wider `b.snp` directly, and a smaller `U` block indirectly, by
leaving less competing variance in `T1`. Raising `Z`'s own coefficients would not
have helped, since `cor(T1, Z)` is already near its ceiling (below).

`Z` is already at its structural ceiling — with two roughly equal parents
`cor(T1, Z) → 1/√2 = 0.707` as the coefficients grow, and `U(1, 1.5)` already
reaches 0.664. `W` is not: a single parent means `cor(T1, W) → 1`, and
`U(0.05, 0.5)` reaches only 0.325. Moving `W` to `Z`'s range would roughly triple
`cor(V1, W)` and is the remaining lever if the intermediate needs to be more
detectable. `W`'s range is the same `c(0.05, 0.5)` the `U` block used before it
was recalibrated, which may be a copy rather than a considered choice.

### The failure is conditional, not universal

Detection depends on the absolute effect sizes and on `n`, not on `b.snp`
relative to `b.med`. Measured over all 3,750 datasets, common-child detection at
the pooled threshold correlates:

| with | `cor` |
| --- | --- |
| `b.snp` | 0.458 |
| `b.med` | 0.422 |
| `b.snp / b.med` | −0.036 |

**`b.med` matters about as much as `b.snp`**, which corrects an earlier reading
of these numbers. `Z` is a common child of *both* genes, and `V1` reaches `T2`
through `b.med` in model1 and model4, so the mediation coefficient does sit on a
`V1 → T2 → Z` path — it is not confined to the `T1 → T2` edge. What the ratio
measures is the *balance* between two coefficients that both help, which is why
it carries no signal.

So the filter works where the trio's effects are strong and `n` is adequate, and
n = 50 remains the weakest case by a wide margin. Rerun the stratum breakdown
after any change to the effect ranges — it is the quickest check on whether the
filter has anything to work with.

Two per-model details, both correct behaviour rather than defects. In **model2**
the intermediate is upstream (`T2 → W → T1`), so `W` and `V1` are both parents of
`T1` and marginally independent — median `|cor(V1, W)| = 0.021` at n = 1000,
de-noising to exactly 0, and the filter can never find it. In **model3** `V1`
reaches `Z` through both genes with independently signed coefficients, so the two
paths partly cancel and the median `|cor(V1, Z)|` is 0.077 at n = 1000, the
lowest of the five models (model1 is next at 0.110, model4 highest at 0.245).
Note the cancellation shows up in the *median*, not in the RMS: model3's
de-noised `rZ` is not correspondingly low, because the cancellation is
sign-dependent and leaves a wide spread rather than a uniformly small effect.

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
> appreciable fraction of `n`. Measured on the current run: median 0.456
> unadjusted against 0.373 adjusted overall, and 0.712 against 0.321 at n = 50
> with median `U_n = 26` — an inflation of 0.39 there, regardless of the true
> value. Compare `verify_simulation.R`'s adjusted figures against the real
> 0.41 / 0.31 instead.

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
- `conf.coef.ranges` is a raw slope while the bound is a correlation (§3), so the
  nominal 0.3 is calibrated only for the current `U_n`, `SD` and `b.snp`.
- `conf.r.squared()` reports unadjusted `R²` (§5).
- The low-MAF resampling distortion (§2) is unaddressed.
- The before/after comparisons against
  `simulated_trios_precalibration.RData` (§3, §4) confound two changes — the `U`
  bound and the `b.snp`/`b.med` ranges moved in the same run. Isolating either
  would need a third run.
- `W`'s `c(0.05, 0.5)` is the range the `U` block used before recalibration and
  may be a copy rather than a considered choice (§4). It is the remaining lever
  if the intermediate needs to be more detectable.
- The n = 50 group's `CSq.filter_int_child` / `CSa.filter_int_child` columns still
  need checking on the next inference run (§4).
