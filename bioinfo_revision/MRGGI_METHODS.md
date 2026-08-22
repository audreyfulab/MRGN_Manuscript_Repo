# MR-GGI on trios — methodology

How MR-GGI (https://github.com/hiows/MRggi) is applied to the simulated trios, why two
adaptations are needed, and what it can and cannot answer. Companion to
[`METHODS.md`](METHODS.md), which covers the simulation and the other three methods.

Every number here is measured. The reproducible checks live in
`simulation_results/mrggi_feasibility.R`; the pipeline code is
`simulation_results/apply_mrggi.R` with `run.mrggi.group()` in `inference_utils.R`.

---

## 1. What MR-GGI estimates

MR-GGI is a **Mendelian randomisation** method. Where MRGN and MRPC decide the T1–T2
edge by conditioning on a set of confounders, MR-GGI uses a genetic variant as an
**instrument** and never adjusts for anything.

For a pair of genes it computes a two-stage least squares estimate. From the installed
source (`MRggi:::.TSLS`):

```r
Bzx_fit  = summary(lm(y1 ~ X1))       # instrument -> exposure   (first stage)
y2.resid = resid(lm(y2 ~ X2))         # strip the outcome's own instruments
Bzy_fit  = summary(lm(y2.resid ~ X1)) # instrument -> outcome    (reduced form)
Bxy      = sum(Bzx*Bzy*IV) / sum(Bzx*Bzx*IV)
```

With one instrument this is the **Wald ratio**:

```
    Bxy  =  beta(V -> outcome) / beta(V -> exposure)
```

The logic is that a genotype is randomised at conception, so `V1` is independent of the
confounders `U`. The part of `T1` that `V1` explains therefore carries no `U`, and an
effect estimated through that part is unconfounded — without `U` ever being measured.

| estimator | mean | bias | sd |
| --- | --- | --- | --- |
| Wald ratio, raw genes | 0.6921 | **−0.0079** | 0.0940 |
| Wald ratio, confounders residualised out first | 0.6974 | −0.0026 | 0.0845 |
| naive `lm(T2 ~ T1)`, no instrument | 0.8419 | **+0.1419** | 0.0510 |

**Table 1. The instrument removes confounding bias that ordinary regression cannot.**
400 replicates, n = 500, one variant, ten confounders acting on both genes, true
`T1 -> T2` effect 0.700. The naive regression overstates the effect by 20% because the
confounders induce association between the genes; the Wald ratio is unbiased. Note
residualising the confounders out beforehand does **not** reduce bias — it reduces the
standard error by 18%, because `V1` is independent of `U` so removing `U` cuts residual
noise without moving either coefficient in the ratio. See §4.

---

## 2. Why a trio needs an adaptation

**MR-GGI assumes every gene carries its own distinct cis-instrument. A trio has one
variant shared by two genes.** That difference breaks the estimator outright.

The `y2.resid = resid(lm(y2 ~ X2))` line strips the *outcome* gene's own instruments
before the reduced-form regression. OLS residuals are orthogonal to their own regressors,
so when the exposure and outcome are handed the **same** instrument, `Bzy` is forced to
exactly zero and the ratio collapses.

| instruments supplied | Bzx | Bzy | **Bxy** |
| --- | --- | --- | --- |
| both genes given `V1` — the trio case | 0.787 | **0.000000** | **0.00000** |
| each gene its own SNP — MR-GGI's intended case | 0.787 | 0.532 | 0.67652 |
| exposure `V1`, outcome none | 0.787 | 0.530 | **0.67397** |

**Table 2. The collapse, on a synthetic `V1 -> T1 -> T2` trio with a true effect of
0.700.** `cor(y2.resid, V1) = -2.24e-16` in the first row — machine zero, not a numerical
accident. The failure is silent: it returns `0.0000` with `p = 1.0000`, which reads as a
confident absence of an edge.

This is not a defect in the package or in the simulation. It is what the two designs
assume:

| | MR-GGI's intended setting | our trio |
| --- | --- | --- |
| genotypes | one per gene (`p_geno = p_pheno`) | 1 |
| genes | many | 2 |
| instrument for the exposure | that gene's cis-SNP | `V1` |
| instrument for the outcome | a **different** cis-SNP | `V1` — the same |

**Table 3. Why the reference script works and a trio does not.** The example MR-GGI
script sets `p_geno = p_pheno = n_t`, giving every gene its own variant, so
`FineMapping()` hands each gene a different instrument and the outcome-side
residualisation removes that gene's own cis-regulation as intended. A trio's *trans* gene
has no cis-variant by construction — that is what makes it a trio rather than a network.

Two consequences worth recording:

- **`FineMapping()` does trigger this**, confirmed rather than assumed. On a trio it
  assigns `V1` to *both* genes (one column each), and running `MRggi()` on its output
  returns `Bg1g2 = 0.0000, p = 1.0000`. It is bypassed for that reason.
- **Pooling would not help.** Running per sample-size group with 300 genotypes and 600
  genes does not fix it: within a trio, `V1` is still the only variant associated with
  either gene, so susie selects it for both.

---

## 3. The adaptation

**Give the exposure its instrument and give the outcome none.** The residualisation then
has nothing to remove, and `MRggi()` computes the standard one-sample Wald ratio —
0.674 against a true 0.700 in Table 2. This is still MR-GGI's own estimator, fed
instruments appropriate to a trio, and it is the textbook single-instrument MR setup.

The "no instrument" entry has to be a **column of zeros**. `MRggi()` runs
`lapply(X, scale)` over every element before its main loop, so:

| outcome-gene entry | result |
| --- | --- |
| `NULL` | error: `'data' must be of a vector type, was 'NULL'` |
| zero-column matrix | error: `subscript out of bounds` |
| **column of zeros** | **works** — `Bg1g2 = 0.674`, reverse direction `NaN` |

**Table 4. What can stand in for "this gene has no instrument".** The `NaN` in the
reverse direction is correct behaviour, not a failure: with no instrument for `T2`, the
`T2 -> T1` effect is not estimable.

Because each call yields only the direction whose exposure holds the instrument, each trio
needs two calls with the roles swapped. Both are stored; only the forward one is used for
edge calls, for the reason in §6.

---

## 4. Confounding — why nothing is adjusted for

`MRggi()` takes no covariate argument, and the estimation is strictly pairwise: the only
data entering the regressions are `y[,g1]`, `y[,g2]` and their instruments. **Adding
confounders to `y` therefore changes the T1–T2 estimate not at all** — it only adds more
gene pairs to the output, which would be discarded.

This differs from MRPC, where confounders as graph nodes genuinely do adjust, because the
PC algorithm's conditional independence tests condition on subsets of the other nodes.
The same phrase — "confounders in the model" — means different things for the two methods.

So the only real option would be residualising `T1` and `T2` on a confounder set before
inference. **It is not done**, on the reasoning that the instrument is what is supposed to
handle confounding, and our `U` block is precisely the mediator–outcome confounding MR
exists to defeat. Table 1 supports this: unadjusted, the Wald ratio's bias is −0.008
against the naive regression's +0.142. Residualising buys an 18% reduction in standard
error — a precision gain, not a bias correction — and would hand MR-GGI an advantage the
other methods do not get for free.

---

## 5. What MR-GGI can and cannot answer

MR-GGI returns directed edges between genes. For a trio that is the T1–T2 edge only. It
never estimates the `V1 -> T1` or `V1 -> T2` edges, so it **cannot** produce an M0–M4
model label:

| generating model | true T1–T2 edge | what MR-GGI can say |
| --- | --- | --- |
| model0 (`V1 -> T1`) | none | no edge |
| model3 (`V1 -> T1`, `V1 -> T2`) | none | no edge |
| model1 (`V1 -> T1 -> T2`) | `T1 -> T2` | `T1 -> T2` |
| model4 (`V1 -> T1 -> T2`, `V1 -> T2`) | `T1 -> T2` | `T1 -> T2` |
| model2 (`V1 -> T1 <- T2`) | `T2 -> T1` | `T2 -> T1` |

**Table 5. Five generating models collapse to three possible answers.** model0 and model3
are indistinguishable to MR-GGI, as are model1 and model4, because each pair differs only
in the `V1 -> T2` edge it never estimates.

It is therefore scored on the **T1–T2 edge**, a metric the manuscript already uses
(`RecallT1`, `MRGN_v8.pdf` line 377), and carries no `correct` flag in the results — the
same convention as GMAC, whose mediation call is likewise not a model label. Cross-tabbing
against the five models belongs in the analysis stage.

---

## 6. Known limitations

### 6.1 The p-value does not test direction

With a single instrument, MR-GGI's test statistic algebraically reduces to

```
    Bxy / se_Bxy  =  Bzy / se_Bzy  *  sign(Bzx)
```

The exposure's first stage **cancels out entirely**. The p-value therefore tests
*instrument → outcome*, and is significant whenever `V1` is associated with the outcome
gene — even when the exposure has no instrument at all and the reported ratio is a
division by approximately zero.

Measured on a known-null trio (`model0`, so `V1` is unrelated to `T2`), asking for
`T2 -> T1`:

| quantity | value |
| --- | --- |
| reported `Bg2g1` | **−38.42** |
| reported p-value | **0.000e+00** |
| plain p for `V1 -> T1` (the *outcome*) | 1.363e−37 — what the p-value is tracking |
| first stage `V1 -> T2` (the *exposure*) | p = 0.728, **F = 0.1** — no instrument |

**Table 6. A confidently significant answer with no instrument behind it.** This is why
edge calls are gated on the exposure's first-stage F rather than on the reported p-value.

The gate is `mrggi.min.F = 10`, the conventional weak-instrument threshold (Staiger &
Stock). Trios failing it are recorded with `mrggi.weak.instrument = TRUE` and no edge.

### 6.2 Only the cis → trans direction is usable

`V1` is the cis gene's eQTL by construction, so `T1` is the only gene in a trio with a
legitimate instrument. The `T2 -> T1` direction is computed and stored but never used for
edge calls.

| model | `T1 -> T2`: B, p, F(T1) | `T2 -> T1`: B, p, F(T2) |
| --- | --- | --- |
| model0 (no edge) | −0.026, p 0.728, F = 178.6 | −38.4, p < 1e−4, **F = 0.1** |
| model1 (`T1 -> T2` = 0.7) | **0.674**, p < 1e−4, F = 178.6 | 1.484, p < 1e−4, F = 43.0 |
| model2 (`T2 -> T1` = 0.7) | −0.028, p 0.728, F = 82.8 | −36.0, p < 1e−4, **F = 0.1** |
| model3 (no edge, `V1 -> T2`) | 0.991, p < 1e−4, F = 178.6 | 1.010, p < 1e−4, F = 175.3 |

**Table 7. The forward direction behaves; the reverse does not.** Synthetic trios,
n = 1000, true effect 0.700 where an edge exists. The forward direction recovers 0.674 for
model1 and is correctly null for model0 and model2. The reverse is uninterpretable
throughout: F = 0.1 exposes model0 and model2 as instrument-free, while model1 and model3
have *strong* reverse instruments that are nonetheless invalid, because `V1` reaches the
outcome other than through the exposure.

**A consequence worth stating plainly: MR-GGI cannot detect model2's edge on these
trios.** Recovering `T2 -> T1` needs an instrument for `T2`, and the design does not
provide one.

The obvious response is to swap the two genes — treat `T2` as the exposure and `T1` as the
outcome, which is structurally model1 mirrored — and it does not work, for a reason worth
recording. MR requires the *exposure* to be anchored by the instrument, and in model2
`V1 -> T1 <- T2` makes `T1` a collider, so `V1` and `T2` are marginally independent.
Measured over the 300 real trios of the n = 1000 group:

| model | \|r(V1,T2)\| | F for `T2` as exposure | trios with F > 10 |
| --- | --- | --- | --- |
| model0 | 0.028 | 0.8 | 0% |
| **model2** | **0.019** | **0.3** | **0%** |
| model1 | 0.120 | 14.6 | 52% |
| model3 | 0.278 | 83.4 | 80% |
| model4 | 0.337 | 128.0 | 85% |

**Table 7b. Whether `V1` can instrument the trans gene, by generating model.** `T2` is
instrumentable only where `V1 -> T2` actually exists, i.e. model3 and model4 — and there
the instrument is *invalid* for a `T2 -> T1` test anyway, since `V1` also reaches `T1`
directly. The swap therefore fails two different ways, and the first-stage F separates
them: **model0 and model2 fail loudly** (F below 1, no instrument at all), while
**model1, model3 and model4 fail silently** — a strong-looking first stage attached to a
violated exclusion restriction, which is how model1's reverse direction returns a
confident 1.484 against a truth of no `T2 -> T1` edge.

Note this makes MR-GGI *incomplete* on model2 rather than wrong: it correctly reports no
`T1 -> T2` edge (Table 8: 0% called), which is true. model2 simply lands in the same
bucket as model0. Resolving it would require the trans gene to carry its own cis-variant,
which is precisely what MR-GGI assumes (Table 3) and what a trio by construction lacks.

### 6.3 Horizontal pleiotropy in model3 and model4

MR requires the instrument to affect the outcome *only* through the exposure — the
exclusion restriction. In model3 and model4, `V1 -> T2` directly, which violates it. For
**model3** in particular there is no true T1–T2 edge, yet `beta(V->T2)/beta(V->T1)` is
non-zero, so a spurious edge is expected.

Measured on 40 real trios per model from the simulation:

| model | true T1–T2 edge | n = 50: called | n = 1000: called | n = 1000: mean B |
| --- | --- | --- | --- | --- |
| model0 | none | 0.0% | 12.5% | −0.071 |
| model1 | **yes** | 0.0% | **62.5%** | 0.401 |
| model2 | `T2 -> T1` | 0.0% | 0.0% | 0.088 |
| model3 | **none** | 25.0% | **87.5%** | 0.965 |
| model4 | **yes** | 12.5% | 75.0% | 1.201 |

**Table 8. Measured edge-call rates by generating model.** `called` is the proportion of
trios where the `T1 -> T2` edge was called at `p < 0.05` with first-stage `F > 10`. At
n = 1000 the method behaves as designed on model0, model1 and model2, but calls a
**spurious edge in 87.5% of model3 trios** — the pleiotropy prediction, confirmed. The
inflated model4 estimate (1.201 against a true 0.7 mediation effect) is the same
violation biasing a real edge rather than manufacturing one.

This is a known limitation of Mendelian randomisation, not an implementation fault, and it
is a genuine finding about the method's applicability: MR-GGI is unreliable exactly where
the variant has a direct effect on the trans gene.

### 6.4 Weak instruments at small n

At n = 50 the F gate rejects most trios — `mrggi.weak.instrument` is TRUE for 62.5% to
100% of them by model, against 0% to 37.5% at n = 1000. That is the gate doing its job:
the ungated estimates at n = 50 are wild (mean B of 1.96 for model1, 4.74 for model3),
because the Wald ratio divides by a first stage that is near zero. Median `|cor(V1,T1)|`
is 0.16 at n = 50 (see `METHODS.md` §3), which is not enough to instrument with.

Practically, MR-GGI has **almost no usable output below n = 300**, and this should be
reported alongside its accuracy rather than folded into it.

---

## 7. Settings and provenance

| setting | value | where |
| --- | --- | --- |
| `mrggi.alpha` | 0.05 | `inference_config.R` |
| `mrggi.cor.thr` | 0 | `inference_config.R` |
| `mrggi.min.F` | 10 | `inference_config.R` |
| instruments | exposure `V1`, outcome a column of zeros | `mrggi.one.trio()` |
| confounders | none | §4 |

`cor.thr = 0` is the package default. It screens gene pairs by correlation before testing,
and with a single pair per trio it only decides whether that trio is tested at all. The
reference script's `0.3` would silently discard 36% of model1, 41% of model2 and 32% of
model4 trios untested (`METHODS.md` §2 correlations), capping recall before any inference
happens — so it is left at 0 and the correlation is stored as `mrggi.GGcor` instead.

**A package bug worth knowing.** `MRggi()` calls
`p.adjust(pval.idx, method = p.adjust.methods)` — `p.adjust.methods` is R's built-in
constant vector of method names, not the function's own `p.adjust.method` argument. The
user's choice is silently ignored and `match.arg` takes the first element, `"holm"`. With
one gene pair per trio this adjusts nothing, so the `FDR_*` columns equal the raw
p-values; the pipeline stores and thresholds `pval_*` directly and does not rely on
`FDR_*`.

### Columns written per trio

`mrggi.B.T1T2`, `mrggi.p.T1T2`, `mrggi.F.T1`, `mrggi.B.T2T1`, `mrggi.p.T2T1`,
`mrggi.F.T2`, `mrggi.GGcor`, `mrggi.edge` (`T1->T2` or `none`),
`mrggi.weak.instrument`, `mrggi.time.seconds`, `mrggi.error`.

### Reproducing

```
Rscript bioinfo_revision/simulation_results/mrggi_feasibility.R   # Tables 2, 4, 6, 7
Rscript bioinfo_revision/simulation_results/apply_mrggi.R         # the run itself
```

`MRggi` is installed from GitHub (`remotes::install_github("hiows/MRggi")`) and depends on
`susieR`. `MRGNgeneral` and `graph`, which appear in the reference script, are **not**
needed — they belong to a different analysis and its `RecallPrecision()` scoring, which
this pipeline does not use.

Tables 1, 2, 4, 6 and 7 come from synthetic trios with a known true effect, generated
inside `mrggi_feasibility.R`. Table 8 comes from 40 real trios per model per group drawn
from `simulation/simulated_data/simulated_trios.RData`. Last reproduced 2026-08-21.
