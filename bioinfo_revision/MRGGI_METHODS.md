# MR-GGI on trios — methodology

How MR-GGI (https://github.com/hiows/MRggi) is applied to the simulated trios, why three
adaptations are needed, what the four covariate arms do and do not change, and what the
method can and cannot answer. Companion to [`METHODS.md`](METHODS.md), which covers the
simulation and the other three methods.

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

### 3.1 What "no instrument" has to look like

The "no instrument" entry has to be a **column of zeros**:

| outcome-gene entry | result |
| --- | --- |
| `NULL` | error: `'data' must be of a vector type, was 'NULL'` |
| zero-**column** matrix (`ncol == 0`) | error: `subscript out of bounds` |
| **column of zeros** (`ncol == 1`, all 0) | **works** — `Bg1g2 = 0.674`, reverse direction `NaN` |
| `V1` (same instrument as the exposure) | `Bg1g2 = 0.0000, p = 1.0000` — the collapse of §2 |

**Table 4. What can stand in for "this gene has no instrument".** The `NaN` in the
reverse direction is correct behaviour, not a failure: with no instrument for `T2`, the
`T2 -> T1` effect is not estimable.

**Why a column of zeros survives is not the obvious reason.** `MRggi()` opens with
`scale.X = lapply(X, scale)`, and `scale()` on a constant column returns all `NaN` —
centring leaves zeros, then it divides by `sd = 0`. That would be fatal if the value were
used. It is not: **`scale.X` is assigned and never read again**, and the main loop takes
the raw `X[[i]]`. Only `NULL` fails at that line, because `lapply` still has to evaluate
it. The zero-column matrix gets past `scale()` too and dies later, inside `.TSLS`.

### 3.2 Two further requirements, both easy to trip over

**`X` must be positionally aligned with `y`.** `X[[i]]` is the instrument set for column
`i` of `y`, so `length(X)` must equal `ncol(y)`; a shorter list gives
`subscript out of bounds`.

**`colnames(y)` must be set.** `MRggi()` builds its output with
`g1 = append(g1, colnames(y)[i])`. With `NULL` colnames `g1` stays empty and the function
dies at its closing `data.frame()` with `arguments imply differing number of rows: 0, 1`
— an error naming nothing involved. The trap is that `cbind(T1 = a, T2 = b)` on unnamed
`n × 1` **matrices** yields `NULL` colnames, because `cbind` ignores the tag for matrix
arguments. `mrggi.one.trio()` sets the names explicitly and asserts them.

---

## 4. Covariates in `y` — what they do and do not change

Each trio is passed to `MRggi()` in **one call**, as `y = (T1, T2, covariates)`. Four
covariate sets are run per trio: `none` (the bare trio), `truth`, `CSq` and `CSa`.

**This is not a confounder adjustment, and must not be read as one.** `MRggi()` takes no
covariate argument, and the estimation is strictly pairwise: the only data entering the
regressions for the T1–T2 row are `y[,T1]`, `y[,T2]`, `X[["T1"]]` and `X[["T2"]]`. Adding
covariates to `y` therefore **changes the T1–T2 estimate not at all**:

| `y` | `Bg1g2` | `pval_Bg1g2` | `FDR_Bg1g2` | rows |
| --- | --- | --- | --- | --- |
| `(T1, T2)` | 0.808 | 0.000 | 0.435 | 1 |
| `(T1, T2, U1…U18)` | **0.808** | **0.000** | 1.000 | 190 |

**Table 4b. The estimate is invariant; the multiplicity correction is not.** Same trio,
measured. What the extra columns change is `FDR_Bg1g2`: `MRggi()` adjusts each `g1`'s
p-values across that gene's pairs, so the T1–T2 p-value is now corrected for T1's pairs
against every covariate too. A wider covariate set is a harsher correction, and that is
the **entire** difference between the four arms.

The results therefore carry two edge calls. `edge` comes from the raw p-value and is
identical in every arm — it is the column comparable with GMAC and MRGN, neither of which
is corrected. `edge.fdr` comes from `FDR_Bg1g2` and is the only column on which the arms
differ. `confusion_mrggi.R` **asserts** the arm-invariance of `edge`; if it ever fails,
`X` has stopped lining up with the columns of `y`.

> **Package bug worth knowing.** `MRggi()`'s `p.adjust.method` argument is ignored. The
> body calls `p.adjust(pval.idx, method = p.adjust.methods)` — note the trailing `s`, base
> R's vector of *all* method names — so `match.arg` silently takes the first and the
> correction is **always holm**. Verified: `p.adjust(c(.01,.02,.03,.04), method =
> p.adjust.methods)` returns `0.04 0.06 0.06 0.06`, which is holm, not bonferroni. Read
> `mrggi.<arm>.FDR.T1T2` as holm-adjusted regardless of `mrggi.p.adjust`.

This differs from MRPC, where confounders as graph nodes genuinely do adjust, because the
PC algorithm's conditional independence tests condition on subsets of the other nodes.
The same phrase — "confounders in the model" — means different things for the two methods.

So the only way to make MR-GGI's *estimate* respond to a confounder set would be
residualising `T1` and `T2` on it before inference. **That is not done**, on the reasoning
that the instrument is what is supposed to handle confounding, and our `U` block is
precisely the mediator–outcome confounding MR exists to defeat. Table 1 supports this:
unadjusted, the Wald ratio's bias is −0.008 against the naive regression's +0.142.
Residualising buys an 18% reduction in standard error — a precision gain, not a bias
correction — and would hand MR-GGI an advantage the other methods do not get for free.

### 4.1 Cost

`cor.thr = 0` means every gene pair is computed, so an arm with `k` covariates costs
O((k+2)²) TSLS fits. Measured, minutes per 300-trio group:

| arm | n=50 | n=150 | n=300 | n=670 | n=1000 |
| --- | --- | --- | --- | --- | --- |
| `none` | ~0 | ~0 | ~0 | ~0 | ~0 |
| `truth` | 6.2 | 7.2 | 9.8 | 10.7 | 9.6 |
| `CSq` | ~0 | ~0 | 0.4 | 3.6 | 8.4 |
| `CSa` | 67.2 | 94.0 | 111.0 | 134.3 | 159.1 |

**Table 4c.** CS-alpha is 9.4 h of the 10.4 h total, because it selects a median of 82–106
covariates — roughly 5,800 pairs per trio, nearly all returning `NaN` since a covariate has
no instrument and cannot be an exposure. This is why `apply_mrggi.R` now builds a cluster
rather than running single threaded.

### 4.2 The correlation screen, `cor.thr`

`MRggi()`'s second argument after the data is a correlation threshold, and it is easy to
mistake for confounder selection. It is not. It gates **which gene pairs are estimated at
all**, on their marginal correlation:

```r
cor.y    <- cor(scale.y)
corMat[which(abs(cor.y) > cor.thr)] <- 1
corMat[lower.tri(corMat, diag = TRUE)] <- 0
calc.idx <- which(corMat == 1, arr.ind = TRUE)
```

It never reaches `.TSLS()`, so for a pair that survives it, `B` and `p` are unchanged. What
it changes is which pairs exist to be corrected over — and, for the trio, whether the trio is
analysed at all.

**Set to `0.1`** (`mrggi.cor.thr`). It was `0` for the earlier runs, which was the package
default and harmless while MR-GGI saw only the bare trio: with one pair, a screen has nothing
to choose between. Once the covariate arms were added, `0` meant every pair in the upper
triangle was computed — 4,950 per trio under CS-alpha, of which we read one.

**Two consequences, and the second is the one to carry into the results.**

**(a) It sets the multiplicity family.** `FDR.T1T2` is `p.adjust()` over the pairs sharing
`g1 == "T1"`, so a higher threshold drops weakly-correlated partners out of that family and
corrects the T1–T2 p-value less harshly. Only `edge.fdr` moves; `edge`, from the raw p, does
not.

**(b) It decides whether the trio is tested at all.** If `|cor(T1, T2)| <= cor.thr` the T1–T2
pair never enters `calc.idx`, there is no estimate, and the trio is a **no-call**. Measured
over all 1,500 trios at `cor.thr = 0.1`:

| generating model | true T1–T2 edge | screened out |
| --- | --- | --- |
| model0 | **none** | **65.7%** |
| model3 | **none** | **43.0%** |
| model1 | yes | 8.3% |
| model2 | yes | 5.0% |
| model4 | yes | 6.3% |

**Table 4d. The screen is not neutral across the generating models.** ~26% of trios overall,
concentrated almost entirely in the two models that have **no** T1–T2 edge — they are the
trios whose T1–T2 correlation is near zero, which is precisely what the screen exists to
remove. By sample size the overall rate is stable (21.7 / 24.0 / 29.3 / 27.7 / 25.7% at
n = 50 → 1000), so this is a property of the models, not of power.

**Read every MR-GGI edge number with that in mind.** MR-GGI is excused from answering on most
of its true negatives, which **raises its edge precision and lowers its edge-absent recall by
construction**. Neither is comparable with a method that answered on every trio unless the
screened-out share is quoted beside it.

This is the package working as designed — MR on an uncorrelated pair is not meaningful — so
the screen is kept and its cost reported, rather than disabled to make the columns look
comparable. Screened trios are recorded as `screened.out = TRUE`, written as
`edge = "screened"`, and tabulated in their own `Screened out` row, distinct from
`Weak instrument`: one is a trio MR-GGI tried to test and could not, the other a trio it
never looked at. Folding either into `Edge Absent` would count a refusal to answer as a
correct rejection.

**One implementation note.** The screen is checked in `mrggi.one.trio()` *before* calling
`MRggi()`, not only by looking for a missing T1–T2 row afterwards. When nothing survives the
screen — which is every screened trio in the `none` arm, having only the one pair — `calc.idx`
is empty and `MRggi()`'s main loop runs `for (k in 1:nrow(calc.idx))` with `nrow = 0`. R
evaluates `1:0` as `c(1, 0)`, the body executes on a zero-row frame, and the function dies.
Without the pre-check a screened trio is indistinguishable from a genuine failure; that is
exactly how the first run produced 8 all-`NA` rows out of 20.

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
(`RecallT1`, `MRGN_v8.pdf` line 377). Cross-tabbing against the five models belongs in the
analysis stage.

A model label **is** produced as well, from a second MRggi call plus the four tests
`class.vec()` needs that a single call does not supply. It is reported because MR-GGI's
performance on it is a result about the method, not because it works — §5.2 measures it,
and it does not.

### 5.1 How the model call is built

`MRGN::class.vec()` takes six binary indicators and returns a model label; `MRGN::get.adj()`
turns the same vector into the 3 × 3 adjacency. Every one is obtainable here:

| indicator | edge | source |
| --- | --- | --- |
| `b11` | `V1 → T1` | `lm(T1 ~ V1 + T2)`, coefficient on `V1` |
| `b12` | `T1 → T2` | **MRggi**, `y = (T1, T2)`, `X = (V1, zeros)` |
| `b21` | `V1 → T2` | `lm(T2 ~ V1 + T1)`, coefficient on `V1` |
| `b22` | `T2 → T1` | **MRggi**, `y = (T2, T1)`, `X = (V1, zeros)` — the *swapped* call |
| `V1:T1` | marginal | `cor.test(V1, T1)` |
| `V1:T2` | marginal | `cor.test(V1, T2)` |

Each is significant at `mrggi.alpha`. A test that cannot run — a constant column, an aliased
coefficient — scores 0 rather than propagating `NA`, so `class.vec()` always receives a
complete vector; the count is kept in `mrggi.n.tests.failed`.

**The reverse direction needs a second call with the genes swapped, not `V1` added to
`X[[2]]`.** Table 4 records what the latter does: `Bg1g2 = 0.0000, p = 1.0000`, and the same
in reverse. Swapping the columns instead keeps the outcome uninstrumented — the §3 setup —
with `T2` now the instrumented exposure, and returns real estimates.

**`b21` is conditional on `T1`, and has to be.** The marginal alternative `lm(T2 ~ V1)` is
exactly the `V1:T2` test already in the vector, and M2 and M4 differ *only* in that marginal
(`MRGN_v8.pdf` Table 1). A marginal `b21` would collapse them together and cap this at four
reachable labels.

**The label is arm-invariant**, so it is computed once per trio and written as
`mrggi.model`, not `mrggi.<arm>.model`. `b12` and `b22` are pairwise TSLS estimates that
covariates in `y` cannot move — the property behind the arm-invariant raw-p `edge` call of
§4 — and the four regression tests condition on the trio alone.

### 5.2 Measured: the causal tests are not independent of the marginals

The classifier is faithful to MR-GGI and still fails, for a reason worth stating precisely.

With a **single** instrument the Wald-ratio p-value reduces to the instrument→**outcome**
t-statistic — the exposure's first stage cancels out of the test entirely (`inference_config.R`,
`mrggi.min.F`). So the forward call's p-value *is* the `V1 → T2` test, and the swapped
call's p-value *is* the `V1 → T1` test. Measured over 100 trios at n = 1000, large effect:

| identity | per-trio agreement |
| --- | --- |
| `b12` = `V1:T2` marginal | **100%** |
| `b22` = `V1:T1` marginal | **100%** |

`b22`'s p-values have zero variance across those trios: `V1 → T1` is strong in all five
models, so the reverse test is maximally significant everywhere. The six-vector therefore
carries four distinct tests, not six, with the two MR estimates duplicating the two
marginals.

The consequence, on the same trios — first-stage `F` median 232, **no** weak instruments, so
this is MR-GGI's best case:

| true | M1.1 | M1.2 | M4 | Other |
| --- | --- | --- | --- | --- |
| model0 | 0 | 0 | 1 | **19** |
| model1 | 14 | 0 | 6 | 0 |
| model2 | 0 | 0 | 0 | **20** |
| model3 | 0 | 0 | **20** | 0 |
| model4 | 0 | 7 | 13 | 0 |

**Table 6. MR-GGI's model call at n = 1000, large effect.** M0 and M2 are never recovered:
a spuriously significant `b22` on every trio breaks the pattern `class.vec()` needs for
either, and both fall through to `Other`. model3 goes to **M4 in 20 of 20** — the §6.3
pleiotropy result arriving in the model call, where a real `V1 → T2` path plus a spurious
`T1 → T2` edge is exactly M4's signature.

**This is a property of single-instrument MR, not of the implementation.** One instrument
identifies one edge. The swap fixed the numerical collapse of Table 4 but could not create a
gene–gene test independent of the instrument–gene associations, because with one instrument
there is none to create. The columns are kept because "MR-GGI cannot classify trios" is a
finding, and one better supported by a scored table than by an assertion.

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
legitimate instrument.

Under the one-call design of §4 the reverse direction is **no longer computed**: `T2` is
given a column of zeros, so `Bg2g1` comes back `NaN` by construction rather than being
calculated and then discarded. The table below is from the earlier two-call
implementation and is kept because it is the evidence for the decision — it shows what
was being thrown away, and why obtaining it was not worth a second `MRggi()` call per
trio per arm. `mrggi_feasibility.R` still reproduces it.

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
user's choice is silently ignored and `match.arg` takes the first element, `"holm"`.

**This used to be harmless and no longer is.** With one gene pair per trio the correction
had nothing to adjust and `FDR_*` equalled the raw p-values. Under the one-call design of
§4 an arm carries `k` covariates, so `FDR_Bg1g2` is now a genuine holm correction over
T1's `k+1` pairs — and it is the only quantity that distinguishes the arms. It is stored
as `mrggi.<arm>.FDR.T1T2` and read as **holm-adjusted regardless of `mrggi.p.adjust`**.
Any other correction has to be recomputed from `mrggi.<arm>.p.T1T2` in the analysis stage.

### Columns written per trio, per arm

Arms are `none`, `truth`, `CSq`, `CSa`; every column is prefixed `mrggi.<arm>.`

`B.T1T2`, `p.T1T2`, `FDR.T1T2`, `GGcor`, `F.T1`, `n.covars`, `n.pairs`,
`edge` (`T1->T2` or `none`, from the raw p — identical in every arm),
`edge.fdr` (same, from `FDR.T1T2` — the only column that varies by arm),
`weak.instrument`, `time.seconds`, `error`.

`B.T2T1`, `p.T2T1` and `F.T2` are **gone**: the reverse direction is not estimable under
the one-call design (§6.2) and was never used for edge calls.

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
