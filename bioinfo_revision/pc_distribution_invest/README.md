# Distribution of confounder (PC) effects on cis and trans genes

Supporting analysis for the revision: how large are the effects of the selected
confounding PCs on cis and trans gene expression in real data? The answer sets
realistic bounds for the confounder effect sizes used in the simulation study.

## Script

`compute_pc_dist_bounds.R` — run from the **repository root**. It records the
root as the working directory, loads `GTEx/data/data.with.PCs.WholeBlood.RData`
from there, switches into this folder to write the two PNGs, and restores the
root working directory on the way out.

For each of the 3,248 Whole Blood trios it regresses the cis gene (and then the
trans gene) on the PCs selected for that trio, and pools the slopes across all
trios into one distribution.

Each element of `data.with.PCs.WholeBlood.RData` is a `670 x p` data frame laid
out by `GTEx/data/PC_LRNA_PC_Selection_manu.R`:

| columns | contents |
| --- | --- |
| 1 | SNP genotype (V) |
| 2 | cis gene expression (T1) |
| 3 | trans gene expression (T2) |
| 4-6 | known covariates: `pcr`, `platform`, `sex` |
| 7+ | the PCs selected for that trio (`PC<k>`) |

Three details matter for the result to be interpretable:

- **Only `PC*` columns are used.** `pcr`, `platform` and `sex` are binary known
  covariates whose slopes are an order of magnitude larger than the PC slopes,
  so including them in the same histogram mixes two different scales.
- **Degenerate PCs are dropped** (`SD.TOL = 1e-08`). `PC670` is the null
  direction left over after centering 670 samples; the PC-selection step retains
  it in 431 trios with SD ~1e-13. It is not exactly constant, so `lm()` does not
  alias it away — its slope instead explodes to ~1e11 and swamps every real
  effect.
- **Effects are standardized** as `b * sd(PC) / sd(Y)`. The residualized
  expression SDs span 8 orders of magnitude across trios (4e-04 to 3.8e+04), so
  raw slopes are not comparable trio to trio. Here that standardized slope is
  *exactly* the Pearson correlation `cor(PC, Y)`: a standardized regression
  coefficient equals the marginal correlation when the predictors are mutually
  uncorrelated, and principal components are orthogonal by construction. Checked
  on the real trios — the largest off-diagonal correlation among a trio's
  retained PC columns is ~5e-15, and the standardized coefficients match
  `cor(PC, Y)` to ~2e-15. Note this is the *marginal* correlation, not the
  partial correlation (which differs by up to ~0.05).

## Result

95,564 PC-gene correlations per target:

| target | 2.5% | median | 97.5% | sd | max abs |
| --- | --- | --- | --- | --- | --- |
| cis | -0.203 | 0.000 | 0.200 | 0.117 | 0.661 |
| trans | -0.189 | 0.002 | 0.191 | 0.104 | 0.653 |

![cis](PC_effects_distribution_cis.png)
![trans](PC_effects_distribution_trans.png)

Both distributions are symmetric about zero and effectively bounded by
`|effect| <~ 0.3`, with a hard ceiling near 0.66. cis and trans genes behave
almost identically.

The dip at zero is not an artifact of the fit: these PCs were *selected* for
significant association with the trio (FDR 0.05, see `get.conf.trios` in
`PC_LRNA_PC_Selection_manu.R`), so near-null effects are filtered out by
construction. This is the distribution of effects **conditional on selection**,
which is the relevant one for choosing simulation parameters. Regressing on the
full `PCs.matrix` instead of the selected subset would give the unconditional
distribution.
