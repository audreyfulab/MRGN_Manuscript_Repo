# `pc_distribution_invest/data/`

Cached output of `../compute_pc_dist_bounds.R`. Not tracked in git — regenerate by
rerunning that script from the repository root.

## `real_pc_effect_pools.RData` (1.4 MB)

Loads one object, `pools`, a list of five elements:

| element | type | contents |
| --- | --- | --- |
| `cis` | `numeric[95564]` | per-PC standardized effect on the cis gene, pooled over all trios |
| `trans` | `numeric[95564]` | same for the trans gene |
| `n.trios` | `integer` | 3248 |
| `source` | `character` | `"GTEx/data/data.with.PCs.WholeBlood.RData"` |
| `scale` | `character` | `"Pearson correlation cor(PC, gene), = b * sd(PC) / sd(Y)"` |

Load it with `MRGN::loadRData()`, which returns the object rather than assigning it:

```r
pools <- MRGN::loadRData("bioinfo_revision/pc_distribution_invest/data/real_pc_effect_pools.RData")
sd(pools$cis)   # 0.117
sd(pools$trans) # 0.104
```

## What it is for

`simulation/verify_simulation.R` reads this file and overlays `pools$cis` /
`pools$trans` against the same quantity measured on the simulated n = 670 trios,
producing `simulation_results/simulated_vs_real_conf_effects.png`. It is the reference
distribution for the `conf.coef.ranges$U` bound.

Two properties are load-bearing:

- **Only the standardized effects are stored.** The correlation scale is unit-free, so
  0.117 means the same thing in a GTEx trio and in a simulated one. The raw slopes are
  not portable — `sd(cis gene)` spans 0.042 to 3.75e+04 across trios, so a raw slope
  carries the units of one particular trio's residualized expression. The writer errors
  out if any stored value falls outside `[-1, 1]`, which would mean a raw slope leaked in.
- **The pools are conditional on selection.** These PCs were retained at FDR 0.05 by
  `get.conf.trios()` during PC selection, which is why the densities dip at zero. That
  is the right reference for choosing simulation parameters, but it means the simulated
  `U` block represents *confounders strong enough to have been found*, not all
  confounders. See [`../../METHODS.md`](../../METHODS.md) §6.
