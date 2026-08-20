# `bioinfo_revision/` — revised simulation study

Everything added for the manuscript revision: the real-data measurements that
set the simulation's parameter bounds, the simulation itself, the calibration check,
and the inference run over the simulated trios.

[`METHODS.md`](METHODS.md) is the methodology record — what each parameter is, why it
has the value it has, and what the run actually realized. This file is the map; the
per-folder READMEs describe the files.

## Pipeline

Run every script from the **repository root** (`MRGN_Manuscript_Repo/`), not from the
folder it lives in. Each script either records the root with `getwd()` and `setwd()`s
back on the way out, or uses root-relative paths throughout.

| # | stage | script | writes |
| --- | --- | --- | --- |
| 1 | real confounder (PC) effects | [`pc_distribution_invest/compute_pc_dist_bounds.R`](pc_distribution_invest/) | `pc_distribution_invest/data/real_pc_effect_pools.RData`, 2 PNGs |
| 2 | real SNP effects | [`pc_distribution_invest/compute_effects_snp_on_gene.R`](pc_distribution_invest/) | 2 PNGs |
| 3 | simulate 3,750 trios | [`simulation/updated_data_simulation.R`](simulation/) | `simulation/simulated_data/simulated_trios.RData` |
| 4 | calibration check | [`simulation/verify_simulation.R`](simulation/) | console report + `simulation_results/simulated_vs_real_conf_effects.png` |
| 5 | inference (MRGN / MRPC / GMAC) | [`simulation_results/updated_simulation_inference.R`](simulation_results/) | `simulation_results/inference_results.RData` + `.csv` |

Stages 1–2 read `GTEx/data/data.with.PCs.WholeBlood.RData`; stage 4 reads the output
of stages 1 and 3; stage 5 reads the output of stage 3. Stages 1 and 2 are
independent of each other, as are 4 and 5.

## Folders

| folder | contents |
| --- | --- |
| [`pc_distribution_invest/`](pc_distribution_invest/) | Real GTEx Whole Blood effect sizes: PC-on-gene and SNP-on-gene distributions, the figures, and the cached effect pools the calibration check compares against |
| [`simulation/`](simulation/) | Data generation: the scenario grid, the generating helpers, the calibration check, and the generated trios |
| [`simulation_results/`](simulation_results/) | Inference over the simulated trios and the calibration figure |

## Requirements

R 4.4.x with `MRGN`, `MRPC`, `mvtnorm`, `ggplot2`, `patchwork`, `gridExtra`,
`R.utils`, `parallel`, and (via `MRGN`) `propagate` and `qvalue`. The inference stage
also sources the adapted GMAC code in [`../adapted_GMAC_func/`](../adapted_GMAC_func/).

## Data files are not in git

The two `data/` and `simulated_data/` folders hold generated `.RData` (up to 368 MB
each) and are untracked. Regenerate them by running stages 1 and 3. `.gitignore` also
excludes the GTEx input `GTEx/data/data.with.PCs.WholeBlood.RData`, so stages 1 and 2
cannot be rerun without obtaining that file separately.
