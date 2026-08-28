# `legacy/second_pass/` — superseded MRPC groups, 120 s timeout, no truth arm

The MRPC checkpoints for `n = 670` and `n = 1000` as they stood before the truth arm was
added. Archived **2026-08-26**, when the large MRPC groups were dropped from the run rather
than recomputed under the current settings.

Unlike [`../first_pass/`](../first_pass/README.md), these were fitted against the CURRENT
simulated trios — `simulation/simulated_data/simulated_trios.RData`, the Aug 22
re-simulation. They are superseded by settings, not by data.

## Why they are not in `data/`

They are 49 columns; the current schema is 55. The six `mrpc.truth.*` fields did not exist
when these ran, so `combine.method("mrpc")` errors on the `rbind` if it finds them —
`inference_config.R` documents this as the reason `--rerun-inference 1` exists. Moving them
here is what unblocks a full-scope combine.

## Why they are kept

They are the only surviving measurement of MRPC's CS-q arm under the **120 s** cap at the
two largest sample sizes, and the timeout rates quoted in
[`../../inference_config.R`](../../inference_config.R) and in the MRPC notes of
[`../../results_scripts/make_all_tables.R`](../../results_scripts/make_all_tables.R) were
read off these files:

| group | CS-q timed out | of |
| --- | --- | --- |
| `n = 670` | 182 (60.7%) | 300 |
| `n = 1000` | 224 (74.7%) | 300 |

The cap was raised to 180 s afterwards, so those figures cannot be reproduced from anything
in `data/`. The CS-alpha arm was not attempted in either group.

## What is here

| file | rows x cols |
| --- | --- |
| `mrpc_group_n670.RData` | 300 x 49 |
| `mrpc_group_n1000.RData` | 300 x 49 |

Nothing in the pipeline reads this folder. To use them, load them directly — do not copy
them back into `data/` unless the truth-arm columns are added first.
