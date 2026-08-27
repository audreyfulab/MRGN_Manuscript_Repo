# `gtex_bootstrapping/` — MRGN with edge-probability bootstrap on real GTEx trios

MRGN applied to the 3,248 GTEx Whole Blood trios, each with the bootstrap that estimates how
often every edge survives resampling. Written by
[`../../gtex_bootstrapping/apply_mrgn_gtex.R`](../../gtex_bootstrapping/apply_mrgn_gtex.R).

This is the real-data counterpart to the simulation stage. The difference that shapes
everything here: **there is no ground truth**. Nothing can be scored `correct`, so the
output is the model call, the bootstrap support behind it, and enough provenance to join
back to the trio. The bootstrap is the only available statement about confidence.

## Running

```
Rscript bioinfo_revision/gtex_bootstrapping/apply_mrgn_gtex.R
```

| flag | default | what |
| --- | --- | --- |
| `--cores N` | `detectCores() - 2` | cluster size |
| `--bootstrap N` | 1000 | resamples per trio |
| `--max-trios N` | all | smoke tests only |
| `--chunk-size N` | 200 | trios per checkpoint |
| `--rerun 0\|1` | 0 | redo chunks that already exist |
| `--gtex-file PATH` | `GTEx/data/data.with.PCs.WholeBlood.RData` | input |
| `--out-dir PATH` | this directory | output |
| `--tissue NAME` | `WholeBlood` | recorded in the `tissue` column |

Chunks are checkpointed to `chunks/chunk_NNNN.RData` and skipped if present, so an
interrupted run resumes rather than restarting. Delete a chunk, or pass `--rerun 1`, to
recompute.

## Input

`GTEx/data/data.with.PCs.WholeBlood.RData` holds `data.sets`, an unnamed list of 3,248
data.frames of 670 samples each, already carrying their selected confounders. Each is laid
out exactly as `MRGN::infer.trio()` expects:

| columns | contents |
| --- | --- |
| 1 | the variant, `chr..._b38` |
| 2 | the *cis* gene `T1`, `ENSG...` |
| 3 | the *trans* gene `T2`, `ENSG...` |
| 4+ | clinical knowns (`pcr`, `platform`, `sex`) then the selected PCs |

Between 12 and 57 columns, median 36 — so 9 to 54 covariates per trio. `infer.trio()` takes
no confounder-count argument: it reads the first three columns as the trio and everything
after as covariates, so the frames go in unmodified.

## Output

| file | what |
| --- | --- |
| `gtex_mrgn_bootstrap.RData` | `gtex.mrgn.results`, one row per trio |
| `gtex_mrgn_bootstrap.csv` | the same, tidy |
| `chunks/chunk_NNNN.RData` | per-chunk checkpoints; safe to delete once combined |
| `progress.log` | one line per trio, written live by the workers |

### Watching a run

`progress.log` gets a line as each trio finishes, from whichever worker ran it:

```
trio    21 complete! inferred as M1.1   | boot M1.1   (agrees) | 27.6s
trio    18 complete! inferred as M4     | boot M4     (agrees) | 28.7s
```

The same lines appear in the console. That works **only** because the cluster is built with
`makeCluster(n.cores, outfile = "")` — a PSOCK cluster without it routes worker stdout to a
null device, so every one of these would be silently discarded while the code still looked
correct. Checkpointing is per chunk, so without this the run appears frozen for 200 trios at
a time.

The log is truncated at the start of each run: it describes the run in progress, while the
chunk files are the durable record. Twelve workers append to it concurrently; each write is
one short line in append mode, and a garbled line would cost a log entry and nothing else,
since results return through the cluster rather than through this file.

### Columns

| column | meaning |
| --- | --- |
| `trio.index` | position in `data.sets`, the join key back to the input |
| `tissue`, `variant`, `cis.gene`, `trans.gene` | identity |
| `n.samples`, `n.covariates`, `n.known`, `n.pcs` | shape of what was fitted |
| `model` | the inferred model, `infer.trio()$Inferred.Model` |
| `time.seconds` | the point inference |
| `boot.model` | the model a majority vote of the resamples supports |
| `boot.agrees` | whether `boot.model` equals `model` |
| `boot.min.edge.prob` | weakest edge in the supported model — the trio's confidence floor |
| `boot.p.V1T1`, `boot.p.T1T2`, `boot.p.V1T2`, `boot.p.T2T1` | per-edge support |
| `boot.n.requested`, `boot.n.used`, `boot.n.dropped` | resample bookkeeping |
| `bootstrap.time.seconds`, `bootstrap.error` | bootstrap cost and failure |
| `error` | the point inference failed; `model` is `NA` |

`boot.n.dropped` counts resamples thrown away for having no genotype variation — a draw that
loses every copy of the minor allele leaves `V1` constant, which says nothing about the `V1`
edges and would otherwise either abort the bootstrap or drag every indicator toward zero.
At n = 670 this is rare; it mattered at n = 50 in the simulation stage.

A failed point inference and a failed bootstrap are recorded separately. The bootstrap is
wrapped in its own handler inside `apply.mrgn()`, so a bad resample costs the `boot.*`
columns and not the model call.

## The parallelism, and why it is arranged this way

`boostrap_edge_probabilities()` accepts a cluster and will spread its resamples across it.
**This script deliberately does not use it that way.** Each resample is one `infer.trio()`
call on 670 rows — a few milliseconds — so dispatching 1,000 of them through `parLapply`
costs more in socket round-trips than it saves in arithmetic.

That is measured, not assumed. On the confounder-structure runs of 2026-08-26, three
concurrent 4-core clusters bootstrapping per-resample held **1.9 core-equivalents of 14**
busy, with the disk at 0.3% and 8 GB of RAM free. The work was not compute-bound; it was
waiting on the cluster.

So the cluster parallelises the **outer** loop: one worker takes one whole trio and runs its
1,000 resamples serially in-process. One dispatch per trio instead of a thousand. Measured
at 6.9–9.6 s per trio serially, i.e. ~7.2 h single-threaded against roughly 35–60 min across
a full cluster.

`clusterMap` is used rather than `parLapply` over indices for a related reason: a closure
referring to `datasets[[i]]` captures the whole 755 MB list and serialises it to every worker
on every chunk.

### Parallel efficiency is real but far from linear

Measured on the full run, 12 workers, 1,000 replicates:

| | per trio |
| --- | --- |
| serial, one trio at a time | 7–10 s |
| under 12-way concurrency | 28–32 s |

So a trio takes roughly **3× longer** when eleven others run beside it, and the throughput
gain is nearer 4× than 12×. The observed rate is ~24 trios/min, i.e. ~8.6 min per 200-trio
chunk and about 2 h for all 3,248 — against ~7.2 h serial.

The cause is almost certainly memory bandwidth: every worker holds its own 670 × N frame and
drives many regressions over it, and the machine's cores are not the binding resource. Do
not read a core count as a speed-up factor when planning a re-run.
