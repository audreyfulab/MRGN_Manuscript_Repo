# Shared configuration for the inference stage.
#
# Sourced by apply_mrgn.R, apply_mrpc.R, apply_gmac.R and run_all_inference.R, so every
# method sees identical paths and settings. Change a value here, not in one of the apply
# scripts, or the three methods stop being comparable.
#
# Paths are relative to the repository root; every script in this stage expects to be run
# from there.
#
# Generated files are segregated from the code: everything this stage writes goes under
# out.dir, which is `data/` rather than the script folder it used to be. logs/ and tables/
# are SIBLINGS of data/, not children -- they hold the run's diagnostics and the
# human-readable report, neither of which is the data. Re-running the simulation is then a
# matter of moving data/ aside, as legacy/first_pass/ records for the previous pass.

# The command line parser is defined here rather than at the foot of the file because
# --sim-file and --out-dir have to be in hand BEFORE the paths below are built: out.dir is
# dir.create()d a few lines down, and the selection cache is keyed on it. See the flag
# documentation in the "command line overrides" section.
inference.args <- commandArgs(trailingOnly = TRUE)
inference.opt <- function(flag, default = NULL) {
    i <- match(flag, inference.args)
    if (is.na(i) || i == length(inference.args)) return(default)
    inference.args[i + 1]
}

results.root  <- "bioinfo_revision/simulation_results"

sim.data.file <- inference.opt(
    "--sim-file", "bioinfo_revision/simulation/simulated_data/simulated_trios.RData")
clinical.file <- "./GTEx/data/kclist_top5_tiss.RData"
out.dir       <- inference.opt("--out-dir", file.path(results.root, "data"))
log.dir       <- file.path(results.root, "logs")
tables.dir    <- file.path(results.root, "tables")   # used by results_scripts/

# out.dir was the script folder and so always existed. It is a generated directory now, and
# the first thing to touch it is a save() that would fail on a fresh clone.
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# ---- per-method inference settings ----
n.bootstrap     <- 1000    # bootstrap replicates per MRGN fit

# Whether MRGN runs its bootstrap at all. THIS IS THE ENTIRE COST OF THE MRGN STAGE.
#
# The bootstrap does not produce the model call. apply.mrgn() gets that from
# MRGN::infer.trio(), which never sees the bootstrap argument and defaults it to FALSE.
# What the bootstrap adds is the boot.* block -- boot.model, boot.min.edge.prob, the four
# boot.p.* edge probabilities and bootstrap.time.seconds -- and mrgn.fields() already
# returns NA for every one of them when res$bootstrap is NULL. So the SCHEMA IS THE SAME
# either way and model / correct / correct.coarse are bit-identical; only the boot.*
# columns go empty.
#
# The cost is n.bootstrap x 3 arms x n.trios replicates per group, dispatched through
# parLapply one replicate at a time. Measured on the confounder-structure runs: three
# concurrent 4-core clusters held 1.9 core-equivalents of 14 busy, with the disk at 0.3%
# and 8 GB of RAM free -- the cluster spends its time on socket round-trips for tasks that
# are individually trivial, not on arithmetic.
#
# Set FALSE for any run whose output is read on the model call -- the confusion matrices,
# and the confounder-structure comparison of METHODS.md section 7b. Set TRUE only when the
# edge probabilities themselves are wanted.
mrgn.bootstrap  <- TRUE
mrpc.timeout    <- 180     # seconds; MRPC has taken hours on trios with many confounders
gmac.nperm      <- 1000
selection.alpha <- 0.05    # significance cutoff for the GMAC mediation calls

# Which confounder sets MRPC is run against. There are exactly two: `truth`, the oracle arm,
# added so the gap between perfect confounders and selected ones can be read for MRPC the
# way METHODS.md section 5 reads it for MRGN, and `CSq`, the CS-q selection.
#
# CS-ALPHA IS PERMANENTLY EXCLUDED FOR MRPC. It is not a budget toggle and not "off by
# default". CS-alpha hands MRPC alpha x pool false-positive covariates by construction -- a
# median of 82-102 under the current simulation -- and MRPC cannot run a PC algorithm over
# that many nodes at any sample size this study uses. There is no setting under which the
# arm returns usable results, so putting "CSa" back in mrpc.arms buys a column of NAs.
#
# This is specific to MRPC. CS-alpha is still a real arm for MRGN and MR-GGI (see
# mrggi.arms), and the CS-alpha SELECTION is still scored for every trio in the CSa.n.*
# columns. Only the MRPC fit is excluded.
#
# COST, for the record. The only CS-alpha fits ever completed were in the first pass, whose
# simulation was invalidated on 2026-08-22 and archived to legacy/first_pass/. These figures
# are therefore indicative, NOT current-simulation results:
#
#   n = 50   CS-a 82 confounders,  0% timeout, median 7.8 s
#   n = 150  CS-a 95 confounders, 79% timeout, ~106 s/trio
#
# The n = 50 fits are fast only because 82 covariates on 50 observations is rank deficient
# and bails early. By n = 150 the fits are estimable, MRPC genuinely attempts a PC algorithm
# over ~98 nodes, and the timeout rate jumps to 79% -- heading to ~100% at larger n, i.e.
# ~10 h per group to produce a column of NAs. The covariate counts do carry over to the
# current simulation (82 at n = 50, 97 at n = 150, 102 at n = 300), which is why the
# conclusion stands even though the timings are from the discarded run.
#
# CS-q costs nothing by comparison, because it selects a median of 2 covariates:
#
#   n = 50   CS-q  2 confounders,  0% timeout, median 0.0 s   |  group took 48 min
#   n = 150  CS-q  2 confounders,  0% timeout, median 0.0 s   |  group took 8.8 h
#
# Every group records CS-alpha as not attempted, with the reason in mrpc.CSa.error, rather
# than silently NA. confusion_mrpc.R reads that column and skips the arm instead of
# tabulating 300 failures.
#
# THE TRUTH ARM IS THE ONE TO WATCH, not CS-alpha. It carries every one of that trio's own
# U confounders -- a median of 25-29 -- against CS-q's 16 at n = 670, so it is strictly
# harder than the arm that already times out 61% of the time there. Smoke-test it before
# committing to a full run:
#
#   Rscript bioinfo_revision/simulation_results/apply_mrpc.R --sizes 670 --max-per-group 20
#
# and read mrpc.truth.timed.out. If it is at or near 100%, drop "truth" from mrpc.arms for
# the large groups rather than spending hours writing a column of NAs -- the same call
# already made for CS-alpha above. A disabled arm still gets its columns, with the reason
# in mrpc.<arm>.error, so the schema is unchanged either way.
mrpc.arms <- c("truth", "CSq", "CSi")   # CSa is excluded, not a toggle

# Largest sample size at which the truth arm is ATTEMPTED. Above this it is recorded as not
# attempted, exactly as a disabled arm is, and the group still gets its columns. Inf runs it
# everywhere.
#
# This is a BUDGET CONTROL, not a claim that the arm is uninformative above the threshold.
#
# MEASURED, truth arm at n = 670, 180 s cap, 10 trios:
#
#   n.conf   4  ->    0.8 s   finished
#   n.conf  18  ->    4.0 s   finished
#   n.conf  20, 21, 33, 34, 37, 38, 47, 53  ->  180.0-180.7 s, ALL timed out
#
#   8/10 timeouts | median 180.1 s | 24.1 min for 10 trios | ~12.1 h per 300-trio group
#
# Two things to read off that. First, THE CAP IS ENFORCED: the worst overrun was 180.7 s
# against a 180 s cap, i.e. 1.00x, so withTimeout() is doing its job and the arm cannot run
# away. Second, the cost is bimodal in the confounder count -- everything at or above ~20
# confounders hit the wall, everything below finished in seconds -- and the truth arm's
# median is 25-29 confounders, so it lands on the wrong side of that split at n = 670.
#
# 12 h per group to write a column that is 80% NA is not worth it by default. At n <= 300
# the CS-q arm saw no timeouts at all, so the arm is cheap and informative there.
#
# Re-measure with:
#   Rscript bioinfo_revision/simulation_results/apply_mrpc.R \
#       --sizes 670 --max-per-group 12 --out-dir /tmp/mrpc_smoke
# then read mrpc.truth.timed.out and mrpc.truth.time.seconds. Raise this value (or set it
# to Inf) if the cost is acceptable, and prefer running the large groups on their own.
mrpc.truth.max.n <- 300

# mrpc.timeout was raised from 120 to 180 s. Measured on the 120 s checkpoints, per group of
# 300 trios, CS-q arm:
#
#   n = 50 / 150 / 300     0 timeouts
#   n = 670              182 timeouts (61%)  completed: median 2.6 s, 90th 58.3 s, max 119.4 s
#   n = 1000             224 timeouts (75%)  completed: median 0.9 s, 90th 48.4 s, max 114.8 s
#
# The completed fits run right up to the old cap, so how much the extra 60 s recovers cannot
# be read off these numbers -- the new timeout rate is itself a result worth reporting. Only
# n = 670 and n = 1000 need rerunning under the new value; n <= 300 never timed out.

# ---- confounder selection, kept in step with run_confounder_selection.R ----
selection_fdr    <- 0.05
filter_fdr       <- 0.1
alpha            <- 0.01
filter_int_child <- TRUE

# FALSE reuses a cached selection when it matches the request; TRUE recomputes and
# overwrites. A cache built against different simulated data is detected and recomputed
# either way, so FALSE is safe across a re-simulation.
rerun.selection <- FALSE

# FALSE skips a sample-size group whose checkpoint for this method already exists, so an
# interrupted run resumes at group granularity. TRUE redoes them.
rerun.inference <- FALSE

# ---- scope ----
sample.sizes  <- NULL    # NULL = every group
max.per.group <- NULL    # trios per group; NULL = all. Smoke tests only.

# ---- MRGN arms ----
# MRGN fits every confounder set that is listed here. A disabled arm still writes its
# columns, with the reason in mrgn.<arm>.error, so the schema does not depend on the
# setting and checkpoints written under different settings still rbind.
#
# CSi is the CS-i selection: get.conf.trios(adjust_by = "individual"), the setting
# GTEx/data/PC_LRNA_PC_Selection_manu.R:127 actually ran and the multiplicity family GMAC
# applies. It needs a CS.i block in the selection cache -- run backfill_csi.R once, which
# derives it from the already-cached reg.pvalues in ~20 s per group rather than re-running
# the ~40 min selection.
mrgn.arms <- c("truth", "CSq", "CSa", "CSi")


# ---- parallelism ----
# Cores for THIS process. run_all_inference.R overrides it per method with --cores, since
# three concurrent processes share one machine and would otherwise each size a cluster
# from detectCores() and oversubscribe it threefold.
n.cores <- max(1, parallel::detectCores() - 2)


# ---------------------------------------------------------------------------------------
# command line overrides
# ---------------------------------------------------------------------------------------
#
# The defaults above apply when a script is run with no arguments. run_all_inference.R
# uses these flags to give each method process its own core budget:
#
#   --cores N            cluster size for this process
#   --sizes 50,150       restrict to these sample sizes
#   --max-per-group N    trios per group, for smoke tests. Must match what built the
#                        selection cache or the cache lookup misses and it recomputes.
#
# The confounder-structure runs (run_structure_sims.R) add three more. These point the whole
# stage at a different simulation without editing this file, so the four structures can be
# run from one checkout:
#
#   --sim-file PATH        which simulated_trios*.RData to read
#   --out-dir PATH         where checkpoints and the selection cache go
#   --filter-int-child 0|1 override filter_int_child for that run
#   --rerun-inference 0|1  redo groups that already have a checkpoint (needed after a
#                          schema change, e.g. a new arm)
#
# --out-dir is not cosmetic. The selection cache is named selection_group_n<size>.RData
# with no other key, so two runs over different simulated data sharing an out.dir would
# fight over one cache file. selection.cache.mismatch() would catch it -- cov.names differ
# -- and silently recompute every time, turning a cache into a 40-minute tax per run.
#
# --sim-file and --out-dir are read at the TOP of this file, not here, because out.dir is
# created before this point.
if (!is.null(inference.opt("--cores"))) {
    n.cores <- max(1L, as.integer(inference.opt("--cores")))
}
if (!is.null(inference.opt("--sizes"))) {
    sample.sizes <- as.numeric(strsplit(inference.opt("--sizes"), ",")[[1]])
}
if (!is.null(inference.opt("--max-per-group"))) {
    max.per.group <- as.integer(inference.opt("--max-per-group"))
}
if (!is.null(inference.opt("--filter-int-child"))) {
    filter_int_child <- as.logical(as.integer(inference.opt("--filter-int-child")))
}
# --rerun-inference 1 redoes groups that already have a checkpoint. Needed whenever the
# result SCHEMA changes -- a new arm, say -- because a checkpoint is trusted on sight and an
# old one would otherwise be combined with new ones and silently produce a master with
# missing columns. Prefer this to deleting checkpoints by hand: it only overwrites the
# groups this invocation actually recomputes.
# --bootstrap 0 turns off the MRGN bootstrap. See mrgn.bootstrap above: it changes no
# model call and no column, it only leaves the boot.* block NA, and it is the difference
# between minutes and hours per group.
if (!is.null(inference.opt("--bootstrap"))) {
    mrgn.bootstrap <- as.logical(as.integer(inference.opt("--bootstrap")))
}
if (!is.null(inference.opt("--rerun-inference"))) {
    rerun.inference <- as.logical(as.integer(inference.opt("--rerun-inference")))
}

# ---- MR-GGI ----
# Edge called at raw p < mrggi.alpha, matching how mrpc.alpha and selection.alpha are
# applied per trio. cor.thr = 0 is the package default: it decides which gene pairs are
# tested at all, and at 0 every pair with a nonzero correlation is.
mrggi.alpha   <- 0.05
mrggi.cor.thr <- 0

# Which covariate sets MR-GGI is run against, one MRggi() call per arm per trio.
#
#   none    the bare trio, no covariates in y. What the pipeline did before, kept so the
#           earlier numbers stay reproducible against the new ones.
#   truth   trio + K + that trio's own U block
#   CSq     trio + K + the CS-q selection
#   CSa     trio + K + the CS-alpha selection
#
# WHAT THE COVARIATES ACTUALLY DO, because it is not what it looks like: MRggi's estimator
# is strictly pairwise. .TSLS() reads only y[,i], y[,j], X[[i]] and X[[j]], so the T1-T2
# Wald ratio is bit-identical across all four arms. Measured on one trio: Bg1g2 = 0.808,
# p = 0.000 both with and without covariates in y. What the extra columns change is
# FDR_Bg1g2 -- MRggi adjusts each g1's p-values across that gene's pairs, so the T1-T2
# p-value is now corrected for T1's pairs with every covariate as well. That multiplicity
# correction is the whole difference between the arms, and it is why the arms are worth
# having rather than a mistake in the setup. See MRGGI_METHODS.md section 4.
#
# COST. cor.thr = 0 means every pair is computed, so an arm with k covariates costs
# O((k+2)^2) TSLS fits. Measured, per 300-trio group:
#
#   arm      n=50    n=150   n=300   n=670    n=1000
#   none     ~0      ~0      ~0      ~0       ~0
#   truth    6 min   7 min   10 min  11 min   10 min
#   CSq      ~0      ~0      0.4 min 4 min    8 min
#   CSa      67 min  94 min  111 min 134 min  159 min
#
# CS-alpha is 90% of the total (9.4 h of 10.4 h across all groups) because it selects a
# median of 82-106 covariates, i.e. ~5,800 pairs per trio. That is why apply_mrggi.R now
# builds a cluster instead of running single threaded.
mrggi.arms    <- c("none", "truth", "CSq", "CSa", "CSi")

# Passed to MRggi(p.adjust.method = ), and IGNORED BY THE PACKAGE. Its body reads
#
#     FDR.idx = p.adjust(pval.idx, method = p.adjust.methods)
#
# with a trailing "s" -- base R's vector of every method name -- so match.arg() silently
# takes the first element and the correction is ALWAYS holm, whatever is passed. Verified:
# p.adjust(c(.01,.02,.03,.04), method = p.adjust.methods) returns 0.04 0.06 0.06 0.06,
# which is holm and not bonferroni. The argument is set here for the day the package is
# fixed; until then read mrggi.<arm>.FDR.T1T2 as holm-adjusted, and recompute from
# mrggi.<arm>.p.T1T2 in the analysis stage if a different correction is wanted.
mrggi.p.adjust <- "bonferroni"

# Minimum first-stage F for the exposure's instrument. MRggi's p-value with a single
# instrument reduces to the instrument->OUTCOME t-statistic, so the exposure's first stage
# cancels out of the test entirely and nothing stops it reporting p ~ 0 for a ratio whose
# denominator is ~0. Measured on a known-null trio: reported B = -38.42 at p = 0 with a
# first-stage F of 0.1. F > 10 is the conventional weak-instrument rule (Staiger & Stock).
# See mrggi_feasibility.R sections 5-6.
mrggi.min.F   <- 10
