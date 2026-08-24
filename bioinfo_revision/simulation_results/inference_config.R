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

results.root  <- "bioinfo_revision/simulation_results"

sim.data.file <- "bioinfo_revision/simulation/simulated_data/simulated_trios.RData"
clinical.file <- "./GTEx/data/kclist_top5_tiss.RData"
out.dir       <- file.path(results.root, "data")     # every generated .RData / .csv
log.dir       <- file.path(results.root, "logs")
tables.dir    <- file.path(results.root, "tables")   # used by results_scripts/

# out.dir was the script folder and so always existed. It is a generated directory now, and
# the first thing to touch it is a save() that would fail on a fresh clone.
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# ---- per-method inference settings ----
n.bootstrap     <- 1000    # bootstrap replicates per MRGN fit
mrpc.timeout    <- 120     # seconds; MRPC has taken hours on trios with many confounders
gmac.nperm      <- 1000
selection.alpha <- 0.05    # significance cutoff for the GMAC mediation calls

# Which confounder sets MRPC is run against. Both were attempted originally; CS-alpha is
# now off because it does not finish.
#
# CS-alpha hands MRPC alpha x pool false-positive covariates by construction -- a median
# of 95 at n = 150 -- and MRPC cannot fit that many. Measured on the completed groups:
#
#   n = 50   CS-q  2 confounders,  0% timeout, median 0.0 s   |  group took 48 min
#   n = 150  CS-q  2 confounders,  0% timeout, median 0.0 s   |  group took 8.8 h
#   n = 50   CS-a 82 confounders,  0% timeout, median 7.8 s
#   n = 150  CS-a 95 confounders, 79% timeout, ~106 s/trio
#
# The n = 50 CS-alpha fits are fast only because 82 covariates on 50 observations is
# rank deficient and bails early. By n = 150 the fits are estimable, MRPC genuinely
# attempts a PC algorithm over ~98 nodes, and the timeout rate jumps to 79% -- heading to
# ~100% at larger n, i.e. ~10 h per group to produce a column of NAs. CS-q costs nothing
# because it selects a median of 2 covariates.
#
# Groups already run with both arms (n = 50, n = 150) keep their real CS-alpha results.
# Groups run under this setting record CS-alpha as not attempted, with the reason in
# mrpc.CSa.error, rather than silently NA.
mrpc.arms <- c("CSq")      # c("CSq", "CSa") to attempt both again

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
inference.args <- commandArgs(trailingOnly = TRUE)
inference.opt <- function(flag, default = NULL) {
    i <- match(flag, inference.args)
    if (is.na(i) || i == length(inference.args)) return(default)
    inference.args[i + 1]
}
if (!is.null(inference.opt("--cores"))) {
    n.cores <- max(1L, as.integer(inference.opt("--cores")))
}
if (!is.null(inference.opt("--sizes"))) {
    sample.sizes <- as.numeric(strsplit(inference.opt("--sizes"), ",")[[1]])
}
if (!is.null(inference.opt("--max-per-group"))) {
    max.per.group <- as.integer(inference.opt("--max-per-group"))
}

# ---- MR-GGI ----
# Edge called at raw p < mrggi.alpha, matching how mrpc.alpha and selection.alpha are
# applied per trio. cor.thr = 0 is the package default: with a single gene pair the screen
# only decides whether the trio is tested at all.
mrggi.alpha   <- 0.05
mrggi.cor.thr <- 0

# Minimum first-stage F for the exposure's instrument. MRggi's p-value with a single
# instrument reduces to the instrument->OUTCOME t-statistic, so the exposure's first stage
# cancels out of the test entirely and nothing stops it reporting p ~ 0 for a ratio whose
# denominator is ~0. Measured on a known-null trio: reported B = -38.42 at p = 0 with a
# first-stage F of 0.1. F > 10 is the conventional weak-instrument rule (Staiger & Stock).
# See mrggi_feasibility.R sections 5-6.
mrggi.min.F   <- 10
