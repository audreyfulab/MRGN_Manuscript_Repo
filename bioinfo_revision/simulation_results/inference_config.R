# Shared configuration for the inference stage.
#
# Sourced by apply_mrgn.R, apply_mrpc.R, apply_gmac.R and run_all_inference.R, so every
# method sees identical paths and settings. Change a value here, not in one of the apply
# scripts, or the three methods stop being comparable.
#
# Paths are relative to the repository root; every script in this stage expects to be run
# from there.

sim.data.file <- "bioinfo_revision/simulation/simulated_data/simulated_trios.RData"
clinical.file <- "./GTEx/data/kclist_top5_tiss.RData"
out.dir       <- "bioinfo_revision/simulation_results"
log.dir       <- file.path(out.dir, "logs")

# ---- per-method inference settings ----
n.bootstrap     <- 1000    # bootstrap replicates per MRGN fit
mrpc.timeout    <- 120     # seconds; MRPC has taken hours on trios with many confounders
gmac.nperm      <- 1000
selection.alpha <- 0.05    # significance cutoff for the GMAC mediation calls

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
