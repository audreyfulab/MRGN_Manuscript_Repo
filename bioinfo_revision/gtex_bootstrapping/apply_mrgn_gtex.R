# Apply MRGN, with the edge-probability bootstrap, to the GTEx Whole Blood trios.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/gtex_bootstrapping/apply_mrgn_gtex.R
#
# Reads GTEx/data/data.with.PCs.WholeBlood.RData -- 3,248 trios, 670 samples each, already
# carrying their selected confounders -- and writes one row per trio to
# bioinfo_revision/reports/gtex_bootstrapping/.
#
# INPUT LAYOUT. data.sets is an unnamed list of data.frames, each laid out exactly as
# MRGN::infer.trio() expects:
#
#   col 1     the variant          (chr..._b38)
#   col 2     the cis gene    T1   (ENSG...)
#   col 3     the trans gene  T2   (ENSG...)
#   col 4+    covariates: the clinical knowns (pcr, platform, sex) then the selected PCs
#
# 12 to 57 columns, median 36, so 9 to 54 covariates per trio. infer.trio() takes no
# confounder-count argument -- it reads the first three columns as the trio and treats the
# rest as covariates -- so the frames go in unmodified.
#
# NO GROUND TRUTH HERE. Unlike the simulation stage there is no generating model to score
# against, so there is no `correct` column and mrgn.fields() is not reused: it exists to
# compare a call against a known truth. What this stage produces is the call itself, the
# bootstrap support for it, and enough provenance to join back to the trio.
#
# ---------------------------------------------------------------------------------------
# PARALLELISM: OVER TRIOS, NOT OVER RESAMPLES
# ---------------------------------------------------------------------------------------
#
# boostrap_edge_probabilities() can take a cluster and will spread its `number_of_samples`
# replicates across it. DO NOT USE IT THAT WAY HERE. Each replicate is one infer.trio() call
# on 670 rows -- a few milliseconds -- so dispatching 1,000 of them through parLapply costs
# more in socket round-trips than it saves in arithmetic.
#
# That was measured, on the confounder-structure runs of 2026-08-26: three concurrent
# 4-core clusters bootstrapping this way held 1.9 core-equivalents of 14 busy, with the disk
# at 0.3% and 8 GB of RAM free. It was not compute-bound; it was waiting on the cluster.
#
# So the cluster parallelises the OUTER loop -- one worker takes one whole trio and runs its
# 1,000 replicates serially in-process (cl = NULL inside apply.mrgn). One dispatch per trio
# instead of 1,000. Measured serially at 6.9-9.6 s per trio, so ~7.2 h single-threaded and
# roughly 35-40 min across a full cluster.

library(MRGN)
library(parallel)

source("bioinfo_revision/simulation_results/inference_utils.R")

# ---------------------------------------------------------------------------------------
# configuration and command line
# ---------------------------------------------------------------------------------------
gtex.args <- commandArgs(trailingOnly = TRUE)
gtex.opt  <- function(flag, default = NULL) {
    i <- match(flag, gtex.args)
    if (is.na(i) || i == length(gtex.args)) return(default)
    gtex.args[i + 1]
}

gtex.file  <- gtex.opt("--gtex-file", "GTEx/data/data.with.PCs.WholeBlood.RData")
out.dir    <- gtex.opt("--out-dir",   "bioinfo_revision/reports/gtex_bootstrapping")
tissue     <- gtex.opt("--tissue",    "WholeBlood")

n.bootstrap <- as.integer(gtex.opt("--bootstrap", 1000))
n.cores     <- as.integer(gtex.opt("--cores", max(1L, parallel::detectCores() - 2L)))
max.trios   <- if (is.null(gtex.opt("--max-trios"))) NULL else
                   as.integer(gtex.opt("--max-trios"))
# Trios per checkpoint. The run is long enough that losing it to a closed window would
# hurt, and small enough that per-chunk files stay tidy. A chunk that already exists is
# skipped, so an interrupted run resumes at chunk granularity.
chunk.size  <- as.integer(gtex.opt("--chunk-size", 200))
rerun       <- identical(gtex.opt("--rerun", "0"), "1")

chunk.dir <- file.path(out.dir, "chunks")
dir.create(chunk.dir, recursive = TRUE, showWarnings = FALSE)

# Per-trio progress, written by the workers. Truncated at the start of a run so it always
# describes the run in progress rather than accumulating across restarts; the chunk files
# are the durable record, this is the live view.
progress.log <- file.path(out.dir, "progress.log")

# ---------------------------------------------------------------------------------------
# per-trio row
# ---------------------------------------------------------------------------------------
#
# Returns a one-row data.frame whatever happens. A trio whose inference fails still gets its
# identity columns and an error string, so the output always has as many rows as the input
# and a failure is visible rather than absent.
gtex.row <- function(trio, index, n.samples.boot) {
    # Force the argument before anything can carry it onto a worker. R would otherwise keep
    # it as a promise pointing at the caller's `n.bootstrap`, and apply.mrgn() does not
    # touch number_of_samples until it reaches the bootstrap -- by which point the promise
    # is being evaluated in a worker where that name does not exist. It surfaces as
    # "object 'n.bootstrap' not found" in bootstrap.error, with the model call unaffected,
    # so every row looks fine except that the whole bootstrap is silently missing.
    force(n.samples.boot)
    nms  <- colnames(trio)
    covs <- nms[-(1:3)]
    # The clinical knowns are carried for every trio; everything else is a selected PC.
    known <- intersect(covs, c("pcr", "platform", "sex"))
    pcs   <- setdiff(covs, known)

    fit <- safely(apply.mrgn(trio, bootstrap = TRUE,
                             number_of_samples = n.samples.boot,
                             cl = NULL, verbose = FALSE))
    r <- fit$value
    b <- if (is.null(r)) NULL else r$bootstrap
    ind <- if (is.null(b)) NULL else b$indicator.means

    # Per-trio progress, emitted from the WORKER.
    #
    # Two destinations, because neither alone is enough. cat() to stdout reaches the console
    # only because the cluster is built with outfile = "" -- without that, PSOCK workers
    # discard their output and this line would vanish. And the console is not readable from
    # outside the window, so the same line is appended to progress.log.
    #
    # Appends from 12 workers race, but each write is one short line opened in append mode,
    # which the OS serialises well enough at this size; a garbled line costs a log entry and
    # nothing else, since the results themselves come back through the cluster.
    msg <- sprintf("trio %5d complete! inferred as %-6s | boot %-6s %s | %4.1fs",
                   index,
                   if (is.null(r)) "FAILED" else as.character(r$model),
                   if (is.null(b)) "--" else as.character(b$boot.model),
                   if (is.null(b) || is.null(r)) " "
                   else if (identical(as.character(b$boot.model), as.character(r$model))) "(agrees)"
                   else "(differs)",
                   if (is.null(r)) NA_real_ else r$time.seconds + r$bootstrap.time.seconds)
    cat(msg, "\n", sep = "")
    try(cat(msg, "\n", sep = "", file = progress.log, append = TRUE), silent = TRUE)

    data.frame(
        trio.index    = index,
        tissue        = tissue,
        variant       = nms[1],
        cis.gene      = nms[2],
        trans.gene    = nms[3],
        n.samples     = nrow(trio),
        n.covariates  = length(covs),
        n.known       = length(known),
        n.pcs         = length(pcs),
        # the call
        model         = if (is.null(r)) NA_character_ else as.character(r$model),
        time.seconds  = if (is.null(r)) NA_real_ else r$time.seconds,
        # bootstrap support. b11 = V1->T1, b12 = T1->T2, b21 = V1->T2, b22 = T2->T1,
        # matching the indicator order infer.trio() returns.
        boot.model         = if (is.null(b)) NA_character_ else as.character(b$boot.model),
        boot.agrees        = if (is.null(b) || is.null(r)) NA
                             else identical(as.character(b$boot.model), as.character(r$model)),
        boot.min.edge.prob = if (is.null(b)) NA_real_ else b$min.edge.prob,
        boot.p.V1T1        = if (is.null(ind)) NA_real_ else unname(ind[1]),
        boot.p.T1T2        = if (is.null(ind)) NA_real_ else unname(ind[2]),
        boot.p.V1T2        = if (is.null(ind)) NA_real_ else unname(ind[3]),
        boot.p.T2T1        = if (is.null(ind)) NA_real_ else unname(ind[4]),
        boot.n.requested   = if (is.null(b)) NA_integer_ else b$number_of_samples,
        # A resample that loses every minor allele leaves V1 constant and is unusable; those
        # are dropped rather than allowed to poison the means, and the count is reported so
        # a trio whose bootstrap rested on few usable draws is identifiable.
        boot.n.used        = if (is.null(b)) NA_integer_ else b$n.used,
        boot.n.dropped     = if (is.null(b)) NA_integer_ else b$n.dropped,
        bootstrap.time.seconds = if (is.null(r)) NA_real_ else r$bootstrap.time.seconds,
        bootstrap.error    = if (is.null(r)) NA_character_ else r$bootstrap.error,
        error              = fit$error,
        stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------------------
cat("=== MRGN + bootstrap on GTEx ", tissue, " ===\n", sep = "")
cat("  input     : ", gtex.file, "\n", sep = "")
cat("  output    : ", out.dir, "\n", sep = "")
cat("  bootstrap : ", n.bootstrap, " replicates per trio\n", sep = "")
cat("  cores     : ", n.cores, " (parallel over TRIOS; each trio bootstraps in-process)\n",
    sep = "")

if (!file.exists(gtex.file)) stop("no GTEx data at ", gtex.file)
env <- new.env(parent = emptyenv())
load(gtex.file, envir = env)
if (!"data.sets" %in% ls(env)) {
    stop(gtex.file, " holds no object named 'data.sets' (found: ",
         paste(ls(env), collapse = ", "), ")")
}
datasets <- get("data.sets", envir = env)
rm(env)

if (!is.null(max.trios)) datasets <- datasets[seq_len(min(max.trios, length(datasets)))]
n.trios <- length(datasets)
cat("  trios     : ", n.trios, "\n", sep = "")
cat("  chunks    : ", ceiling(n.trios / chunk.size), " of up to ", chunk.size, "\n\n", sep = "")

cat(sprintf("[%s] starting\n", format(Sys.time(), "%H:%M:%S")), file = progress.log)

# outfile = "" is what makes the workers' per-trio cat() visible. A PSOCK cluster built
# without it sends worker stdout to a null device, so every progress line would be
# discarded and the run would look silent for a whole chunk at a time.
cl <- parallel::makeCluster(n.cores, outfile = "")
on.exit(parallel::stopCluster(cl), add = TRUE)
invisible(parallel::clusterEvalQ(cl, library(MRGN)))
# By name rather than by re-sourcing inference_utils.R on each worker: that file sources
# adapted_GMAC_func/* at its top level, which is a heavy and pointless dependency for a
# worker that only ever calls apply.mrgn().
parallel::clusterExport(cl, c("apply.mrgn", "boostrap_edge_probabilities", "safely",
                              "gtex.row", "tissue", "n.bootstrap", "progress.log"),
                        envir = environment())

starts <- seq(1, n.trios, by = chunk.size)
started.all <- Sys.time()

for (ci in seq_along(starts)) {
    from <- starts[ci]
    to   <- min(from + chunk.size - 1, n.trios)
    path <- file.path(chunk.dir, sprintf("chunk_%04d.RData", ci))

    if (file.exists(path) && !rerun) {
        cat(sprintf("  chunk %3d/%3d  trios %5d-%-5d  cached, skipping\n",
                    ci, length(starts), from, to))
        next
    }

    t0 <- Sys.time()
    idx <- from:to
    # clusterMap, not parLapply over indices: a closure referring to `datasets[[i]]` would
    # capture the whole 755 MB list and serialise it to every worker on every chunk.
    # clusterMap sends each worker only the trio it is about to run.
    rows <- parallel::clusterMap(cl, function(tr, i) gtex.row(tr, i, n.bootstrap),
                                 datasets[idx], idx, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    chunk <- do.call(rbind, rows)
    save(chunk, file = path)

    mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    done <- sum(!is.na(chunk$model))
    eta  <- mins * (length(starts) - ci)
    cat(sprintf("  chunk %3d/%3d  trios %5d-%-5d  %5.1f min | %3d/%3d inferred | ETA %.0f min\n",
                ci, length(starts), from, to, mins, done, nrow(chunk), eta))
}

# ---------------------------------------------------------------------------------------
# combine
# ---------------------------------------------------------------------------------------
cat("\ncombining:\n")
files <- sort(list.files(chunk.dir, pattern = "^chunk_[0-9]+[.]RData$", full.names = TRUE))
gtex.mrgn.results <- do.call(rbind, lapply(files, function(f) {
    e <- new.env(parent = emptyenv()); load(f, envir = e); get("chunk", envir = e)
}))
gtex.mrgn.results <- gtex.mrgn.results[order(gtex.mrgn.results$trio.index), ]
row.names(gtex.mrgn.results) <- NULL

save(gtex.mrgn.results, file = file.path(out.dir, "gtex_mrgn_bootstrap.RData"))
utils::write.csv(gtex.mrgn.results, file.path(out.dir, "gtex_mrgn_bootstrap.csv"),
                 row.names = FALSE)

cat(sprintf("  gtex_mrgn_bootstrap: %d chunks, %d rows, %d columns\n",
            length(files), nrow(gtex.mrgn.results), ncol(gtex.mrgn.results)))

failed <- sum(is.na(gtex.mrgn.results$model))
if (failed > 0) cat(sprintf("  %d trios have no model call -- see the error column\n", failed))
bad.boot <- sum(!is.na(gtex.mrgn.results$bootstrap.error))
if (bad.boot > 0) cat(sprintf("  %d trios inferred but their bootstrap failed\n", bad.boot))

cat("\n  inferred model distribution:\n")
print(table(gtex.mrgn.results$model, useNA = "ifany"))
cat("\n  bootstrap agrees with the point call: ",
    sprintf("%.1f%%", 100 * mean(gtex.mrgn.results$boot.agrees, na.rm = TRUE)), "\n", sep = "")
cat("\nall trios finished in ",
    round(as.numeric(difftime(Sys.time(), started.all, units = "mins")), 1), " min\n", sep = "")
