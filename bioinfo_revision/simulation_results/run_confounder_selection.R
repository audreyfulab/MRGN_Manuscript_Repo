# Confounder selection for the revised simulation, run as its own stage.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
# Selection is the most expensive step in the pipeline and does not depend on which
# method is fitted afterwards, so it runs once here and is cached per sample-size group.
# The apply_*.R scripts then load the cache instead of recomputing it.
#
# For each sample-size group this writes selection_group_n<size>.RData, and at the end a
# single selection_results.RData / .csv scoring every trio's selected set against the
# confounders that actually generated it.
#
#   Rscript bioinfo_revision/simulation_results/run_confounder_selection.R
#
# Running it again is a no-op unless rerun.selection is TRUE: get.selection() loads each
# cache, checks it against the current request, and only recomputes on a mismatch.

library(MRGN)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

# Paths, the CS-q / CS-alpha settings, sample.sizes and max.per.group all come from
# inference_config.R, so this stage and the apply_*.R scripts cannot drift apart. CS-q is
# the q-value method at FDR 5%; CS-alpha thresholds the same regression p-values at alpha
# with no correction. Both come out of the one get.conf.trios() call.
#
# rerun.selection (also from the config) recomputes every group and overwrites its cache.


# ---------------------------------------------------------------------------------------
# per-trio scoring of a group's two selected sets
# ---------------------------------------------------------------------------------------

score.group <- function(datasets, sel) {

    cov.names <- colnames(group.cov.pool(datasets))

    rows <- lapply(seq_along(datasets), function(i) {
        dat <- datasets[[i]]$data
        params <- datasets[[i]]$params
        index <- params$dataset
        # the U block is the truth; W is an intermediate and Z a common child, so a
        # selection containing either is a false positive with a specific meaning
        true.confs <- own.block.names(colnames(dat), "U", index)

        csq.names <- cov.names[sel$CS.q$sig.asso.covs[[i]]]
        csa.names <- cov.names[sel$CS.alpha$sig.asso.covs[[i]]]

        as.data.frame(c(
            list(dataset = index),
            as.list(params[, setdiff(colnames(params), "dataset")]),
            list(truth.model = unname(TRUTH.LABEL[params$model]),
                 n.true.confs = length(true.confs),
                 true.confs = paste(true.confs, collapse = ";")),
            prefixed(score.selection(csq.names, index, true.confs), "CSq"),
            prefixed(score.selection(csa.names, index, true.confs), "CSa"),
            list(CSq.filter_int_child = sel$CS.q$filter_int_child,
                 CSa.filter_int_child = sel$CS.alpha$filter_int_child,
                 selection.time.seconds = sel$time.seconds)),
            stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}


# ---------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------

run.all.selection <- function(sizes = sample.sizes, rerun = rerun.selection, verbose = TRUE) {

    sim_data <- loadRData(file = sim.data.file)
    if (verbose) cat("loaded", length(sim_data), "simulated datasets\n")

    group.of <- sapply(sim_data, function(x) x$params$sample.size)
    if (is.null(sizes)) sizes <- sort(unique(group.of))

    scored <- list()
    for (size in sizes) {
        idx <- which(group.of == size)
        if (!is.null(max.per.group)) idx <- head(idx, max.per.group)
        datasets <- sim_data[idx]

        if (verbose) cat("\n=== sample size", size, "|", length(datasets), "trios ===\n")
        start <- Sys.time()
        sel <- get.selection(datasets, out.dir = out.dir, sim.file = sim.data.file,
                             rerun = rerun,
                             # never let a subset run overwrite a full-group cache
                             save = is.null(max.per.group),
                             verbose = verbose,
                             selection_fdr = selection_fdr, filter_fdr = filter_fdr,
                             alpha = alpha, filter_int_child = filter_int_child)

        scored[[as.character(size)]] <- score.group(datasets, sel)
        if (verbose) {
            n.csq <- sapply(sel$CS.q$sig.asso.covs, length)
            n.csa <- sapply(sel$CS.alpha$sig.asso.covs, length)
            cat("  CS-q     median", median(n.csq), "selected/trio (range",
                min(n.csq), "-", max(n.csq), ")\n")
            cat("  CS-alpha median", median(n.csa), "selected/trio (range",
                min(n.csa), "-", max(n.csa), ")\n")
            cat("  filter_int_child applied:", sel$CS.q$filter_int_child, "|",
                round(as.numeric(difftime(Sys.time(), start, units = "mins")), 1),
                "min for this group\n")
        }
    }

    selection.results <- do.call(rbind, scored)
    selection.results <- selection.results[order(selection.results$dataset), ]
    row.names(selection.results) <- NULL

    if (is.null(max.per.group)) {
        base::save(selection.results, file = file.path(out.dir, "selection_results.RData"))
        write.csv(selection.results, file.path(out.dir, "selection_results.csv"),
                  row.names = FALSE)
        if (verbose) cat("\nwrote selection_results.RData / .csv |",
                         nrow(selection.results), "rows\n")
    }
    return(invisible(selection.results))
}


selection.results <- run.all.selection()
