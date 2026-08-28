# Add a CS-i block to the cached confounder selections, in place.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/backfill_csi.R
#
# ---------------------------------------------------------------------------------------
# WHY THIS IS A BACKFILL RATHER THAN A RE-RUN
# ---------------------------------------------------------------------------------------
#
# CS-q, CS-alpha and CS-i are three thresholdings of ONE matrix: get.conf.trios()'s
# reg.pvalues, a 2-df F test of lm(covariate ~ T1 + T2) for every (trio, covariate) cell.
# That matrix is what costs ~40 min per group to compute, and select.confounders() already
# stores it -- verified present at 300 x 8,125..8,838 in all five selection_group_n*.RData.
#
# So the whole CS-i arm is recoverable from the cache for ~18 s per group. Re-running
# select.confounders() would cost ~3.3 h across the five groups to reproduce two blocks
# that already exist and add a third that does not need the regression at all.
#
# ---------------------------------------------------------------------------------------
# WHAT THIS TOUCHES, AND WHAT IT MUST NOT
# ---------------------------------------------------------------------------------------
#
# It adds sel$CS.i and rewrites the cache file. CS.q and CS.alpha must come through
# untouched, and that is ASSERTED rather than assumed: their sig.asso.covs are compared
# before and after, and the script stops without writing if either moved. A silent change
# there would repoint every existing MRGN/MRPC/MR-GGI arm at a different confounder set
# while leaving the results files looking unchanged.
#
# Re-running this script is safe. A cache that already has CS.i is recomputed and rewritten
# to the same content, not appended to.

library(MRGN)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

cat("=== CS-i backfill ===  out.dir:", out.dir, "\n")

sim_data <- loadRData(file = sim.data.file)
group.of <- sapply(sim_data, function(x) x$params$sample.size)
sizes <- if (is.null(sample.sizes)) sort(unique(group.of)) else sample.sizes

for (size in sizes) {
    path <- file.path(out.dir, sprintf("selection_group_n%d.RData", size))
    if (!file.exists(path)) {
        cat("\n n =", size, ": no selection cache at", path, "-- skipping\n")
        next
    }
    cat("\n n =", size, "\n")

    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    cache <- env$cache
    sel   <- cache$sel

    if (is.null(sel$selection$reg.pvalues)) {
        stop("n = ", size, ": cache has no selection$reg.pvalues, so CS-i cannot be ",
             "derived from it. Re-run select.confounders() for this group.")
    }

    # The trio and covariate blocks CS-i has to be assembled against. Rebuilt from the
    # simulation exactly as run.method.groups() does, so the CS-i trios.with.confs frames
    # are constructed the same way CS-q's and CS-alpha's were.
    idx <- which(group.of == size)
    if (!is.null(max.per.group)) idx <- head(idx, max.per.group)
    datasets   <- sim_data[idx]
    cov.pool   <- group.cov.pool(datasets)
    known.conf <- group.known.conf(datasets)
    # Built exactly as get.selection() builds it -- same expression, no naming step. The
    # trios list is unnamed there, and CS-q's/CS-alpha's blocks inherited that, so naming
    # it here would give CS-i a different shape from the two it sits beside.
    trios <- lapply(datasets, function(x) x$data[, 1:3])

    reg.p <- as.matrix(sel$selection$reg.pvalues)
    stopifnot(nrow(reg.p) == length(trios), ncol(reg.p) == ncol(cov.pool))

    # Guard the cache against the datasets having moved under it. cov.names is what
    # selection.cache.mismatch() keys on, and a mismatch here means this cache was built
    # from a different simulation than the one just loaded.
    if (!identical(colnames(cov.pool), colnames(sel$selection$reg.pvalues))) {
        stop("n = ", size, ": cov.pool column names do not match the cached reg.pvalues. ",
             "This cache was built from a different simulation -- see --sim-file/--out-dir ",
             "in inference_config.R.")
    }

    before <- list(q = sel$CS.q$sig.asso.covs, a = sel$CS.alpha$sig.asso.covs)

    t0 <- Sys.time()
    csi.covs <- csi.sig.asso.covs(reg.p, selection_fdr = selection_fdr,
                                  trio.names = names(trios))
    sel$CS.i <- assemble.selection(csi.covs, trios, cov.pool, known.conf,
                                   setting = "CS-i",
                                   filter.applied = sel$CS.q$filter_int_child,
                                   elapsed = sel$time.seconds)
    sel$CS.i$adjust_by <- "individual"
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    # The two existing blocks must be bit-for-bit what they were.
    stopifnot(identical(before$q, sel$CS.q$sig.asso.covs),
              identical(before$a, sel$CS.alpha$sig.asso.covs))

    n.sel <- lengths(csi.covs)
    cat(sprintf("   CS-i in %.0f s | selected/trio: mean %.2f, median %d, range %d-%d | %d trios with none\n",
                elapsed, mean(n.sel), as.integer(median(n.sel)), min(n.sel), max(n.sel),
                sum(n.sel == 0)))
    cat(sprintf("   for comparison  CS-q mean %.2f | CS-alpha mean %.2f\n",
                mean(lengths(sel$CS.q$sig.asso.covs)),
                mean(lengths(sel$CS.alpha$sig.asso.covs))))

    cache$sel <- sel
    save(cache, file = path)
    cat("   rewrote", basename(path), "\n")

    rm(env, cache, sel, datasets, cov.pool, trios, reg.p); invisible(gc())
}

cat("\ndone.\n")
