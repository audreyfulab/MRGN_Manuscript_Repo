# Inference utilities for the revised trio simulation.
#
# Sourced by updated_simulation_inference.R. Each apply.* function runs one method on a
# SINGLE trio and attaches the time taken by the inference itself, in seconds, so the
# three methods can be compared directly. Confounder selection is a separate step: MRGN
# and MRPC share MRGN::get.conf.trios() via select.confounders(), while GMAC performs its
# own selection inside run.gmac.all().
#
# Paths below are relative to the repository root, matching updated_simulation_inference.R.
# gmac_one_trio.R and GMAC_moded/R/GMAC.R both define Indirect(), my.solve(),
# nominal.pfun() and get.beta.change() with identical bodies, so sourcing both is
# deliberate and harmless -- gmacOneTrio() lives in the first, gmac() in the second.
source("adapted_GMAC_func/gmac_one_trio.R")
source("adapted_GMAC_func/GMAC_moded/R/GMAC.R")


# ---------------------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------------------
#
# These live here rather than in updated_simulation_inference.R because
# run_confounder_selection.R needs them too.

# the generating model of each scenario, as the label the methods report. The simulation
# always puts the cis gene in T1, so the truth is always the ".1" variant. An explicit
# lookup rather than MRGN::convert.truth(), which maps by sorted position and silently
# mislabels when the input does not contain all five models.
TRUTH.LABEL <- c(model0 = "M0.1", model1 = "M1.1", model2 = "M2.1",
                 model3 = "M3", model4 = "M4")

coarse.model <- function(x) sub("\\..*$", "", x)

# names of the columns of one dataset that belong to that dataset's own confounder blocks.
# name.trio.columns() suffixed every simulated confounder with the dataset index, so a
# column can be attributed to its trio and its block from the name alone.
own.block.names <- function(nms, block, index) {
    grep(paste0("^", block, "[0-9]+\\.", index, "$"), nms, value = TRUE)
}

# trio + known confounders + that trio's true confounders. The U variables are the true
# confounders; W is an intermediate and Z a common child, so both are excluded.
ground.truth.input <- function(dat, K_n, index) {
    nms <- colnames(dat)
    keep <- c(nms[1:3],
              if (K_n > 0) nms[3 + seq_len(K_n)] else character(0),
              own.block.names(nms, "U", index))
    dat[, keep, drop = FALSE]
}

# score one selected confounder set against the truth for that trio
score.selection <- function(selected, index, true.confs) {
    own.suffix <- paste0("\\.", index, "$")
    false.pos <- setdiff(selected, true.confs)
    list(n.selected        = length(selected),
         n.tp              = length(intersect(selected, true.confs)),
         n.fp              = length(false.pos),
         n.fn              = length(setdiff(true.confs, selected)),
         # false positives borrowed from a different trio's confounder block, as opposed
         # to this trio's own intermediate/common child
         n.fp.other.trio   = sum(!grepl(own.suffix, false.pos)),
         has.common.child  = length(own.block.names(selected, "Z", index)) > 0,
         has.intermediate  = length(own.block.names(selected, "W", index)) > 0,
         selected          = paste(selected, collapse = ";"))
}

# prefix a named list so it can be spliced into the wide row
prefixed <- function(x, prefix) {
    names(x) <- paste0(prefix, ".", names(x))
    x
}

# every apply.* call is wrapped: at n = 50 a trio can carry up to 44 confounders, so
# rank-deficient fits are expected and must not take down a whole run
safely <- function(expr) {
    tryCatch(list(value = expr, error = NA_character_),
             error = function(e) list(value = NULL, error = conditionMessage(e)))
}

# the covariate pool and known-confounder block for one sample-size group. Every trio
# contributes its own U/W/Z columns; the K block is shared by all trios in the group (the
# same clinical covariates), so it is passed separately rather than repeated in the pool.
group.cov.pool <- function(datasets) {
    do.call(cbind.data.frame,
            lapply(datasets, function(x) {
                x$data[, (4 + x$params$K_n):ncol(x$data), drop = FALSE]
            }))
}

group.known.conf <- function(datasets) {
    K_n <- datasets[[1]]$params$K_n
    if (K_n > 0) datasets[[1]]$data[, 3 + seq_len(K_n), drop = FALSE] else NULL
}


# ---------------------------------------------------------------------------------------
# confounder selection (MRGN and MRPC)
# ---------------------------------------------------------------------------------------

# Confounder selection for one group of trios via MRGN::get.conf.trios(), returning both
# selection settings side by side:
#
#   CS-q     - q-value method controlling FDR at 5% (adjust_by = "all", selection_fdr = 0.05)
#   CS-alpha - no multiple-testing correction, each test at a fixed type I error rate
#              alpha < 0.01 (equivalent to adjust_by = "none", alpha = 0.01)
#
# Same pairing as the original scripts: MRGN_sim2_filter_int_child.R uses the CS-q call and
# its _liberal.R counterpart the CS-alpha one, identical in every other argument.
#
# get.conf.trios() is called ONCE and both settings are derived from that one result. Every
# expensive step inside it -- propagate::bigcor() over the variant/covariate correlation and
# the per-trio-per-covariate regressions in p.from.reg() -- is identical for the two
# settings; only the final switch(adjust_by, ...) thresholding of reg.pvalues differs, and
# reg.pvalues is returned. Calling it twice therefore doubled the cost for nothing: measured
# at ~0.955 ms per (trio x covariate) test, a 300-trio group with an 8,340 column pool takes
# ~40 min per call.
#
# CS-alpha is reproduced exactly by thresholding reg.pvalues at alpha, which is all
# adjust_by = "none" does internally.
#
# trios:      list of n x 3 trio matrices (V1, T1, T2) that all share the same sample size
# cov.pool:   n x p matrix of candidate covariates pooled over those trios. The confounder
#             columns are suffixed with the dataset index by name.trio.columns(), which is
#             what keeps the pooled column names unique.
# known.conf: n x k matrix of known (clinical) confounders to include with every trio.
#             Only the n = 670 datasets have a K block, so this defaults to none.
#
# Each setting returns trios.with.confs: one data.frame per trio, assembled as
# trio + known confounders + selected confounders, ready to hand to apply.mrgn() or
# apply.mrpc().
#
# get.conf.trios() indexes cov.pool by row, so every trio in the call must have the same
# number of observations as the pool: the driver therefore groups the datasets by sample
# size and calls this once per group. It also drops dimensions and fails on a list of one
# trio, hence the guard below.
select.confounders <- function(trios, cov.pool, known.conf = NULL, blocksize = 500,
                               filter_int_child = TRUE, selection_fdr = 0.05,
                               filter_fdr = 0.1, alpha = 0.01) {

    if (length(trios) < 2) {
        stop("select.confounders() needs at least 2 trios: get.conf.trios() fails on a single trio")
    }
    if (nrow(cov.pool) != nrow(trios[[1]])) {
        stop("cov.pool must have the same number of rows as the trios (group the datasets by sample size)")
    }
    # get.conf.trios() blocks the variant/covariate correlation with propagate::bigcor(),
    # whose rep(1:NBLOCKS, each = size) produces a garbled split when blocksize exceeds the
    # number of trios (NBLOCKS = 0, and 1:0 is not empty). It fails downstream with an
    # opaque "'x' is empty", so catch the misconfiguration here instead.
    if (filter_int_child && blocksize > length(trios)) {
        stop(paste0("blocksize (", blocksize, ") must not exceed the number of trios (",
                    length(trios), "): propagate::bigcor() mis-blocks the correlation otherwise"))
    }

    # get.conf.trios() returns sig.asso.covs from apply(reg.sigmat, 1, which), which
    # simplifies to a matrix on the (unlikely but possible) occasion that every trio
    # selects the same number of covariates. Everything downstream indexes it with [[i]],
    # so normalise it to a list of integer vectors first.
    as.cov.list <- function(x, n.trios, nms) {
        if (!is.list(x)) {
            x <- if (is.matrix(x)) {
                lapply(seq_len(ncol(x)), function(i) as.integer(x[, i]))
            } else {
                lapply(seq_len(n.trios), function(i) as.integer(x[i]))
            }
        }
        x <- lapply(x, function(v) unname(as.integer(v)))
        names(x) <- nms
        return(x)
    }

    # trio + known confounders + selected confounders, one data.frame per trio, plus the
    # per-trio covariate frame the original scripts kept. Shared by both settings.
    assemble <- function(sig.asso.covs, setting, elapsed) {
        conf.list <- lapply(sig.asso.covs, function(x, y) { y[, x] }, y = cov.pool)

        # drop = FALSE keeps a single selected covariate a data.frame rather than a vector,
        # which would otherwise lose its column name in the cbind.
        trios.with.confs <- vector("list", length(trios))
        for (i in seq_along(trios)) {
            assembled <- as.data.frame(trios[[i]])
            if (!is.null(known.conf)) {
                assembled <- cbind.data.frame(assembled, as.data.frame(known.conf))
            }
            selected <- sig.asso.covs[[i]]
            if (length(selected) > 0) {
                assembled <- cbind.data.frame(assembled, cov.pool[, selected, drop = FALSE])
            }
            trios.with.confs[[i]] <- assembled
        }
        names(trios.with.confs) <- names(trios)

        list(trios.with.confs = trios.with.confs,
             conf.list = conf.list,
             sig.asso.covs = sig.asso.covs,
             setting = setting,
             # FALSE when the group fell back to the unfiltered selection below. Both
             # settings share the one call, so they now always agree; both are still
             # reported so the CSq./CSa. results columns stay symmetric.
             filter_int_child = filter.applied,
             time.seconds = elapsed)
    }

    filter.applied <- filter_int_child
    start.time <- Sys.time()
    # get.conf.trios() stops outright when filter_int_child = TRUE and not one covariate in
    # the pool comes back associated with any variant in the group: its apply(sigmat, 1,
    # which) simplifies the all-empty result to a zero-length vector and it errors with
    # "No common child or intermediate variables detected". There is nothing to filter in
    # that case, so fall back to the unfiltered selection and say so.
    #
    # adjust_by = "all" gives CS-q directly; CS-alpha is derived from reg.pvalues below.
    selection <- tryCatch(
        get.conf.trios(trios = trios,
                       cov.pool = cov.pool,
                       blocksize = blocksize,
                       filter_int_child = filter_int_child,
                       selection_fdr = selection_fdr,
                       filter_fdr = filter_fdr,
                       adjust_by = "all",
                       alpha = alpha),
        error = function(e) {
            if (!grepl("No common child or intermediate variables detected", conditionMessage(e))) {
                stop(e)
            }
            message("No common child or intermediate variables detected in this group: ",
                    "re-running get.conf.trios() with filter_int_child = FALSE")
            filter.applied <<- FALSE
            get.conf.trios(trios = trios,
                           cov.pool = cov.pool,
                           blocksize = blocksize,
                           filter_int_child = FALSE,
                           selection_fdr = selection_fdr,
                           filter_fdr = filter_fdr,
                           adjust_by = "all",
                           alpha = alpha)
        })
    elapsed <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))

    trio.names <- names(trios)
    csq.covs <- as.cov.list(selection$sig.asso.covs, length(trios), trio.names)

    # CS-alpha: threshold the same regression p-values at alpha, which is exactly what
    # adjust_by = "none" does inside get.conf.trios(). Row-wise lapply rather than
    # apply(..., 1, which), which would simplify to a matrix when every trio happens to
    # select the same number of covariates. which() drops NA on its own, matching the
    # existing behaviour for the rank-deficient fits that occur at small n.
    reg.p <- as.matrix(selection$reg.pvalues)
    csa.covs <- lapply(seq_len(nrow(reg.p)), function(i) unname(which(reg.p[i, ] < alpha)))
    names(csa.covs) <- trio.names

    CS.q <- assemble(csq.covs, setting = "CS-q", elapsed = elapsed)
    CS.alpha <- assemble(csa.covs, setting = "CS-alpha", elapsed = elapsed)
    CS.q$adjust_by <- "all"
    CS.alpha$adjust_by <- "none"

    # `selection` is stored ONCE at the top level rather than copied into both settings.
    # The two copies were identical apart from sig.asso.covs, which now lives on each
    # setting, and the duplicate cost ~230 MB per group -- get.conf.trios() returns five
    # numeric and two logical matrices of (trios x pool), ~8,400 columns wide here.
    # Read a setting's selected indices as sel$CS.q$sig.asso.covs.
    #
    # one call now serves both settings, so time.seconds is the wall clock for the pair
    return(list(selection = selection,
                CS.q = CS.q,
                CS.alpha = CS.alpha,
                time.seconds = elapsed))
}


# ---------------------------------------------------------------------------------------
# caching the selection
# ---------------------------------------------------------------------------------------
#
# Confounder selection is by far the most expensive step in the pipeline -- ~0.955 ms per
# (trio x covariate) test, so ~30 min for a 300-trio group against an ~8,400 column pool,
# ~2.5 h over the five sample-size groups -- and it does not depend on which method is
# fitted afterwards. Running it once and caching the result lets the inference script be
# re-run freely.
#
# The cached file holds the whole select.confounders() result plus a provenance block, so
# a stale cache can be detected rather than silently used. That matters here: the
# simulated data has been regenerated more than once during this revision, and a cache
# keyed only on sample size would survive a regeneration and quietly pair the wrong
# confounder indices with the wrong trios.
#
# Size: ~1 GB on disk across the five groups. Dominated by `selection`'s (trios x pool)
# matrices and by trios.with.confs.

selection.cache.file <- function(out.dir, size) {
    file.path(out.dir, paste0("selection_group_n", size, ".RData"))
}


# Everything the cache needs to prove it belongs to this request.
selection.provenance <- function(datasets, cov.names, settings, sim.file) {
    list(n = nrow(datasets[[1]]$data),
         datasets = sapply(datasets, function(x) x$params$dataset),
         cov.names = cov.names,
         settings = settings,
         sim.file = sim.file,
         created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}


# Returns NULL when the cache matches the request, or a string describing the first
# mismatch. Checked in the order a regeneration would break them.
selection.cache.mismatch <- function(cached, want) {
    if (is.null(cached$sel) || is.null(cached$datasets)) {
        return("cache predates the provenance block")
    }
    if (!identical(cached$n, want$n)) {
        return(paste0("sample size ", cached$n, " != ", want$n))
    }
    if (!identical(cached$datasets, want$datasets)) {
        return(paste0("covers ", length(cached$datasets), " datasets, request covers ",
                      length(want$datasets),
                      if (length(cached$datasets) == length(want$datasets))
                          " (same count, different dataset indices)" else ""))
    }
    if (!identical(cached$cov.names, want$cov.names)) {
        return("covariate pool column names differ (simulated data regenerated?)")
    }
    changed <- names(want$settings)[!mapply(identical, cached$settings[names(want$settings)],
                                            want$settings)]
    if (length(changed) > 0) {
        return(paste("selection settings changed:", paste(changed, collapse = ", ")))
    }
    return(NULL)
}


# Load the cached selection for one sample-size group, or compute and cache it.
#
# datasets: the simulated datasets for ONE sample-size group
# rerun:    TRUE recomputes and overwrites even when a valid cache exists
# save:     FALSE computes without writing. Pass this when running on a subset of a
#           group (max.per.group), so a smoke test cannot overwrite a good full cache.
get.selection <- function(datasets, out.dir, sim.file = NA_character_, rerun = FALSE,
                          save = TRUE, verbose = TRUE, selection_fdr = 0.05,
                          filter_fdr = 0.1, alpha = 0.01, filter_int_child = TRUE) {

    size <- nrow(datasets[[1]]$data)
    path <- selection.cache.file(out.dir, size)
    cov.pool <- group.cov.pool(datasets)
    known.conf <- group.known.conf(datasets)
    blocksize <- min(500, length(datasets))

    settings <- list(selection_fdr = selection_fdr, filter_fdr = filter_fdr,
                     alpha = alpha, filter_int_child = filter_int_child,
                     blocksize = blocksize)
    want <- selection.provenance(datasets, colnames(cov.pool), settings, sim.file)

    if (!rerun && file.exists(path)) {
        cached <- loadRData(path)
        bad <- selection.cache.mismatch(cached, want)
        if (is.null(bad)) {
            if (verbose) {
                cat("  loaded cached selection:", basename(path),
                    "(built", cached$created, ")\n")
            }
            return(cached$sel)
        }
        message("  cached selection at ", basename(path), " does not match this request (",
                bad, "); recomputing")
    }

    if (verbose) {
        cat("  selecting confounders:", length(datasets), "trios,", ncol(cov.pool),
            "pooled covariates -- this is the slow step\n")
    }
    sel <- select.confounders(lapply(datasets, function(x) x$data[, 1:3]),
                              cov.pool, known.conf = known.conf, blocksize = blocksize,
                              filter_int_child = filter_int_child,
                              selection_fdr = selection_fdr, filter_fdr = filter_fdr,
                              alpha = alpha)

    if (save) {
        cache <- c(list(sel = sel), want)
        base::save(cache, file = path)
        if (verbose) {
            cat("  cached ->", basename(path), "|",
                round(file.size(path)/1024^2), "MB |",
                round(sel$time.seconds/60, 1), "min\n")
        }
    } else if (verbose) {
        cat("  not caching (save = FALSE)\n")
    }
    return(sel)
}


# ---------------------------------------------------------------------------------------
# MRGN
# ---------------------------------------------------------------------------------------

# a function to bootstrap edge probabilities for a single trio
#
# Each replicate keeps the first six elements of the infer.trio() result
# (b11, b12, b21, b22, V1:T1, V1:T2) rather than just the adjacency. The two interaction
# terms are what MRGN::class.vec() uses to separate M2 from M4, and they are needed
# exactly when all four edge indicators are 1, so an averaged adjacency alone cannot name
# the bootstrap-supported model.
#
# cl: optional cluster from parallel::makeCluster(). Workers need MRGN attached. The
#     serial path is kept for cl = NULL.
boostrap_edge_probabilities <- function(trio, number_of_samples=1000, cl = NULL) {

    one.replicate <- function(i, trio) {
        sampled_trio <- trio[sample(nrow(trio), replace=TRUE), ]
        result <- MRGN::infer.trio(
            trio = sampled_trio,
            use.perm = FALSE
        )
        unlist(result[1:6])
    }

    if (is.null(cl)) {
        indicators <- lapply(1:number_of_samples, one.replicate, trio = trio)
    } else {
        indicators <- parallel::parLapply(cl, 1:number_of_samples, one.replicate, trio = trio)
    }

    # proportion of replicates in which each indicator was called significant
    indicator.means <- Reduce("+", indicators) / number_of_samples
    edge.probabilities <- MRGN::get.adj(as.list(indicator.means))

    # the model the bootstrap supports: majority vote on each indicator, then MRGN's own
    # classification of the resulting vector
    supported <- as.numeric(indicator.means >= 0.5)
    boot.model <- MRGN::class.vec(supported)
    present <- indicator.means[supported == 1]

    return(list(edge.probabilities = edge.probabilities,
                indicator.means = indicator.means,
                boot.model = as.character(boot.model),
                min.edge.prob = if (length(present) > 0) min(present) else NA_real_,
                number_of_samples = number_of_samples))
}

# This function applies the MRGN method to a single dataset (trio plus its selected
# confounders). time.seconds covers the infer.trio() call only; the bootstrap is opt-in
# and timed separately so the reported inference time stays comparable to MRPC and GMAC.
apply.mrgn <- function(data, use.perm = FALSE, bootstrap = FALSE, number_of_samples = 1000,
                       cl = NULL, verbose = TRUE) {

    start.time <- Sys.time()
    result <- MRGN::infer.trio(
        trio = data,
        use.perm = use.perm,
        verbose = verbose
    )
    end.time <- Sys.time()

    # by name, not by position: infer.trio() returns Inferred.Model as element 18, so the
    # result[[14]] used by the original scripts now picks up coef11 instead
    inferred_model <- result$Inferred.Model
    inferred_adj <- MRGN::get.adj(result)

    boot <- NULL
    bootstrap.time.seconds <- NA_real_
    if (bootstrap) {
        boot.start <- Sys.time()
        boot <- boostrap_edge_probabilities(data, number_of_samples = number_of_samples, cl = cl)
        boot.end <- Sys.time()
        bootstrap.time.seconds <- as.numeric(difftime(boot.end, boot.start, units = "secs"))
    }

    final_result <- list(model = inferred_model,
                         adj = inferred_adj,
                         edge.probabilities = boot$edge.probabilities,
                         bootstrap = boot,
                         details = result,
                         time.seconds = as.numeric(difftime(end.time, start.time, units = "secs")),
                         bootstrap.time.seconds = bootstrap.time.seconds)
    return(final_result)
}


# ---------------------------------------------------------------------------------------
# MRPC
# ---------------------------------------------------------------------------------------

# The eight trio adjacency matrices MRPC output is compared against, built once and cached:
# this runs on every trio, and rebuilding them each call is wasted work.
mrpc.truth.adjacencies <- local({
    cached <- NULL
    function() {
        if (!is.null(cached)) {
            return(cached)
        }
        MRPCtruth <- NULL # assigned by data(); declared here to keep the lookup explicit
        utils::data("MRPCtruth", package = "MRPC", envir = environment())

        #V1-->T1
        Adj.M0 <- as(MRPCtruth$M0, "matrix")
        #V1-->T2
        Adj.M01 <- matrix(0, nrow = 3, ncol = 3)
        rownames(Adj.M01) <- colnames(Adj.M01) <- colnames(Adj.M0)
        Adj.M01[1, 3] <- 1

        #V1-->T1-->T2
        Adj.M1 <- as(MRPCtruth$M1, "matrix")
        #V1-->T2-->T1
        Adj.M11 <- matrix(0, nrow = 3, ncol = 3)
        rownames(Adj.M11) <- colnames(Adj.M11) <- colnames(Adj.M1)
        Adj.M11[1, 3] <- 1
        Adj.M11[3, 2] <- 1

        #V1-->T1<--T2
        Adj.M2 <- as(MRPCtruth$M2, "matrix")
        #V1-->T2<--T1
        Adj.M21 <- matrix(0, nrow = 3, ncol = 3)
        rownames(Adj.M21) <- colnames(Adj.M21) <- colnames(Adj.M2)
        Adj.M21[1, 3] <- 1
        Adj.M21[2, 3] <- 1

        #V1-->T1, V1-->T2
        Adj.M3 <- as(MRPCtruth$M3, "matrix")
        #V1-->T1, V1-->T2, T1--T2
        Adj.M4 <- as(MRPCtruth$M4, "matrix")

        cached <<- list(M0.1 = Adj.M0, M0.2 = Adj.M01,
                        M1.1 = Adj.M1, M1.2 = Adj.M11,
                        M2.1 = Adj.M2, M2.2 = Adj.M21,
                        M3 = Adj.M3, M4 = Adj.M4)
        return(cached)
    }
})


# classify an inferred 3 x 3 adjacency matrix against the eight trio models
classify.mrpc.adj <- function(adj) {
    truth <- mrpc.truth.adjacencies()
    colnames(adj) <- rownames(adj) <- colnames(truth$M0.2)
    model <- "Other"
    for (nm in names(truth)) {
        if (identical(truth[[nm]], adj)) {
            model <- nm
        }
    }
    return(list(model = model, adj = adj))
}


# This function applies the MRPC method to a single dataset (trio plus its selected
# confounders).
#
# MRPC can take hours on trios with many confounders, so the fit is capped with
# R.utils::withTimeout(). Handling the cap inside the function means a timed-out trio
# still returns a complete record (model = NA, timed.out = TRUE) instead of being lost:
# assigning the fallback from inside a tryCatch handler in a driver loop writes to the
# handler's own frame, not the caller's, and silently leaves a NULL in the results list.
# withTimeout() interrupts between R evaluation steps, which is what MRPC spends its time
# in; it cannot break into long-running compiled code.
apply.mrpc <- function(data, timeout = 120, GV = 1, FDR = 0.05, alpha = 0.01,
                       indepTest = "gaussCItest", FDRcontrol = "ADDIS", verbose = FALSE) {

    X <- as.matrix(data)
    n <- nrow(X)
    V <- colnames(X)     # Column names
    # Classical correlation
    suffStat <- list(C = cor(X, use = "complete.obs"),
                     n = n)

    timed.out <- FALSE
    error.message <- NA_character_
    start.time <- Sys.time()
    MRPC.fit.FDR <- tryCatch(
        R.utils::withTimeout({
            MRPC(X,
                 suffStat,
                 GV = GV,
                 FDR = FDR,
                 alpha = alpha,
                 indepTest = indepTest,
                 labels = V,
                 FDRcontrol = FDRcontrol,
                 verbose = verbose)
        }, timeout = timeout, onTimeout = "error"),
        error = function(e) {
            timed.out <<- inherits(e, "TimeoutException")
            error.message <<- conditionMessage(e)
            return(NULL)
        })
    end.time <- Sys.time()

    if (is.null(MRPC.fit.FDR)) {
        # timed out, or MRPC errored on this trio: record it and let the driver continue
        return(list(model = NA_character_,
                    adj = NA,
                    details = NULL,
                    time.seconds = NA_real_,
                    timed.out = timed.out,
                    timeout.seconds = timeout,
                    error.message = error.message))
    }

    Adj.infe1 <- as(MRPC.fit.FDR@graph, "matrix")
    Adj.infe <- Adj.infe1[1:3, 1:3] #only consider snp, cis, trans
    classified <- classify.mrpc.adj(Adj.infe)

    return(list(model = classified$model,
                adj = classified$adj,
                details = MRPC.fit.FDR,
                time.seconds = as.numeric(difftime(end.time, start.time, units = "secs")),
                timed.out = FALSE,
                timeout.seconds = timeout,
                error.message = NA_character_))
}


# ---------------------------------------------------------------------------------------
# GMAC
# ---------------------------------------------------------------------------------------

# GMAC tests mediation in each direction, so a trio is labelled by which of the two tests
# is significant. Same mapping as replace.tf() in GTEx/scripts/make_mrgn_triotables.R.
# GMAC's own vocabulary is kept rather than the M labels: it only tests mediation and
# cannot separate M0 from M2/M3.
gmac.model.call <- function(p.cis, p.trans, alpha = 0.05) {
    cis.sig <- p.cis < alpha
    trans.sig <- p.trans < alpha
    model <- rep(NA_character_, length(cis.sig))
    model[cis.sig & trans.sig] <- "Undirected"
    model[cis.sig & !trans.sig] <- "Cis Mediated"
    model[!cis.sig & trans.sig] <- "Trans Mediated"
    model[!cis.sig & !trans.sig] <- "No Mediation"
    return(model)
}


# This function times GMAC inference on a single trio with gmacOneTrio(). The results
# reported for the manuscript come from run.gmac.all() / gmac(); this path exists to
# measure the per-trio cost of that same work, so it repeats the same test rather than
# replacing it.
#
# The confounder set is the one GMAC_many_conf.R reports on, Known_sel_pool: the known
# (clinical) confounders plus the pool covariates GMAC selected for this trio, which come
# from run.gmac.all()'s sel.conf.ind. Either part may be empty -- the K block only exists
# for the n = 670 datasets, and a trio can end up with no pool covariate selected -- in
# which case the test is simply run with whatever remains.
#
# trio: n x 3 matrix ordered SNP, cis gene, trans gene. The cis-mediator test uses that
# order; the trans-mediator test swaps columns 2 and 3.
apply.gmac <- function(trio, confounders = NULL, known.conf = NULL, nperm = 1000,
                       nominal.p = TRUE, alpha = 0.05) {

    trio <- as.matrix(trio)
    empty <- matrix(nrow = nrow(trio), ncol = 0)
    confounders <- if (is.null(confounders)) empty else as.matrix(confounders)
    known.conf <- if (is.null(known.conf)) empty else as.matrix(known.conf)
    known.sel.pool <- cbind(known.conf, confounders)

    start.time <- Sys.time()
    cis <- gmacOneTrio(trio = trio, confounders = known.sel.pool,
                       nperm = nperm, nominal.p = nominal.p)
    trans <- gmacOneTrio(trio = trio[, c(1, 3, 2)], confounders = known.sel.pool,
                         nperm = nperm, nominal.p = nominal.p)
    end.time <- Sys.time()

    results <- data.frame(inferred_model = gmac.model.call(cis$pval, trans$pval, alpha = alpha),
                          cispval = cis$pval,
                          transpval = trans$pval,
                          ciseffect = cis$beta_change,
                          transeffect = trans$beta_change,
                          stringsAsFactors = FALSE)

    return(list(results = results,
                n.confounders = ncol(known.sel.pool),
                alpha = alpha,
                time.seconds = as.numeric(difftime(end.time, start.time, units = "secs"))))
}


# This function applies GMAC to a whole group of trios at once, the way
# Simulation/results/GMAC/GMAC_many_conf.R does, and is where GMAC performs its own
# confounder selection: gmac() screens the covariate pool with stratified FDR across the
# trios in the call, which is why this cannot be done one trio at a time.
#
# trios:      list of n x 3 trio matrices sharing a sample size
# cov.pool:   n x p candidate covariate matrix pooled over those trios
# known.conf: n x k matrix of known (clinical) confounders. Only the n = 670 datasets have
#             a K block, so this defaults to no known confounders.
#
# Returns the cis and trans output tables (p-values, effect changes, model call), the
# per-trio selected confounder indicators, and the wall time of the two gmac() calls.
run.gmac.all <- function(trios, cov.pool, known.conf = NULL, nperm = 1000, nominal.p = TRUE,
                         fdr = 0.05, fdr_filter = 0.1, alpha = 0.05, cl = NULL) {

    num.trios <- length(trios)
    sample.size <- nrow(trios[[1]])
    if (nrow(cov.pool) != sample.size) {
        stop("cov.pool must have the same number of rows as the trios (group the datasets by sample size)")
    }
    if (is.null(known.conf)) {
        known.conf <- matrix(nrow = sample.size, ncol = 0)
    }

    # assemble the trios into the gene x sample layout gmac() expects
    trios.mat <- do.call("cbind", lapply(trios, as.matrix))
    snp.idx <- seq(1, ncol(trios.mat) - 2, 3)
    snp.dat.cis <- trios.mat[, snp.idx]
    exp.dat <- trios.mat[, -snp.idx]
    trios.idx <- cbind(c(1:num.trios),
                       matrix(c(1:ncol(exp.dat)), nrow = num.trios, ncol = 2, byrow = TRUE))

    # gmac() transposes known.conf, cov.pool, exp.dat and snp.dat.cis internally, so they
    # go in as covariates/genes by samples
    input.list <- list(kc = t(as.matrix(known.conf)),
                       cov.pool = t(as.matrix(cov.pool)),
                       snp.dat.cis = t(snp.dat.cis),
                       exp.dat = t(exp.dat),
                       trio.indexes = trios.idx,
                       nperm = nperm,
                       use.nominal.p = nominal.p,
                       fdr = fdr,
                       fdr_filter = fdr_filter)

    start.time <- Sys.time()
    print("applying gmac with cis mediators")
    output.cis.med <- gmac(cl = cl, known.conf = input.list$kc, cov.pool = input.list$cov.pool,
                           exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                           trios.idx = input.list$trio.indexes, nperm = input.list$nperm,
                           nominal.p = input.list$use.nominal.p, fdr = input.list$fdr,
                           fdr_filter = input.list$fdr_filter)

    print("applying gmac with trans mediators")
    # columns 2 and 3 swapped: the trans gene becomes the mediator
    output.trans.med <- gmac(cl = cl, known.conf = input.list$kc, cov.pool = input.list$cov.pool,
                             exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
                             trios.idx = input.list$trio.indexes[, c(1, 3, 2)], nperm = input.list$nperm,
                             nominal.p = input.list$use.nominal.p, fdr = input.list$fdr,
                             fdr_filter = input.list$fdr_filter)
    end.time <- Sys.time()

    #reorganize output for saving
    #cis results
    out.table.cis <- cbind.data.frame(output.cis.med[[1]], output.cis.med[[2]])
    colnames(out.table.cis) <- c(paste0("pval_", colnames(output.cis.med[[1]])),
                                 paste0("effect_change_", colnames(output.cis.med[[2]])))

    #trans results
    out.table.trans <- cbind.data.frame(output.trans.med[[1]], output.trans.med[[2]])
    colnames(out.table.trans) <- c(paste0("pval_", colnames(output.trans.med[[1]])),
                                   paste0("effect_change_", colnames(output.trans.med[[2]])))

    # one row per trio, from the Known_sel_pool columns: the mediation test run with the
    # known confounders plus the covariates GMAC selected for that trio, which is what
    # GMAC_many_conf.R reports on
    results <- data.frame(inferred_model = gmac.model.call(out.table.cis$pval_Known_sel_pool,
                                                           out.table.trans$pval_Known_sel_pool,
                                                           alpha = alpha),
                          cispval = out.table.cis$pval_Known_sel_pool,
                          transpval = out.table.trans$pval_Known_sel_pool,
                          ciseffect = out.table.cis$effect_change_Known_sel_pool,
                          transeffect = out.table.trans$effect_change_Known_sel_pool,
                          stringsAsFactors = FALSE)

    return(list(results = results,
                cis = out.table.cis,
                trans = out.table.trans,
                alpha = alpha,
                sel.conf.ind = output.cis.med$sel.conf.ind,
                sel.conf.ind.trans = output.trans.med$sel.conf.ind,
                input.list = input.list,
                time.seconds = as.numeric(difftime(end.time, start.time, units = "secs"))))
}
