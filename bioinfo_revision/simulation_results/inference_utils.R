# Inference utilities for the revised trio simulation.
#
# Sourced by apply_mrgn.R, apply_mrpc.R, apply_gmac.R, run_all_inference.R and
# run_confounder_selection.R, alongside inference_config.R which holds the settings.
#
# Each apply.* function runs one method on a SINGLE trio and attaches the time taken by
# the inference itself, in seconds, so the three methods can be compared directly. Confounder selection is a separate step: MRGN
# and MRPC share MRGN::get.conf.trios() via select.confounders(), while GMAC performs its
# own selection inside run.gmac.all().
#
# Paths below are relative to the repository root, as everywhere in this stage.
# gmac_one_trio.R and GMAC_moded/R/GMAC.R both define Indirect(), my.solve(),
# nominal.pfun() and get.beta.change() with identical bodies, so sourcing both is
# deliberate and harmless -- gmacOneTrio() lives in the first, gmac() in the second.
source("adapted_GMAC_func/gmac_one_trio.R")
source("adapted_GMAC_func/GMAC_moded/R/GMAC.R")


# ---------------------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------------------
#
# Shared by the apply_*.R scripts and by run_confounder_selection.R.

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
# result rows
# ---------------------------------------------------------------------------------------

# What goes in the error column. No new columns are needed for the bootstrap: a row is a
# failure if and only if model is NA, so a note on a row that HAS a model is informational
# rather than fatal. run.mrpc.group() already uses the field this way for a disabled arm.
mrgn.error.note <- function(res, error) {
    if (!is.na(error)) return(error)                       # infer.trio() itself failed
    if (!is.null(res$bootstrap.error) && !is.na(res$bootstrap.error)) {
        return(paste("bootstrap failed, model retained:", res$bootstrap.error))
    }
    n <- res$bootstrap$n.dropped
    if (!is.null(n) && !is.na(n) && n > 0) {
        return(sprintf("bootstrap: %d of %d resamples dropped for no genotype variation",
                       n, res$bootstrap$number_of_samples))
    }
    NA_character_
}

# A NULL res means infer.trio() itself failed and there is no inferred model. That is now
# the only way model comes back NA -- a bootstrap that fails leaves the model intact and
# explains itself in the error column. See apply.mrgn().
mrgn.fields <- function(res, error, truth.label) {
    if (is.null(res)) {
        return(list(model = NA_character_, correct = NA, correct.coarse = NA,
                    time.seconds = NA_real_, boot.model = NA_character_,
                    boot.min.edge.prob = NA_real_, boot.p.V1T1 = NA_real_,
                    boot.p.T1T2 = NA_real_, boot.p.V1T2 = NA_real_,
                    boot.p.T2T1 = NA_real_, bootstrap.time.seconds = NA_real_,
                    error = error))
    }
    ind <- res$bootstrap$indicator.means
    list(model                  = as.character(res$model),
         correct                = res$model == truth.label,
         correct.coarse         = coarse.model(res$model) == coarse.model(truth.label),
         time.seconds           = res$time.seconds,
         boot.model             = if (is.null(res$bootstrap)) NA_character_ else res$bootstrap$boot.model,
         boot.min.edge.prob     = if (is.null(res$bootstrap)) NA_real_ else res$bootstrap$min.edge.prob,
         # b11 = V1->T1, b12 = T1->T2, b21 = V1->T2, b22 = T2->T1
         boot.p.V1T1            = if (is.null(ind)) NA_real_ else unname(ind[1]),
         boot.p.T1T2            = if (is.null(ind)) NA_real_ else unname(ind[2]),
         boot.p.V1T2            = if (is.null(ind)) NA_real_ else unname(ind[3]),
         boot.p.T2T1            = if (is.null(ind)) NA_real_ else unname(ind[4]),
         bootstrap.time.seconds = res$bootstrap.time.seconds,
         error                  = mrgn.error.note(res, error))
}

mrpc.fields <- function(res, error, truth.label) {
    if (is.null(res)) {
        return(list(model = NA_character_, correct = NA, correct.coarse = NA,
                    time.seconds = NA_real_, timed.out = NA, error = error))
    }
    list(model          = as.character(res$model),
         correct        = !is.na(res$model) && res$model == truth.label,
         correct.coarse = !is.na(res$model) && coarse.model(res$model) == coarse.model(truth.label),
         time.seconds   = res$time.seconds,
         timed.out      = res$timed.out,
         error          = if (is.na(error)) res$error.message else error)
}

# GMAC reports a mediation call, not a model label, so no correctness flag is derived here:
# comparing "Cis Mediated" against M0/M1/M2 is a cross-tab decision for the analysis stage.
gmac.fields <- function(batch.row, timing, error) {
    list(model        = if (is.null(batch.row)) NA_character_ else as.character(batch.row$inferred_model),
         cispval      = if (is.null(batch.row)) NA_real_ else batch.row$cispval,
         transpval    = if (is.null(batch.row)) NA_real_ else batch.row$transpval,
         ciseffect    = if (is.null(batch.row)) NA_real_ else batch.row$ciseffect,
         transeffect  = if (is.null(batch.row)) NA_real_ else batch.row$transeffect,
         time.seconds = if (is.null(timing)) NA_real_ else timing$time.seconds,
         error        = error)
}


# ---------------------------------------------------------------------------------------
# one sample-size group
# ---------------------------------------------------------------------------------------

# The columns every per-method file carries, so each is self-contained: which trio, the
# parameters that generated it, and the truth. The selection scores are added on top by
# the MRGN and MRPC builders (which used the CS sets) and by the GMAC builder (which used
# its own), and de-duplicated when the master is joined.
id.columns <- function(dat, params, index) {
    true.confs <- own.block.names(colnames(dat), "U", index)
    c(list(dataset = index),
      as.list(params[, setdiff(colnames(params), "dataset")]),
      list(truth.model = unname(TRUTH.LABEL[params$model]),
           n.true.confs = length(true.confs),
           true.confs = paste(true.confs, collapse = ";")))
}

# CS-q / CS-alpha scores for one trio, shared by the MRGN and MRPC files
cs.score.columns <- function(sel, i, index, dat, cov.names) {
    true.confs <- own.block.names(colnames(dat), "U", index)
    c(prefixed(score.selection(cov.names[sel$CS.q$sig.asso.covs[[i]]], index, true.confs), "CSq"),
      prefixed(score.selection(cov.names[sel$CS.alpha$sig.asso.covs[[i]]], index, true.confs), "CSa"),
      list(CSq.filter_int_child = sel$CS.q$filter_int_child,
           CSa.filter_int_child = sel$CS.alpha$filter_int_child))
}


# ---- MRGN: true confounders, CS-q and CS-alpha ----
run.mrgn.group <- function(datasets, sel, cov.names, cl = NULL, bootstrap = TRUE,
                           number_of_samples = n.bootstrap, verbose = TRUE) {
    n.datasets <- length(datasets)
    rows <- vector("list", n.datasets)
    for (i in seq_len(n.datasets)) {
        dat <- datasets[[i]]$data; params <- datasets[[i]]$params
        index <- params$dataset; K_n <- params$K_n
        truth.label <- unname(TRUTH.LABEL[params$model])

        mrgn.truth <- safely(apply.mrgn(ground.truth.input(dat, K_n, index),
                                        bootstrap = bootstrap,
                                        number_of_samples = number_of_samples,
                                        cl = cl, verbose = FALSE))
        mrgn.csq <- safely(apply.mrgn(sel$CS.q$trios.with.confs[[i]], bootstrap = bootstrap,
                                      number_of_samples = number_of_samples,
                                      cl = cl, verbose = FALSE))
        mrgn.csa <- safely(apply.mrgn(sel$CS.alpha$trios.with.confs[[i]], bootstrap = bootstrap,
                                      number_of_samples = number_of_samples,
                                      cl = cl, verbose = FALSE))

        rows[[i]] <- as.data.frame(c(
            id.columns(dat, params, index),
            cs.score.columns(sel, i, index, dat, cov.names),
            prefixed(mrgn.fields(mrgn.truth$value, mrgn.truth$error, truth.label), "mrgn.truth"),
            prefixed(mrgn.fields(mrgn.csq$value, mrgn.csq$error, truth.label), "mrgn.CSq"),
            prefixed(mrgn.fields(mrgn.csa$value, mrgn.csa$error, truth.label), "mrgn.CSa")),
            stringsAsFactors = FALSE)

        if (verbose && (i %% 25 == 0 || i == n.datasets)) cat("    mrgn", i, "/", n.datasets, "\n")
    }
    do.call(rbind, rows)
}


# ---- MRPC: true confounders, CS-q and CS-alpha ----
#
# The truth arm is the same oracle MRGN gets: trio + K + that trio's own U block, via
# ground.truth.input(). It is what makes the MRPC and MRGN columns readable against each
# other -- the gap between the oracle and CS-q is the cost of selection rather than of the
# method, which is the comparison METHODS.md section 5 is built on.
#
# It is also the arm most likely to time out. At n = 670 it carries a median of 25-29
# confounders against CS-q's 16, and CS-q already times out on 61% of trios there. See the
# mrpc.arms note in inference_config.R for how to smoke-test it before a full run.
run.mrpc.group <- function(datasets, sel, cov.names, timeout = mrpc.timeout,
                           arms = mrpc.arms, truth.max.n = mrpc.truth.max.n,
                           verbose = TRUE) {
    n.datasets <- length(datasets)
    rows <- vector("list", n.datasets)
    sample.size <- nrow(datasets[[1]]$data)

    # An arm that is switched off still produces its columns, so every mrpc_group_*.RData
    # has the same schema and groups run under different settings still rbind. The reason
    # goes in the error field, so "not attempted" is distinguishable from "timed out".
    skipped <- function(arm, why) {
        list(value = NULL, error = paste0(arm, " arm not attempted: ", why))
    }

    # Budget control: the truth arm carries more confounders than CS-q, which already times
    # out on most trios at the large sizes, so a large group would spend
    # n.trios x mrpc.timeout to write a mostly-NA column. See mrpc.truth.max.n in
    # inference_config.R, and set it to Inf to run the arm everywhere.
    do.truth <- "truth" %in% arms && sample.size <= truth.max.n
    if ("truth" %in% arms && !do.truth && verbose) {
        cat("    truth arm skipped at n =", sample.size, "( > mrpc.truth.max.n =",
            truth.max.n, ")\n")
    }

    for (i in seq_len(n.datasets)) {
        dat <- datasets[[i]]$data; params <- datasets[[i]]$params
        index <- params$dataset; K_n <- params$K_n
        truth.label <- unname(TRUTH.LABEL[params$model])

        mrpc.truth <- if (do.truth) {
            safely(apply.mrpc(ground.truth.input(dat, K_n, index), timeout = timeout))
        } else if ("truth" %in% arms) {
            skipped("truth", paste0("n = ", sample.size, " exceeds mrpc.truth.max.n = ",
                                    truth.max.n, "; the oracle confounder set is too large ",
                                    "for MRPC to fit at this sample size within the time ",
                                    "budget (see inference_config.R)"))
        } else {
            skipped("truth", "not listed in mrpc.arms")
        }
        mrpc.csq <- if ("CSq" %in% arms) {
            safely(apply.mrpc(sel$CS.q$trios.with.confs[[i]], timeout = timeout))
        } else skipped("CS-q", "not listed in mrpc.arms")
        mrpc.csa <- if ("CSa" %in% arms) {
            safely(apply.mrpc(sel$CS.alpha$trios.with.confs[[i]], timeout = timeout))
        } else skipped("CS-alpha", paste0("not listed in mrpc.arms: the CS-alpha set is ",
                                          "too large for MRPC to fit (see inference_config.R)"))

        rows[[i]] <- as.data.frame(c(
            id.columns(dat, params, index),
            cs.score.columns(sel, i, index, dat, cov.names),
            prefixed(mrpc.fields(mrpc.truth$value, mrpc.truth$error, truth.label), "mrpc.truth"),
            prefixed(mrpc.fields(mrpc.csq$value, mrpc.csq$error, truth.label), "mrpc.CSq"),
            prefixed(mrpc.fields(mrpc.csa$value, mrpc.csa$error, truth.label), "mrpc.CSa")),
            stringsAsFactors = FALSE)

        if (verbose && (i %% 25 == 0 || i == n.datasets)) cat("    mrpc", i, "/", n.datasets, "\n")
    }
    do.call(rbind, rows)
}


# ---- GMAC: selects its own confounders across the group, then one call per trio ----
# The batch run does the selection and the mediation test together, so unlike CS-q/CS-alpha
# it cannot be split out into the selection stage.
#
# TWO ARMS per trio, as for MRGN and MRPC:
#   gmac.*        GMAC's own selected confounders, scored from the batch table
#   gmac.truth.*  that trio's true U block, scored from its own apply.gmac() call
# The oracle arm costs one extra permutation test per trio (gmac.nperm each), which roughly
# doubles the per-trio part of the group. The batch gmac() call is unaffected -- it is
# GMAC's selection, and the oracle arm bypasses selection entirely.
run.gmac.group <- function(datasets, cov.pool, known.conf, cl = NULL,
                           nperm = gmac.nperm, verbose = TRUE) {
    n.datasets <- length(datasets)
    just.trios <- lapply(datasets, function(x) x$data[, 1:3])
    cov.names <- colnames(cov.pool)

    if (verbose) cat("    running GMAC over the group...\n")
    gmac.batch <- safely(run.gmac.all(just.trios, cov.pool, known.conf = known.conf,
                                      nperm = nperm, alpha = selection.alpha, cl = cl))
    if (!is.null(gmac.batch$error) && !is.na(gmac.batch$error)) {
        warning("GMAC failed for this group: ", gmac.batch$error)
    }
    gmac.out <- gmac.batch$value
    # When the batch call fails, every trio's gmac.* result is NA and gmac.error stays NA
    # too, because that column records the per-trio apply.gmac() call, which can succeed
    # independently. Carry the batch error into every row so a file full of NA models says
    # why. Seen in practice on small groups: qvalue cannot estimate pi0 from few p-values.
    batch.error <- if (is.null(gmac.batch$error)) NA_character_ else gmac.batch$error

    rows <- vector("list", n.datasets)
    for (i in seq_len(n.datasets)) {
        dat <- datasets[[i]]$data; params <- datasets[[i]]$params
        index <- params$dataset
        true.confs <- own.block.names(colnames(dat), "U", index)

        gmac.names <- if (is.null(gmac.out)) character(0) else
            cov.names[which(gmac.out$sel.conf.ind[i, ] == 1)]
        gmac.confs <- if (length(gmac.names) > 0)
            as.matrix(cov.pool[, gmac.names, drop = FALSE]) else NULL
        gmac.one <- safely(apply.gmac(just.trios[[i]], confounders = gmac.confs,
                                      known.conf = known.conf, nperm = nperm,
                                      alpha = selection.alpha))

        # Oracle arm: the same mediation test against that trio's TRUE confounders instead
        # of the ones GMAC selected, so the gap between the two is the cost of GMAC's own
        # selection rather than of the test. This one is scored from its OWN model call --
        # apply.gmac() returns a results row of the same shape as the batch table -- rather
        # than from gmac.out, which only ever describes the selected-confounder run.
        truth.confs <- if (length(true.confs) > 0)
            as.matrix(dat[, true.confs, drop = FALSE]) else NULL
        gmac.truth <- safely(apply.gmac(just.trios[[i]], confounders = truth.confs,
                                        known.conf = known.conf, nperm = nperm,
                                        alpha = selection.alpha))

        rows[[i]] <- as.data.frame(c(
            id.columns(dat, params, index),
            prefixed(score.selection(gmac.names, index, true.confs), "gmac"),
            prefixed(gmac.fields(if (is.null(gmac.out)) NULL else gmac.out$results[i, ],
                                 gmac.one$value, gmac.one$error), "gmac"),
            prefixed(gmac.fields(if (is.null(gmac.truth$value)) NULL else gmac.truth$value$results,
                                 gmac.truth$value, gmac.truth$error), "gmac.truth"),
            list(gmac.batch.error = batch.error)),
            stringsAsFactors = FALSE)

        if (verbose && (i %% 25 == 0 || i == n.datasets)) cat("    gmac", i, "/", n.datasets, "\n")
    }
    results <- do.call(rbind, rows)
    attr(results, "gmac.time.seconds") <-
        if (is.null(gmac.out)) NA_real_ else gmac.out$time.seconds
    results
}




# ---------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------


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
        # A resample can lose every copy of the minor allele -- at n = 50 and MAF 0.01 the
        # three carriers are all missed about 5% of the time -- and a constant V1 breaks
        # this in two different ways. Either the aliased genotype coefficient leaves a test
        # statistic NA and class.vec() dies branching on it, or infer.trio() returns
        # normally but with NA in the two correlation-based indicators, because cor() on a
        # zero-variance column is NA.
        #
        # Both are the same useless draw: a resample with no genotype variation says
        # nothing about the V1 edges. Both are dropped here, so neither can abort the
        # bootstrap (the first case) nor poison every indicator mean through Reduce (the
        # second). Returning NULL also keeps a worker-side error from surfacing as
        # parLapply's "N nodes produced errors".
        result <- tryCatch(MRGN::infer.trio(trio = sampled_trio, use.perm = FALSE),
                           error = function(e) NULL)
        if (is.null(result)) return(NULL)
        stats <- unlist(result[1:6])
        if (any(!is.finite(stats))) NULL else stats
    }

    if (is.null(cl)) {
        indicators <- lapply(1:number_of_samples, one.replicate, trio = trio)
    } else {
        indicators <- parallel::parLapply(cl, 1:number_of_samples, one.replicate, trio = trio)
    }

    unusable <- vapply(indicators, is.null, logical(1))
    n.dropped <- sum(unusable)
    indicators <- indicators[!unusable]
    if (length(indicators) == 0) {
        stop("every one of the ", number_of_samples, " bootstrap resamples failed")
    }

    # proportion of *usable* replicates in which each indicator was called significant, so
    # dropped resamples do not silently deflate every probability toward zero
    indicator.means <- Reduce("+", indicators) / length(indicators)
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
                number_of_samples = number_of_samples,   # requested
                n.used = length(indicators),             # actually contributing
                n.dropped = n.dropped))
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

    # The bootstrap gets its own safely() so that a failure here costs only the bootstrap
    # columns, not the model call above.
    #
    # It used to be a bare call, which meant one bad resample discarded a perfectly good
    # inference. infer.trio() calls class.vec(), which branches on a test statistic; if a
    # resample happens to draw a monomorphic V1 the genotype coefficient is aliased, the
    # statistic is NA, and class.vec() dies on "missing value where TRUE/FALSE needed".
    # parLapply reports that as "N nodes produced errors", the exception unwound all of
    # apply.mrgn(), and the row came back with model = NA. That is how 57 trios -- almost
    # all of them at n = 50, where a low-MAF SNP can lose every minor allele to resampling
    # -- ended up with no inferred model despite infer.trio() having succeeded on the real
    # data. See boostrap_edge_probabilities() for the resample-level handling.
    boot <- NULL
    bootstrap.time.seconds <- NA_real_
    bootstrap.error <- NA_character_
    if (bootstrap) {
        boot.start <- Sys.time()
        attempt <- safely(boostrap_edge_probabilities(data,
                                                      number_of_samples = number_of_samples,
                                                      cl = cl))
        boot <- attempt$value
        bootstrap.error <- attempt$error
        boot.end <- Sys.time()
        bootstrap.time.seconds <- as.numeric(difftime(boot.end, boot.start, units = "secs"))
    }

    final_result <- list(model = inferred_model,
                         adj = inferred_adj,
                         edge.probabilities = boot$edge.probabilities,
                         bootstrap = boot,
                         details = result,
                         time.seconds = as.numeric(difftime(end.time, start.time, units = "secs")),
                         bootstrap.time.seconds = bootstrap.time.seconds,
                         bootstrap.error = bootstrap.error)
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


# ---------------------------------------------------------------------------------------
# result files
# ---------------------------------------------------------------------------------------
#
# Each method writes one checkpoint per sample-size group, so the three never touch the
# same file and can run as concurrent processes. Combining is a separate step.

ALL.METHODS <- c("mrgn", "mrpc", "gmac", "mrggi")

method.checkpoint <- function(out.dir, method, size) {
    file.path(out.dir, paste0(method, "_group_n", size, ".RData"))
}


# Loop one method over its sample-size groups, checkpointing each. Shared by the three
# apply_*.R scripts, which differ only in the `runner` they pass in.
#
# runner: function(datasets, sel, cov.pool, known.conf, cl) -> one data.frame of results.
#         `sel` is NULL for methods that do not use the CS confounder sets (GMAC).
# needs.selection: FALSE skips loading the selection cache entirely, which saves several
#         hundred MB for GMAC.
run.method.groups <- function(method, runner, needs.selection = TRUE, cl = NULL,
                              verbose = TRUE) {

    sim_data <- loadRData(file = sim.data.file)
    group.of <- sapply(sim_data, function(x) x$params$sample.size)
    sizes <- if (is.null(sample.sizes)) sort(unique(group.of)) else sample.sizes
    if (verbose) cat("loaded", length(sim_data), "datasets |", length(sizes), "groups\n")

    for (size in sizes) {
        idx <- which(group.of == size)
        if (!is.null(max.per.group)) idx <- head(idx, max.per.group)
        datasets <- sim_data[idx]
        checkpoint <- method.checkpoint(out.dir, method, size)

        if (!rerun.inference && file.exists(checkpoint)) {
            if (verbose) cat("\n=== n =", size, "=== already checkpointed; skipping\n")
            next
        }
        if (verbose) cat("\n=== n =", size, "|", length(datasets), "trios ===\n")

        cov.pool <- group.cov.pool(datasets)
        known.conf <- group.known.conf(datasets)
        sel <- if (needs.selection) {
            get.selection(datasets, out.dir = out.dir, sim.file = sim.data.file,
                          rerun = rerun.selection,
                          # a subset run must never overwrite a full-group cache
                          save = is.null(max.per.group), verbose = verbose,
                          selection_fdr = selection_fdr, filter_fdr = filter_fdr,
                          alpha = alpha, filter_int_child = filter_int_child)
        } else NULL

        started <- Sys.time()
        results <- runner(datasets, sel, cov.pool, known.conf, cl)
        attr(results, "group.time.seconds") <-
            as.numeric(difftime(Sys.time(), started, units = "secs"))

        base::save(results, file = checkpoint)
        if (verbose) {
            cat("  saved", basename(checkpoint), "|", nrow(results), "rows |",
                round(attr(results, "group.time.seconds") / 60, 1), "min\n")
        }
        rm(sel, cov.pool, results); invisible(gc())
    }
    invisible(sizes)
}


# One method's per-group checkpoints -> inference_<method>.RData / .csv
combine.method <- function(method, sample.sizes = NULL, save.csv = TRUE, verbose = TRUE) {

    if (is.null(sample.sizes)) {
        files <- list.files(out.dir, pattern = paste0("^", method, "_group_n[0-9]+[.]RData$"),
                            full.names = TRUE)
    } else {
        files <- sapply(sample.sizes, function(s) method.checkpoint(out.dir, method, s))
        files <- files[file.exists(files)]
    }
    if (length(files) == 0) return(NULL)

    res <- do.call(rbind, lapply(files, loadRData))
    res <- res[order(res$dataset), ]
    row.names(res) <- NULL

    assign(paste0(method, ".results"), res)
    base::save(list = paste0(method, ".results"),
               file = file.path(out.dir, paste0("inference_", method, ".RData")))
    if (save.csv) {
        write.csv(res, file.path(out.dir, paste0("inference_", method, ".csv")),
                  row.names = FALSE)
    }
    if (verbose) {
        cat("  inference_", method, ": ", length(files), " groups, ", nrow(res),
            " rows, ", ncol(res), " columns\n", sep = "")
    }
    res
}


# Per-method files, then the master joined on `dataset`.
#
# `which.methods` defaults to every method the pipeline knows, so combining reflects what
# exists on disk rather than what one process happened to run. Otherwise a deliberate
# single-method run would rebuild the master from that method alone and silently drop the
# other two methods' columns from a previous full run.
#
# A method with no checkpoints is left out of the master rather than filled with NA
# columns, so the master always describes results that actually exist.
combine.all <- function(sample.sizes = NULL, save.csv = TRUE, verbose = TRUE,
                        which.methods = ALL.METHODS) {

    per.method <- lapply(which.methods, combine.method, sample.sizes = sample.sizes,
                         save.csv = save.csv, verbose = verbose)
    names(per.method) <- which.methods
    per.method <- Filter(Negate(is.null), per.method)
    if (length(per.method) == 0) stop("no method checkpoints found in ", out.dir)

    # the first table supplies the shared identifier and selection-score columns; from the
    # rest take only the columns the master does not already have, matched on dataset
    master <- per.method[[1]]
    for (nm in names(per.method)[-1]) {
        other <- per.method[[nm]]
        new.cols <- setdiff(colnames(other), colnames(master))
        if (length(new.cols) == 0) next
        master <- cbind(master, other[match(master$dataset, other$dataset), new.cols,
                                      drop = FALSE])
    }
    row.names(master) <- NULL
    inference.results <- master

    base::save(inference.results, file = file.path(out.dir, "inference_results.RData"))
    if (save.csv) {
        write.csv(inference.results, file.path(out.dir, "inference_results.csv"),
                  row.names = FALSE)
    }
    if (verbose) {
        cat("  inference_results: ", nrow(inference.results), " rows, ",
            ncol(inference.results), " columns, methods = ",
            paste(names(per.method), collapse = ", "), "\n", sep = "")
    }
    return(inference.results)
}


# ---------------------------------------------------------------------------------------
# MR-GGI
# ---------------------------------------------------------------------------------------
#
# MR-GGI is a Mendelian randomisation method: it uses a genetic instrument to estimate the
# causal effect of one gene on another, via the Wald ratio beta(V->outcome)/beta(V->exposure).
# See MRGGI_METHODS.md for the full derivation and the measurements behind the choices here.
#
# ONE CALL PER TRIO PER ARM. MRggi() takes a whole gene matrix and tests every pair, so a
# trio goes in as y = (T1, T2, covariates) rather than being called once per edge. The arms
# in mrggi.arms differ only in which covariates ride along in y.
#
# WHAT THE COVARIATES DO, since it is not what it looks like. MRggi's estimator is strictly
# pairwise -- .TSLS() reads only y[,i], y[,j], X[[i]] and X[[j]] -- so the T1-T2 Wald ratio
# is IDENTICAL in every arm. Measured on one trio: Bg1g2 = 0.808, p = 0.000 both with and
# without covariates in y. What changes is FDR_Bg1g2: MRggi adjusts each g1's p-values
# across that gene's pairs, so the T1-T2 p-value is now corrected for T1's pairs against
# every covariate too. That correction is the entire difference between the arms.
#
# The arms are therefore NOT a confounder adjustment. MRggi has no covariate argument and
# cannot adjust; the instrument is what is supposed to handle confounding, which is the
# method's premise. See MRGGI_METHODS.md section 4.
#
# THREE ADAPTATIONS ARE REQUIRED, all forced by there being ONE variant.
#
# 1. Only T1 gets the instrument; T2 and every covariate get a COLUMN OF ZEROS. MRggi's
#    .TSLS residualises the outcome on its own instruments before regressing on the
#    exposure's:
#         y2.resid = resid(lm(y2 ~ X2)); Bzy = coef(lm(y2.resid ~ X1))
#    OLS residuals are orthogonal to their own regressors, so if X1 == X2 -- which is what
#    FineMapping() produces for a trio, since V1 is the only variant and it is associated
#    with both genes -- Bzy is forced to exactly zero and every estimate collapses.
#    Measured against a true effect of 0.700:
#
#      T1 = V1, T2 = zeros, covariates = zeros   ->  0.651  (p = 0.000)
#      T1 = V1, T2 = V1,    covariates = zeros   ->  0.000  (p = 1.000)
#      T1 = V1, T2 = V1,    covariates = V1      ->  0.000  (p = 1.000)
#      T1 = V1, T2 = zeros, covariates = V1      ->  0.651  (p = 0.000)
#
#    Rows 1 and 4 are identical: what the COVARIATE columns hold cannot reach the T1-T2 row,
#    because that row reads only X[["T1"]] and X[["T2"]]. Zeros are used throughout so the
#    covariate rows stay estimable rather than collapsing the same way.
#
# 2. X must be POSITIONALLY ALIGNED with y and every element must be a real matrix.
#    X[[i]] is the instrument set for y[, i]. A list shorter than ncol(y) fails with
#    "subscript out of bounds"; a NULL element fails with "'data' must be of a vector type"
#    inside lapply(X, scale). A zero-COLUMN matrix (ncol 0) also fails, in .TSLS. A column
#    of zeros is the only encoding of "no instrument" that works -- not because it survives
#    scale() (scale() turns it into all-NaN) but because MRggi assigns scale.X and then
#    never reads it, so .TSLS gets the raw zeros.
#
# 3. colnames(y) MUST BE SET. MRggi builds its output with g1 = append(g1, colnames(y)[i]),
#    so NULL colnames leave g1 empty and it dies at its final data.frame() with
#    "arguments imply differing number of rows: 0, 1" -- an error that points nowhere near
#    the cause. cbind(T1 = a, T2 = b) on unnamed n x 1 MATRICES yields NULL colnames,
#    because cbind ignores the tag for matrix arguments. Hence the assertion below.
#
# Only the cis -> trans direction is reported. V1 is the cis gene's eQTL by construction, so
# T1 is the only gene with a legitimate instrument, and under this design T2 has none at all
# -- Bg2g1 is NaN by construction rather than merely unused. With a single instrument
# MRggi's test statistic reduces to Bzy/se_Bzy, the instrument->OUTCOME t-statistic, and the
# exposure's first stage cancels out: on a known-null trio it reported B = -38.42 at p = 0
# with a first-stage F of 0.1. Hence the mrggi.min.F gate.

# One trio, one arm. `covars` is a data.frame of covariate columns, or NULL for the bare
# trio. Returns the T1 -> T2 row of the MRggi output plus the diagnostics around it.
mrggi.one.trio <- function(trio, covars = NULL, cor.thr = mrggi.cor.thr,
                           p.adjust.method = mrggi.p.adjust) {

    V1 <- as.matrix(trio[, 1, drop = FALSE])
    genes <- as.data.frame(trio[, 2:3, drop = FALSE])

    # cbind.data.frame(df, NULL) is not a no-op -- NULL counts as a zero-row argument and it
    # fails with "arguments imply differing number of rows: n, 0". The `none` arm passes
    # NULL, so the empty case is branched on rather than relying on cbind to ignore it.
    covars <- if (is.null(covars)) NULL else as.data.frame(covars)
    if (!is.null(covars) && ncol(covars) == 0) covars <- NULL

    T_mat <- as.matrix(if (is.null(covars)) genes else cbind.data.frame(genes, covars))
    colnames(T_mat) <- c("T1", "T2", if (!is.null(covars)) colnames(covars) else NULL)

    # adaptation 3. Cheap, and the alternative is an error 5 frames deep that names none of
    # these variables.
    if (is.null(colnames(T_mat)) || anyNA(colnames(T_mat)) ||
        any(colnames(T_mat) == "") || anyDuplicated(colnames(T_mat))) {
        stop("mrggi.one.trio(): y needs unique, non-missing column names -- MRggi indexes ",
             "colnames(y) to build its output and fails opaquely without them")
    }

    # A constant covariate makes cor() return NA for its row, and MRggi's
    # corMat[which(abs(cor.y) > cor.thr)] <- 1 then subscripts with NA. Drop them; a column
    # with no variance carries no information for a correlation screen anyway.
    if (ncol(T_mat) > 2) {
        sds <- apply(T_mat[, -(1:2), drop = FALSE], 2, stats::sd)
        T_mat <- T_mat[, c(TRUE, TRUE, !is.na(sds) & sds > 0), drop = FALSE]
    }

    zeros <- matrix(0, nrow = nrow(T_mat), ncol = 1)
    X <- c(list(V1), rep(list(zeros), ncol(T_mat) - 1))   # adaptations 1 and 2
    names(X) <- colnames(T_mat)

    res <- MRggi::MRggi(y = T_mat, X = X, cor.thr = cor.thr,
                        p.adjust.method = p.adjust.method)

    # first-stage strength of V1 for T1, the exposure
    s <- summary(stats::lm(T_mat[, "T1"] ~ V1))
    F.T1 <- if (is.null(s$fstatistic)) NA_real_ else unname(s$fstatistic[1])

    # The T1-T2 row is absent when abs(cor(T1, T2)) <= cor.thr, because the pair never
    # enters calc.idx. At cor.thr = 0 that needs an exactly zero correlation, but the arms
    # are meant to survive a raised threshold too, so this is handled rather than indexed
    # past.
    hit <- which(res$g1 == "T1" & res$g2 == "T2")
    if (length(hit) != 1) {
        stop("MRggi returned ", length(hit), " T1-T2 rows (cor.thr = ", cor.thr,
             " screened the pair out?)")
    }

    list(B.T1T2 = res$Bg1g2[hit], p.T1T2 = res$pval_Bg1g2[hit],
         FDR.T1T2 = res$FDR_Bg1g2[hit], GGcor = res$GGcor[hit],
         F.T1 = F.T1, n.covars = ncol(T_mat) - 2, n.pairs = nrow(res))
}


mrggi.fields <- function(res, error) {
    if (is.null(res)) {
        return(list(B.T1T2 = NA_real_, p.T1T2 = NA_real_, FDR.T1T2 = NA_real_,
                    GGcor = NA_real_, F.T1 = NA_real_, n.covars = NA_integer_,
                    n.pairs = NA_integer_, edge = NA_character_,
                    edge.fdr = NA_character_, weak.instrument = NA,
                    time.seconds = NA_real_, error = error))
    }
    # Edge called on the cis -> trans direction only, and only when V1 is a usable
    # instrument for T1. No `correct` flag: MR-GGI resolves the T1-T2 edge, not a model
    # label, so model0/model3 and model1/model4 are indistinguishable to it. Same
    # convention as gmac.fields() -- the cross-tab belongs in the analysis stage.
    #
    # TWO EDGE CALLS, because the arms differ in exactly one thing.
    #
    #   edge      from the RAW p-value. Comparable with mrpc.alpha and selection.alpha,
    #             which are also applied per trio, and with GMAC's and MRGN's edge columns.
    #             IDENTICAL IN EVERY ARM by construction -- MRggi's estimator is pairwise,
    #             so the covariates in y cannot move B or p. This is the column to read
    #             against the other methods.
    #
    #   edge.fdr  from FDR.T1T2, MRggi's own multiplicity-adjusted p-value. This is the ONLY
    #             column that varies across the arms: the adjustment runs over T1's pairs
    #             with each covariate, so a wider covariate set is a harsher correction.
    #             This is the column to read across arms.
    #
    # Reporting only `edge` would make the four arms four identical tables. Reporting only
    # `edge.fdr` would make MR-GGI incomparable with GMAC and MRGN, which get no such
    # correction. Both are kept, and confusion_mrggi.R asserts that `edge` really is
    # arm-invariant -- if it ever is not, X was misaligned with y.
    weak <- !is.na(res$F.T1) && res$F.T1 < mrggi.min.F
    called <- !weak && !is.na(res$p.T1T2) && res$p.T1T2 < mrggi.alpha
    called.fdr <- !weak && !is.na(res$FDR.T1T2) && res$FDR.T1T2 < mrggi.alpha
    list(B.T1T2 = res$B.T1T2, p.T1T2 = res$p.T1T2, FDR.T1T2 = res$FDR.T1T2,
         GGcor = res$GGcor, F.T1 = res$F.T1,
         n.covars = res$n.covars, n.pairs = res$n.pairs,
         edge = if (called) "T1->T2" else "none",
         edge.fdr = if (called.fdr) "T1->T2" else "none",
         weak.instrument = weak,
         time.seconds = res$time.seconds, error = error)
}


# The covariate frame for one arm, or NULL for the bare trio. Mirrors how
# run.mrgn.group() picks its three inputs, but MR-GGI wants the covariates ALONE rather
# than a trio-plus-covariates frame, so the trio columns are dropped back off.
mrggi.arm.covars <- function(arm, dat, sel, i, K_n, index) {
    drop.trio <- function(x) if (ncol(x) > 3) x[, -(1:3), drop = FALSE] else NULL
    switch(arm,
           none  = NULL,
           truth = drop.trio(ground.truth.input(dat, K_n, index)),
           CSq   = drop.trio(sel$CS.q$trios.with.confs[[i]]),
           CSa   = drop.trio(sel$CS.alpha$trios.with.confs[[i]]),
           stop("unknown MR-GGI arm: ", arm))
}


# One trio, every arm, from a SLIM payload: the bare trio and one covariate frame per arm.
# This is the only function that crosses to a cluster worker, and the payload shape is the
# reason -- see the note in run.mrggi.group().
mrggi.run.arms <- function(task, arms = mrggi.arms) {
    fields <- list()
    for (arm in arms) {
        t0 <- Sys.time()
        r <- safely(suppressMessages(
            mrggi.one.trio(task$trio, covars = task$covars[[arm]])))
        if (!is.null(r$value)) {
            r$value$time.seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        }
        fields <- c(fields, prefixed(mrggi.fields(r$value, r$error), paste0("mrggi.", arm)))
    }
    fields
}


# `sel` is required whenever an arm needs it -- only the "none" arm can run without a
# confounder selection, so apply_mrggi.R sets needs.selection = TRUE.
run.mrggi.group <- function(datasets, sel, cov.names, cl = NULL, arms = mrggi.arms,
                            verbose = TRUE) {
    n.datasets <- length(datasets)

    if (any(c("CSq", "CSa") %in% arms) && is.null(sel)) {
        stop("run.mrggi.group(): arms ", paste(arms, collapse = ", "),
             " need a confounder selection, but sel is NULL")
    }

    # ---- why the work is packaged per trio rather than dispatched by index ----
    #
    # The obvious parallel form is parLapply(cl, seq_len(n), function(i) f(datasets[[i]],
    # sel, ...)). It is a memory trap. R serialises the closure's environment to EVERY
    # worker, so each one receives the whole `sel` and the whole `datasets` list. Measured
    # on the n = 1000 group: sel is 728 MB (of which $selection is 119 MB and CS-alpha's
    # $conf.list 247 MB, neither of which MR-GGI reads at all), so an 7-worker cluster would
    # copy ~5 GB before doing any work.
    #
    # Packaging one self-contained task per trio and passing the LIST to parLapply makes it
    # split the payload instead of replicating it: each worker gets only its own chunk, and
    # nothing that MR-GGI does not read is ever sent. ~1.3 MB per trio at n = 1000, so
    # ~55 MB per worker rather than 728 MB.
    #
    # The id and selection-score columns are built on the master. They are cheap, and
    # computing them here keeps `sel` off the workers entirely.
    tasks <- lapply(seq_len(n.datasets), function(i) {
        dat <- datasets[[i]]$data; params <- datasets[[i]]$params
        list(trio = dat[, 1:3],
             covars = setNames(lapply(arms, mrggi.arm.covars, dat = dat, sel = sel,
                                      i = i, K_n = params$K_n, index = params$dataset),
                               arms))
    })

    # The CS-alpha arm is ~5,800 gene pairs per trio and dominates the runtime (measured:
    # 9.4 h of the 10.4 h total across all groups and arms). Load-balanced rather than
    # block-split, because the per-trio cost scales with that trio's covariate count and so
    # varies a lot within a group.
    #
    # Dispatched in BATCHES rather than as one parLapplyLB over the whole group, purely so
    # there is something to watch. A single call returns only when the last trio finishes,
    # which at n = 1000 is roughly an hour of a log that says nothing -- indistinguishable
    # from a hang. Ten batches give a progress line and an ETA at ~10% granularity.
    # Load balancing still applies within each batch, which is where the variation is.
    fields <- if (!is.null(cl)) {
        n.batches <- max(1L, min(10L, n.datasets))
        bounds <- unique(as.integer(round(seq(0, n.datasets, length.out = n.batches + 1))))
        out <- list()
        t0 <- Sys.time()
        for (b in seq_len(length(bounds) - 1L)) {
            idx <- (bounds[b] + 1L):bounds[b + 1L]
            out <- c(out, parallel::parLapplyLB(cl, tasks[idx], mrggi.run.arms, arms = arms))
            if (verbose) {
                done <- bounds[b + 1L]
                el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
                cat(sprintf("    mrggi %d/%d | %.1f min elapsed, ~%.1f min left\n",
                            done, n.datasets, el, el / done * (n.datasets - done)))
                utils::flush.console()
            }
        }
        out
    } else {
        out <- vector("list", n.datasets)
        for (i in seq_len(n.datasets)) {
            out[[i]] <- mrggi.run.arms(tasks[[i]], arms = arms)
            if (verbose && (i %% 25 == 0 || i == n.datasets)) {
                cat("    mrggi", i, "/", n.datasets, "\n")
                utils::flush.console()
            }
        }
        out
    }
    rm(tasks); invisible(gc())

    rows <- lapply(seq_len(n.datasets), function(i) {
        dat <- datasets[[i]]$data; params <- datasets[[i]]$params
        index <- params$dataset
        # The CS scores describe the selected sets, so they only exist when there was a
        # selection. arms = c("none", "truth") is a legitimate configuration that needs
        # none, and asking for the scores anyway would fail on a NULL sel.
        as.data.frame(c(id.columns(dat, params, index),
                        if (!is.null(sel)) cs.score.columns(sel, i, index, dat, cov.names),
                        fields[[i]]),
                      stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}
