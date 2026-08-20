# This script reruns the simulation inference with updated parameters and settings.
# following reviewer comments
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
# For every simulated trio it runs six inferences:
#   MRGN  with the true confounders (trio + K + that trio's own U variables)
#   MRGN  with the CS-q selected confounders
#   MRGN  with the CS-alpha selected confounders
#   MRPC  with the CS-q selected confounders
#   MRPC  with the CS-alpha selected confounders
#   GMAC  with the confounders GMAC selected for itself
# and consolidates them, with the simulation parameters and a scoring of each confounder
# selection, into one wide data.frame with one row per trio.
#
# Selection and GMAC both index a covariate pool by row, so the datasets are processed in
# groups sharing a sample size. Each group is checkpointed as it finishes.


# load libraries
library(MRGN)
library(MRPC)
library(parallel)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_utils.R")


# ---------------------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------------------

sim.data.file   <- "bioinfo_revision/simulation/simulated_data/simulated_trios.RData"
clinical.file   <- "./GTEx/data/kclist_top5_tiss.RData"
out.dir         <- "bioinfo_revision/simulation_results"
n.bootstrap     <- 1000    # bootstrap replicates per MRGN fit
mrpc.timeout    <- 120     # seconds; MRPC has taken hours on trios with many confounders
gmac.nperm      <- 1000
selection.alpha <- 0.05    # significance cutoff for the GMAC mediation calls
n.cores         <- max(1, parallel::detectCores() - 2)
max.per.group   <- NULL    # set to a small number to smoke test the pipeline

# the generating model of each scenario, as the label the methods report. The simulation
# always puts the cis gene in T1, so the truth is always the ".1" variant. An explicit
# lookup rather than MRGN::convert.truth(), which maps by sorted position and silently
# mislabels when the input does not contain all five models.
TRUTH.LABEL <- c(model0 = "M0.1", model1 = "M1.1", model2 = "M2.1",
                 model3 = "M3", model4 = "M4")

coarse.model <- function(x) sub("\\..*$", "", x)


# ---------------------------------------------------------------------------------------
# per-dataset helpers
# ---------------------------------------------------------------------------------------

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
# rank-deficient fits are expected and must not take down a 3750-trio run
safely <- function(expr) {
    tryCatch(list(value = expr, error = NA_character_),
             error = function(e) list(value = NULL, error = conditionMessage(e)))
}


# ---------------------------------------------------------------------------------------
# result rows
# ---------------------------------------------------------------------------------------

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
         error                  = error)
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

run.group <- function(datasets, cl = NULL, bootstrap = TRUE,
                      number_of_samples = n.bootstrap, timeout = mrpc.timeout,
                      nperm = gmac.nperm, verbose = TRUE) {

    n.datasets <- length(datasets)
    indices <- sapply(datasets, function(x) x$params$dataset)
    K_n <- datasets[[1]]$params$K_n

    just.trios <- lapply(datasets, function(x) x$data[, 1:3])

    # every dataset's own confounder block. The K block is shared by all datasets in the
    # group (the same clinical covariates), so it is passed as known.conf rather than
    # pooled, which would repeat identical columns once per dataset.
    cov.pool <- do.call(cbind.data.frame,
                        lapply(datasets, function(x) {
                            x$data[, (4 + x$params$K_n):ncol(x$data), drop = FALSE]
                        }))
    known.conf <- if (K_n > 0) datasets[[1]]$data[, 3 + seq_len(K_n), drop = FALSE] else NULL

    if (verbose) {
        cat("group n =", nrow(just.trios[[1]]), "|", n.datasets, "datasets |",
            ncol(cov.pool), "pooled covariates |",
            if (K_n > 0) paste(K_n, "known confounders") else "no known confounders", "\n")
    }

    # ---- confounder selection: CS-q and CS-alpha for MRGN and MRPC ----
    if (verbose) cat("selecting confounders (CS-q and CS-alpha)...\n")
    sel <- select.confounders(just.trios, cov.pool, known.conf = known.conf,
                              blocksize = min(500, n.datasets))

    # ---- GMAC, which selects its own confounders across the group ----
    if (verbose) cat("running GMAC over the group...\n")
    gmac.batch <- safely(run.gmac.all(just.trios, cov.pool, known.conf = known.conf,
                                      nperm = nperm, alpha = selection.alpha, cl = cl))
    if (!is.null(gmac.batch$error) && !is.na(gmac.batch$error)) {
        warning("GMAC failed for this group: ", gmac.batch$error)
    }
    gmac.out <- gmac.batch$value

    # ---- per-trio inference ----
    rows <- vector("list", n.datasets)
    for (i in seq_len(n.datasets)) {

        index <- indices[i]
        dat <- datasets[[i]]$data
        params <- datasets[[i]]$params
        truth.label <- unname(TRUTH.LABEL[params$model])
        true.confs <- own.block.names(colnames(dat), "U", index)

        # the confounder set each method was given
        csq.names <- colnames(cov.pool)[sel$CS.q$selection$sig.asso.covs[[i]]]
        csa.names <- colnames(cov.pool)[sel$CS.alpha$selection$sig.asso.covs[[i]]]
        gmac.names <- if (is.null(gmac.out)) character(0) else
            colnames(cov.pool)[which(gmac.out$sel.conf.ind[i, ] == 1)]

        # ---- MRGN, three confounder settings ----
        mrgn.truth <- safely(apply.mrgn(ground.truth.input(dat, K_n, index),
                                        bootstrap = bootstrap,
                                        number_of_samples = number_of_samples,
                                        cl = cl, verbose = FALSE))
        mrgn.csq <- safely(apply.mrgn(sel$CS.q$trios.with.confs[[i]],
                                      bootstrap = bootstrap,
                                      number_of_samples = number_of_samples,
                                      cl = cl, verbose = FALSE))
        mrgn.csa <- safely(apply.mrgn(sel$CS.alpha$trios.with.confs[[i]],
                                      bootstrap = bootstrap,
                                      number_of_samples = number_of_samples,
                                      cl = cl, verbose = FALSE))

        # ---- MRPC, two confounder settings ----
        mrpc.csq <- safely(apply.mrpc(sel$CS.q$trios.with.confs[[i]], timeout = timeout))
        mrpc.csa <- safely(apply.mrpc(sel$CS.alpha$trios.with.confs[[i]], timeout = timeout))

        # ---- GMAC: results from the batch run, timing from the per-trio path ----
        gmac.confs <- if (length(gmac.names) > 0) as.matrix(cov.pool[, gmac.names, drop = FALSE]) else NULL
        gmac.one <- safely(apply.gmac(just.trios[[i]], confounders = gmac.confs,
                                      known.conf = known.conf, nperm = nperm,
                                      alpha = selection.alpha))

        rows[[i]] <- as.data.frame(c(
            list(dataset = index),
            as.list(params[, setdiff(colnames(params), "dataset")]),
            list(truth.model = truth.label,
                 n.true.confs = length(true.confs),
                 true.confs = paste(true.confs, collapse = ";")),
            prefixed(score.selection(csq.names, index, true.confs), "CSq"),
            prefixed(score.selection(csa.names, index, true.confs), "CSa"),
            prefixed(score.selection(gmac.names, index, true.confs), "gmac"),
            list(CSq.filter_int_child = sel$CS.q$filter_int_child,
                 CSa.filter_int_child = sel$CS.alpha$filter_int_child),
            prefixed(mrgn.fields(mrgn.truth$value, mrgn.truth$error, truth.label), "mrgn.truth"),
            prefixed(mrgn.fields(mrgn.csq$value, mrgn.csq$error, truth.label), "mrgn.CSq"),
            prefixed(mrgn.fields(mrgn.csa$value, mrgn.csa$error, truth.label), "mrgn.CSa"),
            prefixed(mrpc.fields(mrpc.csq$value, mrpc.csq$error, truth.label), "mrpc.CSq"),
            prefixed(mrpc.fields(mrpc.csa$value, mrpc.csa$error, truth.label), "mrpc.CSa"),
            prefixed(gmac.fields(if (is.null(gmac.out)) NULL else gmac.out$results[i, ],
                                 gmac.one$value, gmac.one$error), "gmac")),
            stringsAsFactors = FALSE)

        if (verbose && (i %% 25 == 0 || i == n.datasets)) {
            cat("  ", i, "/", n.datasets, "trios\n")
        }
    }

    results <- do.call(rbind, rows)
    attr(results, "selection.time.seconds") <- sel$time.seconds
    attr(results, "gmac.time.seconds") <- if (is.null(gmac.out)) NA_real_ else gmac.out$time.seconds
    return(results)
}


# ---------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------

run.all.inference <- function(sample.sizes = NULL, bootstrap = TRUE,
                              number_of_samples = n.bootstrap, verbose = TRUE) {

    sim_data <- loadRData(file = sim.data.file)
    if (verbose) cat("loaded", length(sim_data), "simulated datasets\n")

    group.of <- sapply(sim_data, function(x) x$params$sample.size)
    all.sizes <- sort(unique(group.of))
    if (is.null(sample.sizes)) sample.sizes <- all.sizes

    # one cluster for the whole run, shared by GMAC and the MRGN bootstrap. gmac()
    # dispatches child.p/conf.fdr/getp.func through parLapply and those helpers live in a
    # sourced file rather than a package, so the workers have to source it themselves.
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    root <- normalizePath(getwd(), winslash = "/")
    parallel::clusterExport(cl, "root", envir = environment())
    invisible(parallel::clusterEvalQ(cl, {
        setwd(root)
        library(MRGN)
        source("adapted_GMAC_func/gmac_one_trio.R")
        source("adapted_GMAC_func/GMAC_moded/R/GMAC.R")
        NULL
    }))
    parallel::clusterSetRNGStream(cl, 234)
    if (verbose) cat("cluster started with", n.cores, "workers\n")

    for (size in sample.sizes) {
        idx <- which(group.of == size)
        if (!is.null(max.per.group)) idx <- head(idx, max.per.group)
        if (verbose) cat("\n=== sample size", size, "===\n")

        group.start <- Sys.time()
        results <- run.group(sim_data[idx], cl = cl, bootstrap = bootstrap,
                             number_of_samples = number_of_samples, verbose = verbose)
        group.end <- Sys.time()
        attr(results, "group.time.seconds") <-
            as.numeric(difftime(group.end, group.start, units = "secs"))

        checkpoint <- file.path(out.dir, paste0("inference_group_n", size, ".RData"))
        save(results, file = checkpoint)
        if (verbose) {
            cat("saved", checkpoint, "|", nrow(results), "rows |",
                round(attr(results, "group.time.seconds") / 60, 1), "min\n")
        }
    }

    invisible(sample.sizes)
}


# combine the per-group checkpoints into the single results data.frame
combine.groups <- function(sample.sizes = NULL, save.csv = TRUE) {

    files <- list.files(out.dir, pattern = "^inference_group_n[0-9]+\\.RData$",
                        full.names = TRUE)
    if (!is.null(sample.sizes)) {
        files <- file.path(out.dir, paste0("inference_group_n", sample.sizes, ".RData"))
        files <- files[file.exists(files)]
    }
    if (length(files) == 0) stop("no group checkpoints found in ", out.dir)

    groups <- lapply(files, loadRData)
    inference.results <- do.call(rbind, groups)
    inference.results <- inference.results[order(inference.results$dataset), ]
    row.names(inference.results) <- NULL

    save(inference.results, file = file.path(out.dir, "inference_results.RData"))
    if (save.csv) {
        write.csv(inference.results, file.path(out.dir, "inference_results.csv"),
                  row.names = FALSE)
    }
    cat("combined", length(files), "groups into", nrow(inference.results), "rows\n")
    return(inference.results)
}


run.all.inference()
inference.results <- combine.groups()
