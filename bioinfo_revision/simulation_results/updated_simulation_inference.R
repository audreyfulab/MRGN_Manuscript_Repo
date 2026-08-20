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

# Which methods to run. Each writes its own per-group checkpoint and its own combined
# results file, so a single method can be re-run without disturbing the others.
methods <- c("mrgn", "mrpc", "gmac")

# Confounder selection is cached per sample-size group by run_confounder_selection.R.
# FALSE reuses the cache when it matches this request; TRUE recomputes and overwrites it.
# A cache built against different simulated data is detected and recomputed either way.
rerun.selection <- FALSE

# FALSE skips a method whose checkpoint for that group already exists. Note this applies
# to GMAC too: a completed GMAC group is not re-run when you later run MRGN alone. Set
# TRUE (or delete the checkpoints) to force everything to run again.
rerun.inference <- FALSE

# selection settings, kept in step with run_confounder_selection.R
selection_fdr    <- 0.05
filter_fdr       <- 0.1
alpha            <- 0.01
filter_int_child <- TRUE

# TRUTH.LABEL, coarse.model(), own.block.names(), ground.truth.input(), score.selection(),
# prefixed() and safely() now live in inference_utils.R, which run_confounder_selection.R
# also sources.


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


# ---- MRPC: CS-q and CS-alpha ----
run.mrpc.group <- function(datasets, sel, cov.names, timeout = mrpc.timeout, verbose = TRUE) {
    n.datasets <- length(datasets)
    rows <- vector("list", n.datasets)
    for (i in seq_len(n.datasets)) {
        dat <- datasets[[i]]$data; params <- datasets[[i]]$params
        index <- params$dataset
        truth.label <- unname(TRUTH.LABEL[params$model])

        mrpc.csq <- safely(apply.mrpc(sel$CS.q$trios.with.confs[[i]], timeout = timeout))
        mrpc.csa <- safely(apply.mrpc(sel$CS.alpha$trios.with.confs[[i]], timeout = timeout))

        rows[[i]] <- as.data.frame(c(
            id.columns(dat, params, index),
            cs.score.columns(sel, i, index, dat, cov.names),
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

        rows[[i]] <- as.data.frame(c(
            id.columns(dat, params, index),
            prefixed(score.selection(gmac.names, index, true.confs), "gmac"),
            prefixed(gmac.fields(if (is.null(gmac.out)) NULL else gmac.out$results[i, ],
                                 gmac.one$value, gmac.one$error), "gmac"),
            list(gmac.batch.error = batch.error)),
            stringsAsFactors = FALSE)

        if (verbose && (i %% 25 == 0 || i == n.datasets)) cat("    gmac", i, "/", n.datasets, "\n")
    }
    results <- do.call(rbind, rows)
    attr(results, "gmac.time.seconds") <-
        if (is.null(gmac.out)) NA_real_ else gmac.out$time.seconds
    results
}


method.checkpoint <- function(out.dir, method, size) {
    file.path(out.dir, paste0(method, "_group_n", size, ".RData"))
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
        datasets <- sim_data[idx]
        if (verbose) {
            cat("\n=== sample size", size, "|", length(datasets), "trios ===\n")
        }

        # which methods still need running for this group
        todo <- Filter(function(m) rerun.inference ||
                           !file.exists(method.checkpoint(out.dir, m, size)), methods)
        if (length(todo) == 0) {
            if (verbose) cat("  all requested methods already checkpointed; skipping\n")
            next
        }
        if (verbose && length(todo) < length(methods)) {
            cat("  already done:", paste(setdiff(methods, todo), collapse = ", "), "\n")
        }

        # Selection is needed by MRGN and MRPC but not by GMAC, which picks its own
        # confounders. Loading it costs several hundred MB, so only pay for it if one of
        # those two is actually going to run.
        sel <- NULL
        cov.pool <- group.cov.pool(datasets)
        known.conf <- group.known.conf(datasets)
        if (any(c("mrgn", "mrpc") %in% todo)) {
            sel <- get.selection(datasets, out.dir = out.dir, sim.file = sim.data.file,
                                 rerun = rerun.selection,
                                 save = is.null(max.per.group), verbose = verbose,
                                 selection_fdr = selection_fdr, filter_fdr = filter_fdr,
                                 alpha = alpha, filter_int_child = filter_int_child)
        }

        for (m in todo) {
            m.start <- Sys.time()
            results <- switch(m,
                mrgn = run.mrgn.group(datasets, sel, colnames(cov.pool), cl = cl,
                                      bootstrap = bootstrap,
                                      number_of_samples = number_of_samples,
                                      verbose = verbose),
                mrpc = run.mrpc.group(datasets, sel, colnames(cov.pool),
                                      timeout = mrpc.timeout, verbose = verbose),
                gmac = run.gmac.group(datasets, cov.pool, known.conf, cl = cl,
                                      nperm = gmac.nperm, verbose = verbose),
                stop("unknown method: ", m))
            attr(results, "group.time.seconds") <-
                as.numeric(difftime(Sys.time(), m.start, units = "secs"))

            checkpoint <- method.checkpoint(out.dir, m, size)
            base::save(results, file = checkpoint)
            if (verbose) {
                cat("  saved", basename(checkpoint), "|", nrow(results), "rows |",
                    round(attr(results, "group.time.seconds") / 60, 1), "min\n")
            }
        }
        rm(sel, cov.pool); invisible(gc())
    }

    invisible(sample.sizes)
}


# ---------------------------------------------------------------------------------------
# combining
# ---------------------------------------------------------------------------------------

# One method's per-group checkpoints -> inference_<method>.RData / .csv
combine.method <- function(method, sample.sizes = NULL, save.csv = TRUE, verbose = TRUE) {

    if (is.null(sample.sizes)) {
        files <- list.files(out.dir, pattern = paste0("^", method, "_group_n[0-9]+\\.RData$"),
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


# every method this pipeline knows how to run, independent of which are configured now
ALL.METHODS <- c("mrgn", "mrpc", "gmac")

# Per-method files, then the master joined on `dataset`.
#
# `which.methods` defaults to every method with checkpoints on disk, NOT to the `methods`
# config. Those control what gets *run*; combining is about what exists. Otherwise a
# deliberate single-method run (methods <- "mrpc") would rebuild the master from MRPC
# alone and silently drop the MRGN and GMAC columns from a previous full run.
#
# A method with no checkpoints is left out of the master rather than filled with NA
# columns, so the master always describes results that actually exist.
combine.groups <- function(sample.sizes = NULL, save.csv = TRUE, verbose = TRUE,
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


run.all.inference()
cat("\ncombining:\n")
inference.results <- combine.groups()
