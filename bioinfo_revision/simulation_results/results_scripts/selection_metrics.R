# Confounder-selection performance: CS-q against CS-alpha.
#
# NOTE: Set your working directory to the repository root before running anything in this
# folder. e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
# Reads data/selection_results.csv -- 1,500 trios, written by run_confounder_selection.R --
# and scores each trio's selection against the confounders that trio actually has. Writes
# nothing; make_selection_report.R renders it.
#
# Deliberately base R, matching confusion_utils.R, so the selection numbers can be
# regenerated without a working MRGN install.
#
# THE TWO SELECTIONS. Both come from the same get.conf.trios()/select.confounders() call and
# differ only in the threshold applied to the regression p-values:
#
#   CS-q       Benjamini-Hochberg FDR at q < 0.05   (selection_fdr)
#   CS-alpha   the raw per-test cutoff alpha < 0.01 (alpha)
#
# That single difference is the whole of this comparison, and it is why CS-alpha's false
# positive count is a near-deterministic function of the pool: an uncorrected cutoff over a
# pool of N nulls returns alpha x N of them by construction. See POOL WIDTH below.

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

SELECTION.METHODS <- c(CSq = "CS-q", CSa = "CS-alpha")


load.selection <- function(results.dir = out.dir) {
    rdata <- file.path(results.dir, "selection_results.RData")
    csv   <- file.path(results.dir, "selection_results.csv")
    if (file.exists(rdata)) {
        env <- new.env(parent = emptyenv())
        load(rdata, envir = env)
        nm <- ls(env)
        if (length(nm) != 1) {
            stop(rdata, " holds ", length(nm), " objects; expected exactly one")
        }
        return(get(nm, envir = env))
    }
    if (file.exists(csv)) return(utils::read.csv(csv, stringsAsFactors = FALSE))
    stop("no selection results: neither ", rdata, " nor ", csv, " exists. ",
         "Run run_confounder_selection.R first.")
}


# ---------------------------------------------------------------------------------------
# POOL WIDTH
# ---------------------------------------------------------------------------------------
#
# Every trio contributes its own U/W/Z columns to the group's shared pool, and no covariate
# is a candidate for more than the trio that produced it (METHODS.md section 6, Table 21).
# So the pool a single trio is scored against is the WHOLE group's covariate count -- the
# column sum below -- and not that trio's own U_n.
#
# This is the denominator behind the true-negative cell, and it is what makes CS-alpha's
# behaviour arithmetic rather than statistical: at n = 670 the pool is 8,125 columns, of
# which a median of 25 are that trio's real confounders, so alpha = 0.01 returns ~81 false
# positives whatever the data does. Measured medians of CSa.n.fp are 80-90 across the five
# groups, within ~2% of alpha x pool at every sample size.
#
# K is excluded: the known covariates are passed to every model rather than selected, so
# they are never candidates and do not belong in the denominator.
selection.pool.width <- function(sel) {
    w <- vapply(split(sel, sel$sample.size),
                function(x) sum(x$U_n) + sum(x$W_n) + sum(x$Z_n), numeric(1))
    w[order(as.numeric(names(w)))]
}


# ---------------------------------------------------------------------------------------
# per-trio precision and recall
# ---------------------------------------------------------------------------------------
#
#   precision = tp / (tp + fp)      of what was selected, the share that was real
#   recall    = tp / (tp + fn)      of what was real, the share that was selected
#
# matching the definitions carried with the published study, in
# Manuscript/scripts/6_6_2023_all_conf_types_simulation_updated_2025.Rmd.
#
# BOTH CAN BE UNDEFINED, and neither case is exotic here.
#
# Precision is 0/0 when a trio selected nothing at all. At n = 50 the MEDIAN CS-q selection
# is empty, so this is the typical trio in that group rather than an edge case -- confounder
# selection has no power at n = 50 (METHODS.md section 4). Returning NA and reporting the
# undefined share beside the mean is the only honest treatment: scoring an empty selection
# as precision 0 says it selected wrongly, and as precision 1 that it selected perfectly,
# and it did neither.
#
# Recall is 0/0 for a trio with no true confounders at all. Rare, but present.
#
# n.selected is carried through unchanged rather than recomputed as tp + fp: the two differ
# when a selection picks up a covariate belonging to ANOTHER trio's block, which is scored
# in n.fp.other.trio and is a different kind of error from selecting a null column of the
# trio's own. Keeping both lets the report separate them.
selection.per.trio <- function(sel) {
    pools <- selection.pool.width(sel)
    out <- lapply(names(SELECTION.METHODS), function(p) {
        g  <- function(suffix) sel[[paste0(p, ".", suffix)]]
        tp <- g("n.tp"); fp <- g("n.fp"); fn <- g("n.fn")
        data.frame(
            dataset          = sel$dataset,
            sample.size      = sel$sample.size,
            effect_size      = sel$effect_size,
            truth.model      = sel$truth.model,
            n.true.confs     = sel$n.true.confs,
            pool.width       = as.numeric(pools[as.character(sel$sample.size)]),
            selection        = p,
            selection.lab    = unname(SELECTION.METHODS[p]),
            n.selected       = g("n.selected"),
            tp = tp, fp = fp, fn = fn,
            fp.other.trio    = g("n.fp.other.trio"),
            has.common.child = g("has.common.child"),
            has.intermediate = g("has.intermediate"),
            precision = ifelse(tp + fp > 0, tp / (tp + fp), NA_real_),
            recall    = ifelse(tp + fn > 0, tp / (tp + fn), NA_real_),
            stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}


# ---------------------------------------------------------------------------------------
# the selection confusion matrix, both views
# ---------------------------------------------------------------------------------------
#
# A 2x2 over covariates rather than over trios: every candidate column in the pool either is
# or is not a true confounder of the trio, and either was or was not selected for it. Counts
# are summed over the trios in the group, so the unit is covariate-trio pairs, which is what
# the selection procedure actually makes a decision about.
#
#   TN = pool.width - tp - fp - fn
#
# READ THIS TABLE WITH THE POOL WIDTH IN VIEW. TN dominates it by two orders of magnitude --
# ~8,000 of ~8,100 cells per trio are correctly-rejected nulls -- so specificity sits above
# 99% for both selections and separates them not at all. That is a property of a 0.3% signal
# density (METHODS.md Table 21), not a finding about either method. Precision and recall are
# the columns carrying information, which is why the report leads with them and presents
# this table as context rather than as the result.
SELECTION.TRUTH.LEVELS <- c("True confounder", "Not a confounder")
SELECTION.PRED.LEVELS  <- c("Selected", "Not selected")

selection.confusion.one <- function(rows) {
    tp <- sum(rows$tp); fp <- sum(rows$fp); fn <- sum(rows$fn)
    tn <- sum(rows$pool.width) - tp - fp - fn
    matrix(c(tp, fp, fn, tn), nrow = 2, ncol = 2,
           dimnames = list(truth     = SELECTION.TRUTH.LEVELS,
                           predicted = SELECTION.PRED.LEVELS))
}

# view = "size"   one matrix per sample size, effect treatments pooled
# view = "effect" one matrix per effect treatment, summed across sample sizes
selection.confusion <- function(per.trio, selection, view = c("size", "effect")) {
    view <- match.arg(view)
    x    <- per.trio[per.trio$selection == selection, , drop = FALSE]
    col  <- if (view == "size") "sample.size" else "effect_size"
    keys <- if (view == "size") SAMPLE.SIZES[SAMPLE.SIZES %in% x$sample.size]
            else                EFFECT.SIZES[EFFECT.SIZES %in% x$effect_size]
    out <- lapply(keys, function(k)
        selection.confusion.one(x[x[[col]] == k, , drop = FALSE]))
    names(out) <- as.character(keys)
    out
}


# ---------------------------------------------------------------------------------------
# summary scorecard
# ---------------------------------------------------------------------------------------
#
# Means are taken over the trios where the metric is DEFINED, with the undefined share
# reported beside them rather than folded in. n.undef.prec counts the trios that selected
# nothing; under CS-q at n = 50 that is most of the group, and a mean precision read without
# it is meaningless.
selection.scorecard <- function(per.trio, view = c("size", "effect")) {
    view <- match.arg(view)
    col  <- if (view == "size") "sample.size" else "effect_size"
    keys <- if (view == "size") SAMPLE.SIZES[SAMPLE.SIZES %in% per.trio$sample.size]
            else                EFFECT.SIZES[EFFECT.SIZES %in% per.trio$effect_size]

    rows <- list()
    for (p in names(SELECTION.METHODS)) {
        for (k in keys) {
            x <- per.trio[per.trio$selection == p & per.trio[[col]] == k, , drop = FALSE]
            if (nrow(x) == 0) next
            rows[[length(rows) + 1]] <- data.frame(
                selection    = unname(SELECTION.METHODS[p]),
                group        = as.character(k),
                trios        = nrow(x),
                pool.width   = round(mean(x$pool.width)),
                med.true     = stats::median(x$n.true.confs),
                med.selected = stats::median(x$n.selected),
                mean.prec    = mean(x$precision, na.rm = TRUE),
                mean.recall  = mean(x$recall,    na.rm = TRUE),
                n.undef.prec = sum(is.na(x$precision)),
                med.fp       = stats::median(x$fp),
                med.fp.other = stats::median(x$fp.other.trio),
                stringsAsFactors = FALSE)
        }
    }
    do.call(rbind, rows)
}
