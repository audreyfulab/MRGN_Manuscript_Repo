# Shared helpers for the confusion-matrix tables.
#
# NOTE: Set your working directory to the repository root before running anything in this
# folder. e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
# This file has no top-level side effects beyond sourcing inference_config.R, so it can be
# sourced from anywhere in the stage.
#
# Deliberately pure base R -- no MRGN, no ggplot2, no kable/xtable. Two reasons. First,
# inference_utils.R sources adapted_GMAC_func/* at its top level, which is a heavy and
# pointless dependency for a script that only counts labels, so the one helper we need from
# it (coarse.model) is re-declared below instead of sourced. Second, these tables are the
# artefact a reader is most likely to want to regenerate without a working MRGN install.
#
# inference_config.R IS sourced, for out.dir, so paths stay in step with the rest of the
# stage rather than being hardcoded a second time.

source("bioinfo_revision/simulation_results/inference_config.R")

tables.dir <- file.path(out.dir, "tables")


# ---------------------------------------------------------------------------------------
# label sets
# ---------------------------------------------------------------------------------------
#
# Truth rows are the five generating models (METHODS.md Table 2):
#
#   model0 -> M0.1   V1 -> T1, T2 independent
#   model1 -> M1.1   V1 -> T1 -> T2          (cis gene mediates)
#   model2 -> M2.1   V1 -> T1 <- T2          (trans gene mediates)
#   model3 -> M3     V1 -> T1, V1 -> T2
#   model4 -> M4     V1 -> T1 -> T2, V1 -> T2
#
# Truth is never a ".2" variant -- the cis gene is always T1 -- so collapsing the truth
# labels is a no-op and TRUTH.LEVELS is the same list either way. The MRGN predictions do
# use the .2 variants, and these tables report them coarsely, matching the collapse behind
# the mrgn.*.correct.coarse column.

TRUTH.LEVELS <- c("M0", "M1", "M2", "M3", "M4")

# MRGN returns MRGN::infer.trio()$Inferred.Model, which is one of the eight M-labels or
# "Other". There is no column for a failed fit, because there are none: every trio now has
# a model call. 114 of them used to be NA, not because MRGN could not classify the trio but
# because the bootstrap that runs after it threw and took the finished call with it; those
# were recovered by backfill_mrgn_models.R.
#
# If an NA ever shows up again, confusion() maps it to "Failed", which is not in this list,
# so the run stops with a named error rather than quietly dropping the trio. That is the
# intent -- a missing model is a bug to investigate, not a category to report.
MRGN.LEVELS  <- c("M0", "M1", "M2", "M3", "M4", "Other")

# GMAC reports a mediation call, not a model label -- see gmac.model.call() at
# inference_utils.R:885. The four levels are the four sign combinations of the cis and
# trans mediation p-values against selection.alpha; "Undirected" is the both-significant
# cell. There is no failure level here either: gmac.model has never been NA. As with MRGN,
# an NA would stop the run rather than quietly become a column.
GMAC.LEVELS  <- c("Cis Mediated", "Trans Mediated", "No Mediation", "Undirected")

# ...but the tables do not report those four. GMAC's statistic is the Wald test on the
# mediator coefficient in `outcome ~ mediator + treatment + confounders`, and that
# regression is symmetric in T1 and T2: under V1 -> T1 -> T2 the reverse coefficient is
# nonzero too, so the trans test fires whenever the cis one does. GMAC cannot orient the
# T1-T2 edge, and splitting "Cis Mediated" from "Trans Mediated" reports noise as a
# direction call -- across all 1500 trios those two cells hold 2 and 3 trios.
#
# So the three edge-present calls are pooled and GMAC is scored as what it can actually
# do: decide whether a T1-T2 edge is present. This is the framing the manuscript already
# used -- Manuscript/other/tablescraps/MRGN.GMAC.class.inference.50conf reports GMAC as
# "T1 - T2 Edge Absent" / "T1 - T2 Edge Present" against the same five truth models.
EDGE.LEVELS <- c("T1 - T2 Edge Absent", "T1 - T2 Edge Present")

gmac.edge <- function(call) {
    call <- as.character(call)
    out <- ifelse(call == "No Mediation", EDGE.LEVELS[1], EDGE.LEVELS[2])
    out[is.na(call)] <- NA_character_   # confusion() turns these into "Failed" and stops
    out
}

# MRGN is re-scored on the same two rows so the two methods can be read against each other
# on the one question they both answer. The edge status of an M-label follows from the
# topology and is the same mapping as for the truth side: M0 (V1 -> T1, T2 unconnected)
# and M3 (V1 -> T1, V1 -> T2) have no T1-T2 edge; M1, M2 and M4 do. The collapse is on the
# COARSE label, so M0.2 and M1.2 -- the mirrored variants, where the trans gene takes T1's
# place -- land on the same edge status as their .1 counterparts, which is right: the
# mirror swaps the two genes, it does not remove the edge between them.
#
# "Other" is kept as a third row rather than folded into either. It means infer.trio()
# fitted the trio but its four edge indicators matched none of the eight topologies, so it
# is not a wrong edge call -- it is no edge call. Assigning it one would be inventing
# output MRGN did not produce, and there are far too many to hide: 300 of 300 trios in the
# CS-alpha arm at n = 50, and 209 of 300 even in the oracle arm. The row makes that
# visible, and it is exactly the structural difference from GMAC, which always answers.
MRGN.EDGE.LEVELS <- c(EDGE.LEVELS, "Other")

mrgn.edge <- function(model) {
    model <- coarse.model(as.character(model))
    out <- rep("Other", length(model))
    out[model %in% c("M0", "M3")] <- EDGE.LEVELS[1]
    out[model %in% c("M1", "M2", "M4")] <- EDGE.LEVELS[2]
    out[is.na(model)] <- NA_character_  # confusion() turns these into "Failed" and stops
    out
}

# MR-GGI also scores on the T1-T2 edge, and for the same reason as GMAC: it resolves the
# edge, not a model label. Its estimator is a two-stage least squares fit of T2 on T1
# instrumented by V1, so model0/model3 and model1/model4 are indistinguishable to it --
# the pairs differ only in whether V1 also acts on T2 directly, which is a property of the
# trio topology and not of the T1-T2 relationship being estimated. Same two rows, same
# EDGE.CORRECT, so MR-GGI reads off the GMAC and MRGN edge tables directly.
#
# The third row is the structural difference. mrggi.fields() (inference_utils.R:1253)
# gates the edge call on the first-stage F for V1 -> T1 and records the trio as
# `edge = "none"` when F < mrggi.min.F, so the raw `none` level conflates two different
# things: the estimator ran and found no edge, and the estimator was never trusted to run
# at all. Those cannot share a row. Weak-instrument trios are 570 of 1,500 overall and 209
# of 300 at n = 50 -- filing them under edge-absent would hand MR-GGI an edge-absent
# recall it did not earn, most of it at the sample size where the claim matters most.
#
# So `weak.instrument` is broken out, exactly as MRGN's "Other" is: not a wrong edge call,
# no edge call. It counts against accuracy, because a trio whose edge was never called is
# a trio whose edge was not identified.
MRGGI.EDGE.LEVELS <- c(EDGE.LEVELS, "Weak instrument")

# Takes the two columns rather than one, because the gate lives in a separate column from
# the call. `weak` wins over `edge`: mrggi.fields() has already forced edge to "none" for
# those trios, so reading edge alone would silently score them as edge-absent.
mrggi.edge <- function(edge, weak) {
    edge <- as.character(edge)
    out <- rep(NA_character_, length(edge))
    out[edge == "none"]   <- EDGE.LEVELS[1]
    out[edge == "T1->T2"] <- EDGE.LEVELS[2]

    unknown <- setdiff(unique(edge[!is.na(edge)]), c("none", "T1->T2"))
    if (length(unknown)) {
        stop("unrecognised mrggi.edge value(s): ", paste(unknown, collapse = ", "),
             ". Only \"none\" and \"T1->T2\" are written -- see mrggi.fields().")
    }

    out[!is.na(weak) & weak] <- MRGGI.EDGE.LEVELS[3]
    out[is.na(edge)] <- NA_character_   # confusion() turns these into "Failed" and stops
    out
}

SAMPLE.SIZES <- c(50, 150, 300, 670, 1000)
EFFECT.SIZES <- c("small", "medium", "large")

# Same collapse as inference_utils.R:32, which is what mrgn.*.correct.coarse is built on.
# M0.1/M0.2 -> M0, M3 -> M3, Other -> Other. A no-op on the GMAC labels, none of which
# contain a dot.
coarse.model <- function(x) sub("\\..*$", "", x)


# ---------------------------------------------------------------------------------------
# what counts as correct
# ---------------------------------------------------------------------------------------
#
# Precision and recall both need to know, for each generating model, which inferred label
# is the right answer for it. That is a per-method mapping, and making it explicit is what
# lets one scoring function serve both a 6-row model table and a 2-row edge table.
#
# MRGN names the model, so the mapping is the identity on the five truth labels. "Other"
# is never the right answer for anything, so its row gets no precision.
MRGN.CORRECT <- c(M0 = "M0", M1 = "M1", M2 = "M2", M3 = "M3", M4 = "M4")

# On the edge tables the right answer is the edge status of the generating model
# (METHODS.md Table 2):
#
#   M0  V1 -> T1, T2 independent      no T1-T2 edge
#   M3  V1 -> T1, V1 -> T2            no T1-T2 edge
#   M1  V1 -> T1 -> T2                edge
#   M2  V1 -> T1 <- T2                edge
#   M4  V1 -> T1 -> T2, V1 -> T2      edge
#
# This is a property of the simulation, not of a method, so GMAC and the MRGN edge tables
# share it -- which is the whole point of scoring MRGN this way: same columns, same right
# answers, same margins, so the two are read off each other directly.
EDGE.CORRECT <- c(M0 = EDGE.LEVELS[1], M3 = EDGE.LEVELS[1],
                  M1 = EDGE.LEVELS[2], M2 = EDGE.LEVELS[2], M4 = EDGE.LEVELS[2])


# ---------------------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------------------

# Reads the combined results for one method. Prefers the .RData over the .csv: read.csv
# cannot tell a real NA (a fit that errored) from the literal string "NA", and that
# distinction is exactly what the last column of every matrix reports.
#
# The object inside inference_<method>.RData is named <method>.results -- see
# combine.method() at inference_utils.R:1072.
load.method <- function(method, results.dir = out.dir) {
    rdata <- file.path(results.dir, sprintf("inference_%s.RData", method))
    csv   <- file.path(results.dir, sprintf("inference_%s.csv", method))
    obj   <- paste0(method, ".results")

    if (file.exists(rdata)) {
        env <- new.env(parent = emptyenv())
        load(rdata, envir = env)
        if (!exists(obj, envir = env, inherits = FALSE)) {
            stop(sprintf("%s holds no object named '%s' (found: %s)",
                         rdata, obj, paste(ls(env), collapse = ", ")))
        }
        return(get(obj, envir = env))
    }

    if (file.exists(csv)) {
        cat(sprintf("  %s missing, falling back to %s\n", basename(rdata), basename(csv)))
        return(utils::read.csv(csv, stringsAsFactors = FALSE))
    }

    stop(sprintf("no results for method '%s': neither %s nor %s exists. ",
                 method, rdata, csv),
         "Run the inference stage first -- see ../README.md.")
}


# ---------------------------------------------------------------------------------------
# the confusion matrix itself
# ---------------------------------------------------------------------------------------

# Cross-tabulates true generating model against inferred label.
#
# Both margins are factors with FIXED levels, so every matrix this function returns has the
# same shape whether or not a category was observed. That matters: the tables are meant to
# be read across sample sizes, and a matrix that silently dropped an all-zero column would
# not line up with its neighbours. A NA prediction -- a fit that errored -- is folded into
# the literal "Failed" level rather than dropped, so the rows always sum to the number of
# trios.
#
# coarse.pred = FALSE for GMAC, whose labels have no sub-types to collapse.
confusion <- function(truth, pred, pred.levels, coarse.pred = TRUE) {
    truth <- factor(coarse.model(as.character(truth)), levels = TRUTH.LEVELS)
    if (any(is.na(truth))) {
        stop("unrecognised truth label(s); expected one of ",
             paste(TRUTH.LEVELS, collapse = ", "))
    }

    pred <- as.character(pred)
    if (coarse.pred) pred <- coarse.model(pred)
    pred[is.na(pred)] <- "Failed"
    unknown <- setdiff(unique(pred), pred.levels)
    if (length(unknown)) {
        stop("unrecognised inferred label(s): ", paste(unknown, collapse = ", "),
             if ("Failed" %in% unknown)
                 paste0(". A NA model means the fit never produced a call -- investigate",
                        " it rather than tabulating it; see backfill_mrgn_models.R.")
             else "")
    }
    pred <- factor(pred, levels = pred.levels)

    m <- table(truth, pred)
    dimnames(m) <- list(truth = TRUTH.LEVELS, predicted = pred.levels)
    m
}

# ---------------------------------------------------------------------------------------
# scoring: precision down the last column, recall along the bottom row
# ---------------------------------------------------------------------------------------
#
# confusion() returns truth in the rows because that is the natural way to build it. The
# tables are read the other way round -- inferred label in the rows, generating model in
# the columns -- which is the layout of the manuscript tables in
# Manuscript/other/tablescraps/MRGN.GMAC.class.inference.50conf. scored.table() transposes
# and adds the four margins:
#
#   Total column       trios given this inferred label
#   Precision column   of those, the share whose generating model maps to it
#   Total row          trios generated under this model
#   Recall row         of those, the share given the label that model maps to
#   bottom-right       overall accuracy, i.e. every correct cell over every trio
#
# correct.pred is the truth -> inferred mapping described above (MRGN.CORRECT /
# GMAC.CORRECT). Precision is defined per row as the share of that row falling in ANY
# column that maps to it, which is what makes the two-row GMAC table score correctly: the
# "T1 - T2 Edge Present" row is right for M1, M2 and M4 at once. An inferred label that is
# never correct for any model ("Other") has no precision and is left blank rather than
# scored 0 -- 0 would read as "always wrong", when the truth is that the question does not
# apply.
#
# Returns a character matrix, formatted and ready to write; counts stay integers and the
# rates are rounded to `digits`.
scored.table <- function(m, correct.pred, digits = 4) {
    d <- t(m)                       # rows = inferred label, columns = generating model
    truth.levels <- colnames(d)
    pred.levels  <- rownames(d)

    unmapped <- setdiff(truth.levels, names(correct.pred))
    if (length(unmapped)) {
        stop("correct.pred has no entry for truth label(s): ",
             paste(unmapped, collapse = ", "))
    }

    row.total <- rowSums(d)
    col.total <- colSums(d)

    # precision: for each inferred row, the columns it is the right answer for
    precision <- vapply(pred.levels, function(p) {
        cols <- truth.levels[correct.pred[truth.levels] == p]
        if (length(cols) == 0 || row.total[[p]] == 0) NA_real_
        else sum(d[p, cols]) / row.total[[p]]
    }, numeric(1))

    # recall: for each generating model, the single row it maps to
    recall <- vapply(truth.levels, function(t) {
        p <- unname(correct.pred[[t]])
        if (!(p %in% pred.levels) || col.total[[t]] == 0) NA_real_ else d[p, t] / col.total[[t]]
    }, numeric(1))

    correct <- sum(vapply(truth.levels, function(t) {
        p <- unname(correct.pred[[t]])
        if (p %in% pred.levels) d[p, t] else 0
    }, numeric(1)))
    accuracy <- if (sum(d) == 0) NA_real_ else correct / sum(d)

    rate <- function(x) ifelse(is.na(x), "", formatC(round(x, digits), format = "f",
                                                     digits = digits))
    out <- rbind(
        cbind(apply(d, 2, as.character), Total = as.character(row.total),
              Precision = rate(precision)),
        Total  = c(as.character(col.total), as.character(sum(d)), ""),
        Recall = c(rate(recall), "", rate(accuracy)))
    dimnames(out) <- list(inferred = c(pred.levels, "Total", "Recall"),
                          truth = c(truth.levels, "Total", "Precision"))
    out
}


# Long format, one row per cell. This is what make_all_tables.R stacks into
# confusion_counts_long.csv -- the single file to re-plot or re-aggregate from.
#
# `level` separates the two ways the same trios are tabulated: "model" is the five/six-way
# model call, "edge" the two/three-way T1-T2 edge call. Both exist for MRGN over the same
# rows, so without this column a (method, arm, sample_size, effect_size) group would hold
# two overlapping cross-tabs and any aggregation over it would double-count.
confusion.long <- function(m, method, arm, sample.size, effect.size, level = "model") {
    data.frame(method      = method,
               arm         = arm,
               level       = level,
               sample_size = sample.size,
               effect_size = effect.size,
               truth       = rep(rownames(m), times = ncol(m)),
               predicted   = rep(colnames(m), each  = nrow(m)),
               n           = as.vector(m),
               stringsAsFactors = FALSE)
}


# The margins of one edge table as numbers rather than formatted strings, for the
# method-by-method comparison table. m has truth in the rows and edge calls in the columns;
# the no-call column counts against accuracy, because a trio with no edge call is a trio
# whose edge was not identified.
#
# `no.call.level` names that column, since the two methods that have one name it
# differently and for different reasons -- MRGN's "Other" is a fitted trio matching no
# topology, MR-GGI's "Weak instrument" is a trio whose instrument failed the F gate. GMAC
# has no such column and passes the default, which simply never matches.
edge.scores <- function(m, no.call.level = "Other") {
    absent <- EDGE.LEVELS[1]; present <- EDGE.LEVELS[2]
    truth.absent  <- names(EDGE.CORRECT)[EDGE.CORRECT == absent]
    truth.present <- names(EDGE.CORRECT)[EDGE.CORRECT == present]
    called <- function(lvl) if (lvl %in% colnames(m)) sum(m[, lvl]) else 0
    hit <- function(rows, lvl) if (lvl %in% colnames(m)) sum(m[rows, lvl]) else 0

    safe <- function(num, den) if (den == 0) NA_real_ else num / den
    list(n              = sum(m),
         accuracy       = safe(hit(truth.absent, absent) + hit(truth.present, present), sum(m)),
         absent.prec    = safe(hit(truth.absent, absent), called(absent)),
         absent.recall  = safe(hit(truth.absent, absent), sum(m[truth.absent, ])),
         present.prec   = safe(hit(truth.present, present), called(present)),
         present.recall = safe(hit(truth.present, present), sum(m[truth.present, ])),
         no.call        = safe(called(no.call.level), sum(m)))
}


# ---------------------------------------------------------------------------------------
# writing
# ---------------------------------------------------------------------------------------

ensure.dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
}

CORNER <- "inferred \\ true"

# One CSV per arm rather than per sample size, with the five sample-size tables stacked
# inside it under `n = <size>` captions and separated by a blank line. Five files that have
# to be opened side by side to be compared is the wrong unit: the whole point of these
# tables is reading a method down the sample sizes, so the sample sizes belong in one file.
#
# Hand-rolled CSV rather than write.csv() because a stacked file is not one rectangle and
# write.csv() only knows how to write one. Fields are quoted only when they need it, which
# keeps the counts unquoted and diffable.
csv.field <- function(x) {
    if (grepl('[",\n]', x)) paste0('"', gsub('"', '""', x, fixed = TRUE), '"') else x
}
csv.row <- function(x) paste(vapply(as.character(x), csv.field, character(1)), collapse = ",")

# st: a scored.table() character matrix. Returns the block's lines, caption included.
scored.csv.block <- function(st, caption) {
    width <- ncol(st) + 1
    c(csv.row(c(caption, rep("", width - 1))),
      csv.row(c(CORNER, colnames(st))),
      vapply(seq_len(nrow(st)), function(i) csv.row(c(rownames(st)[i], st[i, ])),
             character(1)))
}

# blocks: named list of scored.table() matrices; the names become the captions.
write.scored.csv <- function(blocks, path) {
    lines <- unlist(lapply(seq_along(blocks), function(i) {
        c(scored.csv.block(blocks[[i]], names(blocks)[i]),
          if (i < length(blocks)) "" else character(0))
    }), use.names = FALSE)
    ensure.dir(dirname(path))
    writeLines(lines, path)
    invisible(path)
}

# Renders one scored table as a GitHub Markdown table. Hand-rolled because nothing else in
# bioinfo_revision pulls in knitr or xtable, and one table formatter is not worth breaking
# that. The Total and Recall rows are emphasised so the margins do not read as two more
# inferred labels.
md.scored.table <- function(st, caption = NULL, note = NULL) {
    header <- c(CORNER, colnames(st))
    margin <- rownames(st) %in% c("Total", "Recall")
    rows <- vapply(seq_len(nrow(st)), function(i) {
        cells <- c(rownames(st)[i], st[i, ])
        if (margin[i]) cells <- ifelse(nzchar(cells), paste0("**", cells, "**"), "")
        paste("|", paste(cells, collapse = " | "), "|")
    }, character(1))

    c(if (!is.null(caption)) c(paste("####", caption), ""),
      paste("|", paste(header, collapse = " | "), "|"),
      paste("|", paste(rep("---", length(header)), collapse = " | "), "|"),
      rows,
      "",
      if (!is.null(note)) c(note, ""))
}
