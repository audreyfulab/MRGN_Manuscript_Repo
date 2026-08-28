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
# inference_config.R IS sourced, for out.dir and tables.dir, so paths stay in step with the
# rest of the stage rather than being hardcoded a second time. Note that tables.dir is a
# sibling of out.dir, not a child of it: out.dir is `data/`, which holds only the generated
# .RData/.csv, and these tables are the report rather than the data.

source("bioinfo_revision/simulation_results/inference_config.R")


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

# MRPC returns the same eight M-labels plus "Other" -- classify.mrpc.adj()
# (inference_utils.R) matches the fitted 3x3 adjacency against the same topologies -- so it
# shares MRGN's levels with ONE addition.
#
# "Failed" IS A REAL LEVEL FOR MRPC, and this is the difference from MRGN. apply.mrpc()
# caps each fit at mrpc.timeout seconds and returns model = NA when it expires, and that is
# not rare: at 120 s it was 182 of 300 trios at n = 670 and 224 of 300 at n = 1000, in the
# CS-q arm alone. For MRGN an NA model means a bug to chase; for MRPC it means the fit did
# not finish in the time allowed, which is a property of the method on this design and one
# of the things the table is for. So it gets a column instead of stopping the run.
MRPC.LEVELS  <- c(MRGN.LEVELS, "Failed")

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

# MRPC is re-scored on the same rows, with the timeout column carried through. mrgn.edge()
# does the mapping for both -- the M-labels mean the same thing whichever method produced
# them -- and only the level set differs, because MRPC can fail to produce a label at all.
MRPC.EDGE.LEVELS <- c(EDGE.LEVELS, "Other", "Failed")

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
MRGGI.EDGE.LEVELS <- c(EDGE.LEVELS, "Weak instrument", "Screened out")

# Takes the two columns rather than one, because the gate lives in a separate column from
# the call. `weak` wins over `edge`: mrggi.fields() has already forced edge to "none" for
# those trios, so reading edge alone would silently score them as edge-absent.
mrggi.edge <- function(edge, weak) {
    edge <- as.character(edge)
    out <- rep(NA_character_, length(edge))
    out[edge == "none"]   <- EDGE.LEVELS[1]
    out[edge == "T1->T2"] <- EDGE.LEVELS[2]

    unknown <- setdiff(unique(edge[!is.na(edge)]), c("none", "T1->T2", "screened"))
    if (length(unknown)) {
        stop("unrecognised mrggi.edge value(s): ", paste(unknown, collapse = ", "),
             ". Only \"none\", \"T1->T2\" and \"screened\" are written -- see mrggi.fields().")
    }

    # Two distinct no-calls, kept apart because they say different things. "Weak instrument"
    # is a trio MR-GGI tried to test and could not -- the first-stage F gate. "Screened out"
    # is a trio it never looked at, because |cor(T1,T2)| fell below mrggi.cor.thr. Folding
    # either into "Edge Absent" would count a refusal to answer as a correct rejection, and
    # the screen removes 65.7% of model0 and 43.0% of model3 -- the two models whose right
    # answer IS "absent" -- so that error would flatter MR-GGI badly.
    out[edge == "screened"] <- MRGGI.EDGE.LEVELS[4]
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
# the full trio skeleton, for Figure 4
# ---------------------------------------------------------------------------------------
#
# EDGE.CORRECT above scores ONE edge -- T1-T2 -- because that is the only edge GMAC and
# MR-GGI resolve, and the two-row tables exist so those methods can be read against MRGN.
# Figure 4 is MRGN-only, so it has no such constraint and scores all three trio edges:
# V1-T1, V1-T2 and T1-T2.
#
# ADJACENCY, NOT ORIENTATION. An edge is a true positive when the pair is adjacent in both
# the generating topology and the inferred one, however either orients it. That is the
# consistent extension of the T1-T2 metric, which has always asked "is there an edge"
# rather than "which way does it point". The cost is that M1.1 (V1 -> T1 -> T2) and M2.1
# (V1 -> T1 <- T2) have the SAME skeleton, so a mediation trio called a v-structure scores
# a full 2/2 here. Orientation error is not what Figure 4 measures; it is what the model
# rows of Table 2 measure, and the two figures should be read together.
#
# The mapping is on the FINE label, not the coarse one -- unlike EDGE.CORRECT, which can
# collapse M0.1/M0.2 because the mirror does not change T1-T2 adjacency. Here it does:
# M0.1 is V1-T1 and M0.2 is V1-T2, which are different edges.
TRIO.EDGES <- list(
    "M0.1" = c("V1-T1"),
    "M0.2" = c("V1-T2"),
    "M1.1" = c("V1-T1", "T1-T2"),
    "M1.2" = c("V1-T2", "T1-T2"),
    "M2.1" = c("V1-T1", "T1-T2"),
    "M2.2" = c("V1-T2", "T1-T2"),
    "M3"   = c("V1-T1", "V1-T2"),
    "M4"   = c("V1-T1", "V1-T2", "T1-T2"))

# "Other" predicts the EMPTY skeleton, which is the same convention the two-row tables use
# for it: infer.trio() fitted the trio but matched no topology, so it named no edges. It
# therefore contributes a false negative for every true edge and never a false positive --
# a no-call costs recall and leaves precision untouched, exactly as at the model level.
# An NA model (the fit errored outright) is treated the same way.
trio.skeleton <- function(model) {
    model <- as.character(model)
    lapply(model, function(m) {
        if (is.na(m) || is.null(TRIO.EDGES[[m]])) character(0) else TRIO.EDGES[[m]]
    })
}


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
# The same grid as NUMBERS, unrounded, NA in the two cells that have no meaning. Split out
# from scored.table() because the Excel workbooks need real numeric cells: a sheet built
# from the character version puts text in every cell, so nothing in it can be summed, and
# the rates arrive already rounded to `digits`.
#
# scored.table() below is now only the formatter for this, so precision, recall and accuracy
# have one definition rather than two that can drift apart.
scored.table.values <- function(m, correct.pred) {
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

    # [Total, Precision] and [Recall, Total] are NA: a precision of the column totals and a
    # recall of the row totals are not quantities.
    out <- rbind(
        cbind(d, Total = row.total, Precision = precision),
        Total  = c(col.total, sum(d), NA_real_),
        Recall = c(recall, NA_real_, accuracy))
    dimnames(out) <- list(inferred = c(pred.levels, "Total", "Recall"),
                          truth = c(truth.levels, "Total", "Precision"))
    out
}

scored.table <- function(m, correct.pred, digits = 4) {
    v <- scored.table.values(m, correct.pred)
    nr <- nrow(v); nc <- ncol(v)

    # A cell is a RATE if it sits in the Precision column or the Recall row, and a COUNT
    # otherwise. The two NA cells fall in the rate region and format to "" either way.
    is.rate <- matrix(FALSE, nr, nc)
    is.rate[, nc] <- TRUE
    is.rate[nr, ] <- TRUE

    rate <- function(x) ifelse(is.na(x), "", formatC(round(x, digits), format = "f",
                                                     digits = digits))
    out <- ifelse(is.rate, rate(v), as.character(v))
    dim(out) <- dim(v)
    dimnames(out) <- dimnames(v)
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
# `no.call.level` names that column, since the methods that have one name it differently and
# for different reasons -- MRGN's "Other" is a fitted trio matching no topology, MR-GGI's
# "Weak instrument" is a trio whose instrument failed the F gate. GMAC has no such column
# and passes the default, which simply never matches.
#
# It takes a VECTOR because MRPC has two: "Other" (fitted, matched no topology) and "Failed"
# (the fit hit mrpc.timeout and never returned a graph). Both are no-calls and both count
# against accuracy, but they are different events and are kept as separate columns rather
# than pooled -- a timeout is a statement about cost, "Other" a statement about the fit.
edge.scores <- function(m, no.call.level = "Other") {
    absent <- EDGE.LEVELS[1]; present <- EDGE.LEVELS[2]
    truth.absent  <- names(EDGE.CORRECT)[EDGE.CORRECT == absent]
    truth.present <- names(EDGE.CORRECT)[EDGE.CORRECT == present]
    called <- function(lvl) {
        lvl <- intersect(lvl, colnames(m))
        if (length(lvl) == 0) 0 else sum(m[, lvl])
    }
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

# There is deliberately no per-arm CSV writer here. Every count these scripts produce goes
# into tables/confusion_counts_long.csv, one row per cell, and the pooled matrices are
# rendered for reading in tables/confusion_matrices.md. A stacked per-arm CSV is neither:
# it is not one rectangle, so read.csv() cannot load it, and it duplicates numbers the long
# file already holds. Reconstruct any matrix from the long file instead -- see
# matrix.from.long() in make_all_tables.R, which is what the Markdown report itself uses.

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


# ---------------------------------------------------------------------------------------
# the two views: by sample size, and by effect size
# ---------------------------------------------------------------------------------------
#
# confusion_counts_long.csv is keyed method x arm x level x sample_size x effect_size, with
# effect_size taking "all" alongside the three treatments. That single file already holds
# both views of every matrix; they differ only in which key is held fixed and which is
# summed over.
#
#   "size"    effect_size == "all", one matrix per sample size. The pooled view
#             make_all_tables.R renders.
#   "effect"  one matrix per effect treatment, summed ACROSS sample sizes. Not stored
#             anywhere, because it is an aggregation rather than a measurement -- every
#             cell is a raw count, so summing is exact and nothing is re-derived.
#
# Both are needed to read a result: a method can look stable across n while hiding a
# collapse at one effect size, and vice versa.

matrix.from.long <- function(counts, pred.levels, truth.levels = TRUTH.LEVELS) {
    m <- matrix(0L, nrow = length(truth.levels), ncol = length(pred.levels),
                dimnames = list(truth = truth.levels, predicted = pred.levels))
    if (nrow(counts) == 0) return(m)
    m[cbind(as.character(counts$truth), as.character(counts$predicted))] <-
        as.integer(counts$n)
    m
}

# Returns a named list of matrices: names are the sample sizes ("50", "150", ...) for
# view = "size" and the treatments ("small", "medium", "large") for view = "effect".
# Groups absent from `counts` are dropped rather than returned empty, so MRPC -- which
# stops at n = 300 -- yields three entries under "size" and not five.
#
# `check` verifies the identity the effect view depends on: for every sample size, the
# three treatment matrices must sum to the "all" matrix. If that fails the long file is
# internally inconsistent and every number downstream is suspect, so it stops rather than
# warns.
view.matrices <- function(counts, method, arm, level,
                          view = c("size", "effect"),
                          pred.levels, truth.levels = TRUTH.LEVELS, check = TRUE) {
    view <- match.arg(view)
    x <- counts[counts$method == method & counts$arm == arm & counts$level == level, ,
                drop = FALSE]
    if (nrow(x) == 0) return(list())

    if (check) {
        for (s in unique(x$sample_size)) {
            xs <- x[x$sample_size == s, , drop = FALSE]
            all.m <- matrix.from.long(xs[xs$effect_size == "all", , drop = FALSE],
                                      pred.levels, truth.levels)
            parts <- xs[xs$effect_size %in% EFFECT.SIZES, , drop = FALSE]
            if (nrow(parts) == 0) next
            sum.m <- matrix.from.long(
                stats::aggregate(n ~ truth + predicted, data = parts, FUN = sum),
                pred.levels, truth.levels)
            if (!identical(as.vector(all.m), as.vector(sum.m))) {
                stop(sprintf(paste("%s/%s/%s n=%s: the three effect-size matrices do not",
                                   "sum to the pooled one (%d vs %d trios). The long",
                                   "counts are inconsistent."),
                             method, arm, level, s, sum(sum.m), sum(all.m)))
            }
        }
    }

    if (view == "size") {
        x <- x[x$effect_size == "all", , drop = FALSE]
        keys <- SAMPLE.SIZES[SAMPLE.SIZES %in% x$sample_size]
        out <- lapply(keys, function(s)
            matrix.from.long(x[x$sample_size == s, , drop = FALSE],
                             pred.levels, truth.levels))
        names(out) <- as.character(keys)
    } else {
        x <- x[x$effect_size %in% EFFECT.SIZES, , drop = FALSE]
        keys <- EFFECT.SIZES[EFFECT.SIZES %in% x$effect_size]
        out <- lapply(keys, function(e)
            matrix.from.long(
                stats::aggregate(n ~ truth + predicted,
                                 data = x[x$effect_size == e, , drop = FALSE], FUN = sum),
                pred.levels, truth.levels))
        names(out) <- keys
    }
    out
}

# How a view's groups are labelled in a caption. Kept next to view.matrices() so the two
# cannot drift apart.
view.label <- function(view, key) {
    if (view == "size") paste0("n = ", key) else paste0(key, " effect")
}

# The inferred-label space for one method at one level. The report writers all need this and
# would otherwise each re-derive it from the data, which would silently accept a truncated
# level set whenever a label happens not to occur in some group -- a matrix missing its
# "Failed" column reads as a method that never failed.
#
# GMAC is the one that is NOT the obvious constant: it has no "Other", because gmac.edge()
# maps every mediation call onto one of the two edge states, so it uses the bare
# EDGE.LEVELS rather than MRGN.EDGE.LEVELS.
pred.levels.for <- function(method, level) {
    switch(paste(method, level),
           "mrgn model"     = MRGN.LEVELS,
           "mrgn edge"      = MRGN.EDGE.LEVELS,
           "mrpc model"     = MRPC.LEVELS,
           "mrpc edge"      = MRPC.EDGE.LEVELS,
           "gmac edge"      = EDGE.LEVELS,
           # MR-GGI's model call comes from MRGN::class.vec() and so uses MRGN's label set.
           "mrggi model"    = MRGN.LEVELS,
           "mrggi edge"     = MRGGI.EDGE.LEVELS,
           "mrggi edge.fdr" = MRGGI.EDGE.LEVELS,
           stop("no level set defined for method '", method, "' at level '", level, "'"))
}

# The truth -> correct-inferred-label map that scores a given level. Model levels score
# against the eight-topology map, edge levels against the edge map.
correct.for <- function(level) if (level == "model") MRGN.CORRECT else EDGE.CORRECT

# ---------------------------------------------------------------------------------------
# the confounder-structure runs
# ---------------------------------------------------------------------------------------
#
# The four covariate structures of METHODS.md section 7b. Shared by confusion_structures.R
# and make_structure_report.R so the directory map has one definition -- in particular that
# u_w_z IS the main simulation and reads from out.dir rather than from data_structures/,
# which is the detail most easily got wrong when copied.
STRUCTURE.SIZE <- 670

STRUCTURES <- list(
    u_only = list(dir = file.path(results.root, "data_structures", "u_only"),
                  label = "Confounders only",
                  detail = "W_n = 0, Z_n = 0, filter_int_child = FALSE"),
    u_w    = list(dir = file.path(results.root, "data_structures", "u_w"),
                  label = "Confounders + 1 intermediate",
                  detail = "W_n = 1, Z_n = 0, filter_int_child = TRUE"),
    u_z    = list(dir = file.path(results.root, "data_structures", "u_z"),
                  label = "Confounders + 1 common child",
                  detail = "W_n = 0, Z_n = 1, filter_int_child = TRUE"),
    u_w_z  = list(dir = out.dir,
                  label = "Confounders + intermediate + common child",
                  detail = "W_n = 1, Z_n = 1, filter_int_child = TRUE -- the main simulation"))

# The per-group checkpoint, not the combined inference_mrgn.RData: three of the four cases
# have only this one group, so there is nothing to combine and combine.method() was never
# run for them.
load.structure <- function(dir) {
    path <- file.path(dir, sprintf("mrgn_group_n%d.RData", STRUCTURE.SIZE))
    if (!file.exists(path)) return(NULL)
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    get("results", envir = env)
}


# ---------------------------------------------------------------------------------------
# markdown: table of contents
# ---------------------------------------------------------------------------------------
#
# The reports are long enough that a reader arriving from a link needs a way to move around
# them. Building the contents list FROM the assembled `lines` rather than maintaining it by
# hand is the whole point: rename a section, add one, reorder them, and the next run picks
# the change up. A hand-written contents block in a file stamped "do not edit by hand" would
# rot on the first heading change and nothing would catch it.
#
# md.anchor() reproduces GitHub's heading-slug rule, which is also what the VS Code preview
# and most static-site renderers use: take the RENDERED text (so markdown syntax and link
# targets are gone), lowercase it, drop every character that is not alphanumeric, space,
# hyphen or underscore, then turn spaces into hyphens. Runs of punctuation therefore leave
# runs of hyphens -- "Figure 3 -- selection" slugs to "figure-3----selection" -- which looks
# wrong and is correct.
md.anchor <- function(heading) {
    s <- sub("^#+[[:space:]]*", "", heading)
    s <- gsub("\\[([^]]*)\\]\\([^)]*\\)", "\\1", s)   # links render as their label
    s <- gsub("[`*_]", "", s)                         # code, bold and italic markers
    s <- tolower(s)
    s <- gsub("[^a-z0-9 _-]", "", s, perl = TRUE)     # everything else is dropped
    gsub(" ", "-", s)
}

# Headings inside fenced code blocks are shell comments and R comments, not headings, so the
# fence state has to be tracked rather than matching "^#" line by line.
md.headings <- function(lines, min.depth = 2, max.depth = 3) {
    fenced <- FALSE
    keep <- character(0)
    for (ln in lines) {
        if (grepl("^[[:space:]]*```", ln)) { fenced <- !fenced; next }
        if (fenced) next
        hashes <- regmatches(ln, regexpr("^#+", ln))
        if (!length(hashes)) next
        depth <- nchar(hashes)
        if (depth < min.depth || depth > max.depth) next
        keep <- c(keep, ln)
    }
    keep
}

md.toc <- function(lines, min.depth = 2, max.depth = 3) {
    heads <- md.headings(lines, min.depth, max.depth)
    if (!length(heads)) return(character(0))
    vapply(heads, function(h) {
        depth <- nchar(regmatches(h, regexpr("^#+", h)))
        label <- sub("^#+[[:space:]]*", "", h)
        label <- gsub("`", "", label)          # backticks inside a link label render badly
        label <- gsub("\\*\\*", "", label)
        sprintf("%s- [%s](#%s)", strrep("  ", depth - min.depth), label, md.anchor(h))
    }, character(1), USE.NAMES = FALSE)
}

# Splices a contents block in immediately above the first top-level section, leaving the
# title and the "generated by" preamble above it. Returns `lines` unchanged when the report
# has no headings to list, so a short report does not get an empty stub.
with.contents <- function(lines, heading = "## Contents", max.depth = 3) {
    toc <- md.toc(lines, 2, max.depth)
    if (!length(toc)) return(lines)
    at <- which(grepl("^## ", lines))
    if (!length(at)) return(lines)
    at <- at[1]
    c(lines[seq_len(at - 1)], heading, "", toc, "", "---", "", lines[at:length(lines)])
}


# ---------------------------------------------------------------------------------------
# markdown: orientation primer and shared glossary
# ---------------------------------------------------------------------------------------
#
# These reports are read by people arriving from a manuscript or a review, not only by
# whoever built the pipeline. One canonical definition per term, kept here and selected per
# report, is what stops five files drifting into five slightly different accounts of what
# "arm" or "no-call" means.
GLOSSARY <- list(
    trio = paste("The unit of analysis: one genetic variant `V1` and two genes, `T1` (the",
                 "*cis* gene, near the variant) and `T2` (the *trans* gene, far from it).",
                 "Every method is asked the same question about a trio -- how are these",
                 "three connected?"),
    models = paste("The five ways a trio can be wired, and the five right answers a method",
                   "is scored against. **M0** null, no `T1`-`T2` link. **M1** mediation,",
                   "`V1 -> T1 -> T2`. **M2** v-structure, `V1 -> T1 <- T2`. **M3**",
                   "conditional independence, `T1 <- V1 -> T2`. **M4** fully connected.",
                   "M1, M2 and M4 have a real `T1`-`T2` edge; M0 and M3 do not."),
    coarse = paste("Some models have two configurations (`M0.1` and `M0.2` differ in which",
                   "gene the variant acts on). Scoring is **coarse**: both count as `M0`."),
    levels = paste("Two questions, scored separately. The **model level** asks a method to",
                   "name which of M0-M4 generated the trio. The **edge level** asks only",
                   "whether the `T1`-`T2` link is there. GMAC and MR-GGI cannot name a",
                   "model, so the edge level is the only place all four methods compare."),
    recall = paste("Of the trios that really were generated under a model, the share the",
                   "method labelled correctly. Power, in other words."),
    precision = paste("Of the trios a method gave a label, the share that really were",
                      "generated that way. A method that answers rarely but well has high",
                      "precision and low recall, so the two are read together."),
    confounder = paste("A variable that affects both genes and so creates association",
                       "between them that is not causal. Adjusting for the right ones is",
                       "what lets a method tell a real edge from an induced one."),
    selection = paste("The step before inference: choosing, from a pool of several thousand",
                      "candidate covariates, which ones to adjust for. It is a separate",
                      "problem from inference and is scored separately."),
    arms = paste("The same method run on different confounder sets. Comparing arms",
                 "separates the cost of choosing confounders badly from the cost of the",
                 "method itself."),
    cs = paste("The three selection rules, differing only in how harshly they correct for",
               "testing many covariates at once. **CS-q** corrects across the whole",
               "trio-by-covariate matrix and picks fewest. **CS-alpha** applies a raw",
               "per-test cutoff with no correction and picks the most, nearly all wrong.",
               "**CS-i** corrects per covariate across trios; it is the rule the published",
               "GTEx analysis used."),
    oracle = paste("The **truth** or **oracle** arm is a method handed the real confounders.",
                   "It is a ceiling showing what the method could do if selection were",
                   "perfect, not a competitor, so it is normally excluded when marking the",
                   "best result in a row."),
    nocall = paste("A method can decline to answer. MRGN returns `Other` when no topology",
                   "matched, MRPC `Failed` when it timed out, MR-GGI `Weak instrument` when",
                   "its instrument was too weak to trust and `Screened out` when it never",
                   "looked at the trio. **A no-call is not a wrong call**: it lowers recall",
                   "without inflating anyone's precision, and a method that abstains often",
                   "can look precise for the wrong reason."),
    strata = paste("Every trio is drawn in one of three **effect-size strata** -- small,",
                   "medium or large -- which set how strong the variant-to-gene and",
                   "gene-to-gene effects are, and at one of five **sample sizes**",
                   "(n = 50, 150, 300, 670, 1000). n = 670 matches GTEx whole blood."),
    wz = paste("Two covariates that must **not** be adjusted for. The **intermediate** `W`",
               "sits on the causal path between the genes, so conditioning on it blocks the",
               "effect being estimated. The **common child** `Z` sits below both, so",
               "conditioning on it is a collider adjustment and opens a path that was not",
               "there."),
    methods = paste("**MRGN** is the method under study. **MRPC** infers a graph with a PC",
                    "algorithm. **GMAC** returns a mediation call rather than a model label.",
                    "**MR-GGI** is a Mendelian randomisation method using the variant as an",
                    "instrument."))

md.glossary <- function(keys, heading = "## Terms used here") {
    missing <- setdiff(keys, names(GLOSSARY))
    if (length(missing)) stop("no glossary entry for: ", paste(missing, collapse = ", "))
    labels <- c(trio = "Trio, `V1` / `T1` / `T2`", models = "M0 - M4",
                coarse = "Coarse labels", levels = "Model level / edge level",
                recall = "Recall", precision = "Precision", confounder = "Confounder",
                selection = "Confounder selection", arms = "Arm", cs = "CS-q / CS-alpha / CS-i",
                oracle = "Truth (oracle) arm", nocall = "No-call", strata = "Effect strata, sample sizes",
                wz = "Intermediate `W`, common child `Z`", methods = "The four methods")
    c(heading, "",
      "| Term | What it means |", "| --- | --- |",
      vapply(keys, function(k) sprintf("| **%s** | %s |", labels[[k]], GLOSSARY[[k]]),
             character(1), USE.NAMES = FALSE),
      "")
}

# Splices an orientation section and a glossary in above the first numbered section, so a
# reader who has never seen the pipeline can start at the top of any report. Call this
# BEFORE with.contents() so the new sections appear in the contents list.
with.primer <- function(lines, what, terms, what.heading = "## What this report is") {
    at <- which(grepl("^## ", lines))
    if (!length(at)) return(lines)
    at <- at[1]
    block <- c(what.heading, "", what, "", md.glossary(terms), "---", "")
    c(lines[seq_len(at - 1)], block, lines[at:length(lines)])
}
