# The aggregated confusion matrices for MRGN (CS-q) and MRPC (CS-q), as one workbook.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/make_aggregated_confusion_workbook.R
#
# Writes bioinfo_revision/Final_tables/Aggregated_confusion_matrices.xlsx. Reads
# tables/confusion_counts_long.csv, so run make_all_tables.R first.
#
# WHAT "AGGREGATED" MEANS HERE. make_matrix_workbooks.R writes one workbook per
# method-arm-GROUP: five sample sizes and three effect sizes, each its own matrix. This
# script writes the pooled matrix instead -- every sample size and every effect size summed
# into one table, which is the same pooling the "Aggregated" sheet of
# Final_tables/Final_rec_prec_tables.xlsx reports rates from. The counts are summed and the
# margins computed once from the pooled matrix (micro-averaging), NOT averaged across the
# per-cell rates; those two differ whenever the cells hold different numbers of trios, and
# the pooled matrix is the one whose margins can be read straight off the cells.
#
# The pooled matrix is built from the effect_size == "all" rows only. Those already hold
# each sample size's three treatments summed, so adding the "small"/"medium"/"large" rows on
# top would double-count every trio. view.matrices() checks that identity on the way past.
#
# WHY MRGN APPEARS TWICE. MRPC was run only at n = 50, 150 and 300 -- apply.mrpc() caps each
# fit at mrpc.timeout seconds and the larger designs were not run -- so its pooled matrix
# covers 900 trios where MRGN's covers 1,500. The two are not on the same footing: MRPC's
# 900 are the three SMALLEST sample sizes, which is where every method is weakest, so
# comparing them straight across understates MRPC. The workbook therefore carries a third
# MRGN block restricted to n <= 300, which is the like-for-like comparison; the full-run MRGN
# block is kept because it, not the restricted one, is the number the rest of the manuscript
# reports.
#
# THE THREE FORMATS, as in make_matrix_workbooks.R:
#
#   model      inferred model  x generating model
#   edge       T1-T2 edge call x generating model
#   edge_2x2   T1-T2 edge call x true T1-T2 edge state
#
# The first two are the `model` and `edge` levels of confusion_counts_long.csv. The third is
# derived here by collapsing the truth axis of the second through EDGE.CORRECT. Rows are the
# inferred label and columns the truth in all three -- scored.table.values() transposes -- so
# a Total/Precision column and a Total/Recall row come with each.

library(writexl)

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

final.tables.dir <- "bioinfo_revision/Final_tables"
out.file <- file.path(final.tables.dir, "Aggregated_confusion_matrices.xlsx")

counts <- utils::read.csv(file.path(tables.dir, "confusion_counts_long.csv"),
                          stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------------------
# pooling
# ---------------------------------------------------------------------------------------

# One pooled matrix. `sizes` restricts which sample-size groups go in; the default is all of
# them, and the matched MRGN block passes the three MRPC ran at. Only the effect_size ==
# "all" rows are summed -- see the header note on double-counting.
pooled.matrix <- function(method, arm, level, sizes = SAMPLE.SIZES) {
    x <- counts[counts$method == method & counts$arm == arm & counts$level == level &
                counts$effect_size == "all" & counts$sample_size %in% sizes, , drop = FALSE]
    if (nrow(x) == 0) {
        stop(sprintf("no rows for %s/%s at level '%s' in n = {%s}",
                     method, arm, level, paste(sizes, collapse = ", ")))
    }
    matrix.from.long(stats::aggregate(n ~ truth + predicted, data = x, FUN = sum),
                     pred.levels.for(method, level))
}

# The same identity view.matrices() checks, applied to what this script actually sums: for
# every sample size in the block, the three effect-size matrices must add up to the "all"
# one. If they do not, the long file is inconsistent and every number below is suspect.
# Also refuses to pool a block whose sample sizes are not all present, so a silently
# truncated run cannot be reported as the full one.
check.pooling <- function(method, arm, level, sizes) {
    invisible(view.matrices(counts, method, arm, level, view = "size",
                            pred.levels = pred.levels.for(method, level), check = TRUE))
    x <- counts[counts$method == method & counts$arm == arm & counts$level == level, ,
                drop = FALSE]
    missing <- setdiff(sizes, unique(x$sample_size))
    if (length(missing)) {
        stop(sprintf("%s/%s has no n = %s group", method, arm,
                     paste(missing, collapse = ", ")))
    }
}

# Collapses the generating-model rows of an edge table into the two true edge states.
# Counts, so summing is exact. Same helper as make_matrix_workbooks.R.
collapse.truth.to.edge <- function(m) {
    absent  <- names(EDGE.CORRECT)[EDGE.CORRECT == EDGE.LEVELS[1]]
    present <- names(EDGE.CORRECT)[EDGE.CORRECT == EDGE.LEVELS[2]]
    pick <- function(lv) {
        rows <- intersect(lv, rownames(m))
        if (length(rows) == 0) rep(0, ncol(m)) else colSums(m[rows, , drop = FALSE])
    }
    out <- rbind(pick(absent), pick(present))
    dimnames(out) <- list(truth = EDGE.LEVELS, predicted = colnames(m))
    out
}

# For the edge-by-edge format the truth axis and the correct inferred label are the same two
# states, so the map is the identity.
EDGE.IDENTITY <- stats::setNames(EDGE.LEVELS, EDGE.LEVELS)

# A scored grid as a data.frame ready for a sheet: numeric cells, row labels promoted to a
# real column because row names do not survive into a worksheet.
sheet.of <- function(m, correct.pred) {
    v <- scored.table.values(m, correct.pred)
    d <- as.data.frame(v, stringsAsFactors = FALSE, check.names = FALSE)
    cbind(inferred = rownames(v), d, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------------------
# the three blocks
# ---------------------------------------------------------------------------------------

BLOCKS <- list(
    list(key = "MRGN_CSq", method = "mrgn", arm = "CSq", sizes = SAMPLE.SIZES,
         label = "MRGN (CS-q)",
         scope = "all five sample sizes (n = 50, 150, 300, 670, 1000), all effect sizes"),
    list(key = "MRPC_CSq", method = "mrpc", arm = "CSq", sizes = c(50, 150, 300),
         label = "MRPC (CS-q)",
         scope = "n = 50, 150, 300 -- the only sizes MRPC was run at; all effect sizes"),
    list(key = "MRGN_CSq_n300", method = "mrgn", arm = "CSq", sizes = c(50, 150, 300),
         label = "MRGN (CS-q), matched to MRPC",
         scope = "n = 50, 150, 300 -- restricted to match MRPC; all effect sizes"))

sheets <- list()
summary.rows <- list()

for (b in BLOCKS) {
    check.pooling(b$method, b$arm, "model", b$sizes)

    mm <- pooled.matrix(b$method, b$arm, "model", b$sizes)
    em <- pooled.matrix(b$method, b$arm, "edge",  b$sizes)
    e2 <- collapse.truth.to.edge(em)

    if (sum(mm) != sum(em)) {
        stop(sprintf("%s: the model and edge matrices cover different trios (%d vs %d)",
                     b$key, sum(mm), sum(em)))
    }

    sheets[[paste0(b$key, "_model")]]    <- sheet.of(mm, correct.for("model"))
    sheets[[paste0(b$key, "_edge")]]     <- sheet.of(em, EDGE.CORRECT)
    sheets[[paste0(b$key, "_edge_2x2")]] <- sheet.of(e2, EDGE.IDENTITY)

    ms <- scored.table.values(mm, correct.for("model"))
    es <- scored.table.values(e2, EDGE.IDENTITY)
    summary.rows[[b$key]] <- data.frame(
        block          = b$label,
        trios          = sum(mm),
        scope          = b$scope,
        model_accuracy = ms["Recall", "Precision"],
        no_call_rate   = 1 - sum(mm[, intersect(TRUTH.LEVELS, colnames(mm))]) / sum(mm),
        edge_accuracy  = es["Recall", "Precision"],
        edge_recall    = es["Recall", EDGE.LEVELS[2]],
        edge_precision = es[EDGE.LEVELS[2], "Precision"],
        stringsAsFactors = FALSE)

    cat(sprintf("  %-30s %5d trios | model acc %.3f | edge acc %.3f\n",
                b$label, sum(mm), ms["Recall", "Precision"], es["Recall", "Precision"]))
}

summary.tbl <- do.call(rbind, summary.rows)
rownames(summary.tbl) <- NULL

# ---------------------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------------------

stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

readme <- data.frame(
    field = c(
        "What this workbook is",
        "Pooling",
        "Layout",
        "Sheet suffix _model",
        "Sheet suffix _edge",
        "Sheet suffix _edge_2x2",
        "Total column / Total row",
        "Precision column",
        "Recall row",
        "Bottom-right cell",
        "Blank precision",
        "Inferred label: Other",
        "Inferred label: Failed",
        "MRGN blocks",
        "Arm CS-q",
        "summary sheet",
        "Generated",
        "Source",
        "Script"),
    value = c(
        paste("Aggregated (pooled) confusion matrices for MRGN and MRPC over the simulated",
              "trio run, one matrix per sheet."),
        paste("Counts are summed over every sample-size x effect-size cell in the block and",
              "the margins computed once from the pooled matrix (micro-averaged), not",
              "averaged across per-cell rates."),
        paste("Rows are the INFERRED label, columns the GENERATING model -- the layout of",
              "the manuscript tables. Column A holds the row labels."),
        "Inferred model (M0-M4, Other, and Failed for MRPC) x generating model.",
        "T1-T2 edge call x generating model.",
        paste("T1-T2 edge call x true T1-T2 edge state: the truth axis of the _edge sheet",
              "collapsed, M0 and M3 to edge-absent, M1/M2/M4 to edge-present. The no-call",
              "rows (Other, Failed) are left alone -- they are not edge calls and must not",
              "be folded into either edge state."),
        paste("Trios given that inferred label (column) / generated under that model (row).",
              "Every cell is a count, so the sheets can be summed and re-aggregated."),
        "Of the trios given that inferred label, the share whose generating model maps to it.",
        "Of the trios generated under that model, the share given the label it maps to.",
        "Overall accuracy: every correct cell over every trio in the sheet.",
        paste("An inferred label that is never the right answer for any model (Other,",
              "Failed) has no precision. Blank, not 0 -- 0 would read as 'always wrong'",
              "when the truth is that the question does not apply."),
        paste("The method fitted the trio but its edge indicators matched none of the eight",
              "topologies. A no-call, not a wrong call: it costs recall and leaves every",
              "precision untouched."),
        paste("MRPC only. The fit hit mrpc.timeout and never returned a graph. Kept apart",
              "from Other because a timeout is a statement about cost, not about the fit."),
        paste("Two, and they are not interchangeable. MRGN_CSq pools all 1,500 trios and is",
              "the number the rest of the manuscript reports. MRGN_CSq_n300 pools only the",
              "900 trios at n <= 300, which is the like-for-like comparison against MRPC:",
              "MRPC was run only at n = 50, 150 and 300, so its 900 trios are the three",
              "smallest designs, where every method is weakest. Comparing MRPC against the",
              "full MRGN block understates MRPC."),
        paste("Confounders selected by the q-value rule -- the strictest of the three",
              "selection arms, and the one the manuscript reports. The oracle (truth) arm",
              "is not in this workbook."),
        paste("One row per block: trio count, what it pools, and the headline margins",
              "(model accuracy, no-call rate, and the edge accuracy / recall / precision",
              "from the 2x2 sheet)."),
        stamp,
        "bioinfo_revision/simulation_results/tables/confusion_counts_long.csv",
        paste("bioinfo_revision/simulation_results/results_scripts/",
              "make_aggregated_confusion_workbook.R", sep = "")),
    stringsAsFactors = FALSE)

invisible(ensure.dir(final.tables.dir))
writexl::write_xlsx(c(list(README = readme, summary = summary.tbl), sheets), out.file)

cat(sprintf("\nwrote %s (%d sheets)\n", out.file, length(sheets) + 2L))
