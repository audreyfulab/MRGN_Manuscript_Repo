# One Excel workbook per confusion matrix, organised by method and by view.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/make_matrix_workbooks.R
#
# Writes tables/<method>/by_sample_size/*.xlsx and tables/<method>/by_effect_size/*.xlsx,
# one workbook per method-arm-group, each holding that group's matrix in three formats plus
# the scored margins. Reads tables/confusion_counts_long.csv, so run make_all_tables.R first.
#
# WHY .xlsx AND NOT .csv. make_all_tables.R removed an earlier layout of per-arm CSVs
# because they stacked several matrices into one file: "They were unreadable by read.csv() --
# a stacked file is not one rectangle." Three matrices per group cannot be one rectangle
# either, so a CSV here would repeat that mistake. Sheets remove the constraint that forced
# the stacking: each format gets a sheet, each sheet is one rectangle.
#
# THE THREE FORMATS.
#
#   model_pred_x_model_true   inferred model    x generating model
#   edge_pred_x_model_true    T1-T2 edge call   x generating model
#   edge_pred_x_edge_true     T1-T2 edge call   x true T1-T2 edge state
#
# The first two are the `model` and `edge` levels of confusion_counts_long.csv. The third is
# not stored anywhere and is derived here, by collapsing the truth axis of the second through
# EDGE.CORRECT: M0 and M3 are edge-absent, M1/M2/M4 edge-present. The predicted axis is
# untouched, so the no-call columns (Other, Failed, Weak instrument) survive into it -- they
# are not edge calls and must not be folded into either edge state.
#
# GMAC AND MR-GGI HAVE NO MODEL SHEET. GMAC returns a mediation call and MR-GGI only a T1-T2
# edge call; neither names one of the five models, so the first format does not exist for
# them. Their workbooks carry the two edge sheets, and the info sheet records why.

library(writexl)

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

METHOD.LABEL <- c(mrgn = "MRGN", mrpc = "MRPC", gmac = "GMAC", mrggi = "MR-GGI")

counts <- utils::read.csv(file.path(tables.dir, "confusion_counts_long.csv"),
                          stringsAsFactors = FALSE)

# Collapses the generating-model rows of an edge table into the two true edge states.
# Counts, so summing is exact.
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
# states, so the map is the identity. scored.table.values() needs one entry per truth level.
EDGE.IDENTITY <- stats::setNames(EDGE.LEVELS, EDGE.LEVELS)

# A scored grid as a data.frame ready for a sheet: numeric cells, row labels promoted to a
# real column because row names do not survive into a worksheet.
sheet.of <- function(m, correct.pred) {
    v <- scored.table.values(m, correct.pred)
    d <- as.data.frame(v, stringsAsFactors = FALSE, check.names = FALSE)
    cbind(inferred = rownames(v), d, stringsAsFactors = FALSE)
}

# Which levels a given method/arm actually has. MR-GGI splits: its `mrggi` arm carries the
# arm-invariant raw-p call at level `edge`, while none/truth/CSq/CSa carry the FDR-adjusted
# call at level `edge.fdr`.
levels.for <- function(method, arm) {
    lv <- unique(counts$level[counts$method == method & counts$arm == arm])
    list(model = if ("model" %in% lv) "model" else NA_character_,
         edge  = if ("edge" %in% lv) "edge" else if ("edge.fdr" %in% lv) "edge.fdr" else NA_character_)
}

stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

n.written <- 0L
per.method <- integer(0)

for (method in names(METHOD.LABEL)) {
    arms <- unique(counts$arm[counts$method == method])
    if (length(arms) == 0) next
    method.count <- 0L

    for (view in c("size", "effect")) {
        sub <- if (view == "size") "by_sample_size" else "by_effect_size"
        dir <- file.path(tables.dir, method, sub)
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)

        for (arm in arms) {
            lv <- levels.for(method, arm)
            if (is.na(lv$edge)) next          # nothing to tabulate

            edge.mats <- view.matrices(counts, method, arm, lv$edge, view,
                                       pred.levels = pred.levels.for(method, lv$edge))
            model.mats <- if (is.na(lv$model)) list() else
                view.matrices(counts, method, arm, lv$model, view,
                              pred.levels = pred.levels.for(method, lv$model))

            for (grp in names(edge.mats)) {
                em <- edge.mats[[grp]]
                sheets <- list()

                if (length(model.mats) && grp %in% names(model.mats)) {
                    sheets[["model_pred_x_model_true"]] <-
                        sheet.of(model.mats[[grp]], correct.for("model"))
                }
                sheets[["edge_pred_x_model_true"]] <- sheet.of(em, EDGE.CORRECT)
                sheets[["edge_pred_x_edge_true"]]  <-
                    sheet.of(collapse.truth.to.edge(em), EDGE.IDENTITY)

                sheets[["info"]] <- data.frame(
                    field = c("method", "arm", "view", "group", "trios",
                              "model level", "edge level",
                              "model sheet present", "generated", "source"),
                    value = c(unname(METHOD.LABEL[method]), arm,
                              if (view == "size") "by sample size" else "by effect size",
                              view.label(view, grp), as.character(sum(em)),
                              ifelse(is.na(lv$model), "none -- see below", lv$model),
                              lv$edge,
                              if (is.na(lv$model))
                                  paste("no --", unname(METHOD.LABEL[method]),
                                        "does not return a model call, so it cannot be",
                                        "cross-tabulated against the generating model")
                              else "yes",
                              stamp, "tables/confusion_counts_long.csv"),
                    stringsAsFactors = FALSE)

                fname <- sprintf("%s_%s_%s.xlsx", method, arm,
                                 if (view == "size") paste0("n", grp) else grp)
                writexl::write_xlsx(sheets, file.path(dir, fname))
                n.written <- n.written + 1L
                method.count <- method.count + 1L
            }
        }
    }
    per.method[unname(METHOD.LABEL[method])] <- method.count
    cat(sprintf("  %-7s %3d workbooks\n", METHOD.LABEL[method], method.count))
}

cat(sprintf("\nwrote %d workbooks under %s\n", n.written, tables.dir))
print(per.method)
