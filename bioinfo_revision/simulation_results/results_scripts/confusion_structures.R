# MRGN across the four confounder structures, at n = 670.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/confusion_structures.R
#
# The main simulation only ever generates one covariate structure: confounders (U) plus one
# intermediate (W) plus one common child (Z). Reviewers asked what MRGN does under the
# others, so simulation/confounder_structure_simulation.R generates three more and
# run_structure_sims.R runs MRGN over them. This script puts all four side by side.
#
#   u_only   confounders only                      data_structures/u_only/
#   u_w      confounders + 1 intermediate          data_structures/u_w/
#   u_z      confounders + 1 common child          data_structures/u_z/
#   u_w_z    confounders + both -- THE MAIN RUN    data/   (mrgn_group_n670.RData)
#
# All four are 300 trios at n = 670, 5 models x 3 effect sizes x 20 replicates, generated
# with identical effect-size strata, minor allele frequencies, U_n range and coefficient
# ranges. The ONLY thing that differs is W_n and Z_n, so a difference between the tables has
# one possible cause.
#
# MRGN ONLY, three arms (truth / CS-q / CS-alpha). The other methods are not run over these
# simulations.
#
# WHAT THE COMPARISON IS FOR. The intermediate and the common child are the two covariates
# that must NOT be adjusted for: W sits on the causal path and conditioning on it blocks the
# effect being estimated, Z is a collider and conditioning on it opens a path that was not
# there. get.conf.trios()'s filter_int_child step exists to keep both out of the selected
# set. These four tables measure what that filter is worth, one structure at a time --
# u_only is the case where there is nothing to filter, u_w and u_z isolate one hazard each,
# and u_w_z is both at once, which is the design every other table in this report uses.
#
# filter_int_child is FALSE for u_only and TRUE for the other three; see the note in
# run_structure_sims.R for why.
#
# Writes tables/confusion_structures.md and tables/structure_comparison.csv.

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

# STRUCTURE.SIZE, STRUCTURES and load.structure() now live in confusion_utils.R, so this
# script and make_structure_report.R read one definition of the directory map.

MRGN.ARMS <- c("truth", "CSq", "CSa")


# load.structure() is in confusion_utils.R.


# ---------------------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------------------

long <- list()
rates <- list()

cat("=== MRGN by confounder structure, n =", STRUCTURE.SIZE, "===\n")

for (nm in names(STRUCTURES)) {
    st <- STRUCTURES[[nm]]
    res <- load.structure(st$dir)
    if (is.null(res)) {
        cat("  ", nm, ": no mrgn_group_n", STRUCTURE.SIZE, ".RData in ", st$dir,
            " -- run run_structure_sims.R first; skipping\n", sep = "")
        next
    }
    res <- res[res$sample.size == STRUCTURE.SIZE, , drop = FALSE]
    cat("  --", nm, "--", st$label, "|", nrow(res), "trios\n")

    for (arm in MRGN.ARMS) {
        model.col <- paste0("mrgn.", arm, ".model")
        if (!model.col %in% names(res)) next

        m <- confusion(res$truth.model, res[[model.col]], MRGN.LEVELS)
        long[[length(long) + 1]] <-
            confusion.long(m, "mrgn", arm, STRUCTURE.SIZE, "all", paste0("model.", nm))

        e <- confusion(res$truth.model, mrgn.edge(res[[model.col]]),
                       MRGN.EDGE.LEVELS, coarse.pred = FALSE)
        long[[length(long) + 1]] <-
            confusion.long(e, "mrgn", arm, STRUCTURE.SIZE, "all", paste0("edge.", nm))

        correct <- sum(diag(m[, TRUTH.LEVELS, drop = FALSE]))
        es <- edge.scores(e)

        # How often selection picked up the covariate it was supposed to reject. These
        # columns are written by score.selection() for whichever set the arm used, so they
        # only mean something for the CS arms -- the truth arm's set is the U block by
        # construction and can contain neither.
        has.int <- if (arm != "truth" && paste0(arm, ".has.intermediate") %in% names(res))
            mean(res[[paste0(arm, ".has.intermediate")]], na.rm = TRUE) else NA_real_
        has.child <- if (arm != "truth" && paste0(arm, ".has.common.child") %in% names(res))
            mean(res[[paste0(arm, ".has.common.child")]], na.rm = TRUE) else NA_real_

        rates[[length(rates) + 1]] <- data.frame(
            structure = nm, structure_label = st$label, arm = arm,
            n = sum(m),
            model_accuracy = correct / sum(m),
            edge_accuracy = es$accuracy,
            edge_present_recall = es$present.recall,
            edge_absent_recall = es$absent.recall,
            no_edge_call = es$no.call,
            selected_intermediate = has.int,
            selected_common_child = has.child,
            stringsAsFactors = FALSE)

        cat(sprintf("    %-5s model %5.1f%% | edge %5.1f%% | no call %4.1f%%\n",
                    arm, 100 * correct / sum(m), 100 * es$accuracy, 100 * es$no.call))
    }
}

if (length(long) == 0) {
    stop("no structure results found. Run:\n",
         "  Rscript bioinfo_revision/simulation/confounder_structure_simulation.R\n",
         "  Rscript bioinfo_revision/simulation_results/run_structure_sims.R")
}

structure.confusion.long <- do.call(rbind, long)
structure.rates <- do.call(rbind, rates)


# ---------------------------------------------------------------------------------------
# write
# ---------------------------------------------------------------------------------------

invisible(ensure.dir(tables.dir))

csv.path <- file.path(tables.dir, "structure_comparison.csv")
utils::write.csv(cbind(structure.rates[, 1:4],
                       round(structure.rates[, -(1:4)], 4)),
                 csv.path, row.names = FALSE)
cat(sprintf("\n  wrote %s | %d rows\n", basename(csv.path), nrow(structure.rates)))

present <- unique(structure.rates$structure)

lines <- c(
    "# MRGN by confounder structure -- simulated trios, n = 670",
    "",
    paste("Generated by `results_scripts/confusion_structures.R`. Do not edit by hand;",
          "rerun the script."),
    "",
    paste("Four covariate structures, 300 trios each, identical in every other respect --",
          "same effect-size strata, same minor allele frequencies, same `U_n` range, same",
          "coefficient ranges. Only `W_n` and `Z_n` differ, so a difference between these",
          "tables has one possible cause."),
    "",
    paste("The intermediate `W` and the common child `Z` are the two covariates that must",
          "**not** be adjusted for. `W` sits on the causal path, so conditioning on it",
          "blocks the very effect being estimated; `Z` is a collider, so conditioning on it",
          "opens a path that was not there. `get.conf.trios()`'s `filter_int_child` step",
          "exists to keep both out of the selected set, and these tables measure what it is",
          "worth one hazard at a time."),
    "",
    "| case | structure | detail |",
    "| --- | --- | --- |",
    vapply(present, function(nm) sprintf("| `%s` | %s | %s |", nm,
                                         STRUCTURES[[nm]]$label, STRUCTURES[[nm]]$detail),
           character(1)),
    "",
    "## Summary rates",
    "",
    paste("`selected_intermediate` and `selected_common_child` are the share of trios whose",
          "selected set contained that trio's own `W` or `Z`. They are blank for the `truth`",
          "arm, whose set is the `U` block by construction and can contain neither, and for",
          "any structure that has no such covariate to select."),
    "")

for (arm in MRGN.ARMS) {
    sub <- structure.rates[structure.rates$arm == arm, , drop = FALSE]
    if (nrow(sub) == 0) next
    lines <- c(lines, sprintf("### %s arm", arm), "",
        "| structure | model accuracy | edge accuracy | edge-present recall | edge-absent recall | no edge call | picked W | picked Z |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |")
    fmt <- function(x) if (is.na(x)) "--" else sprintf("%.3f", x)
    for (i in seq_len(nrow(sub))) {
        lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
                                  sub$structure[i],
                                  fmt(sub$model_accuracy[i]), fmt(sub$edge_accuracy[i]),
                                  fmt(sub$edge_present_recall[i]),
                                  fmt(sub$edge_absent_recall[i]),
                                  fmt(sub$no_edge_call[i]),
                                  fmt(sub$selected_intermediate[i]),
                                  fmt(sub$selected_common_child[i])))
    }
    lines <- c(lines, "")
}

# ---- the full matrices ----
lines <- c(lines, "## Confusion matrices", "",
           paste("Columns are the generating model, rows the inferred label. Last column is",
                 "precision, bottom row recall, bottom-right cell overall accuracy."), "")

matrix.from.long <- function(counts, pred.levels) {
    m <- matrix(0L, nrow = length(TRUTH.LEVELS), ncol = length(pred.levels),
                dimnames = list(truth = TRUTH.LEVELS, predicted = pred.levels))
    m[cbind(counts$truth, counts$predicted)] <- as.integer(counts$n)
    m
}

for (nm in present) {
    lines <- c(lines, sprintf("### %s -- %s", nm, STRUCTURES[[nm]]$label), "",
               STRUCTURES[[nm]]$detail, "")
    for (arm in MRGN.ARMS) {
        sel <- structure.confusion.long[
            structure.confusion.long$arm == arm &
                structure.confusion.long$level == paste0("model.", nm), , drop = FALSE]
        if (nrow(sel) == 0) next
        m <- matrix.from.long(sel, MRGN.LEVELS)
        correct <- sum(diag(m[, TRUTH.LEVELS, drop = FALSE]))
        lines <- c(lines,
                   md.scored.table(scored.table(m, MRGN.CORRECT),
                                   sprintf("%s arm", arm),
                                   sprintf("Accuracy: **%.1f%%** (%d of %d).",
                                           100 * correct / sum(m), correct, sum(m))))
    }
}

md.path <- file.path(tables.dir, "confusion_structures.md")
writeLines(lines, md.path)
cat(sprintf("  wrote %s | %d lines\n", basename(md.path), length(lines)))
cat("\ndone.\n")
