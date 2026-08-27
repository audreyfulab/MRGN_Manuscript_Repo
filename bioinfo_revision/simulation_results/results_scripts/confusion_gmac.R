# Confusion matrices for GMAC: true generating model vs T1-T2 edge call, by sample size.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/confusion_gmac.R
#
# GMAC is scored on the T1-T2 edge, not on the trio topology, and the two-row table is the
# point rather than a simplification. gmac.model.call() (inference_utils.R:885) returns one
# of four labels from the two mediation p-values:
#
#   Cis Mediated     cis significant, trans not
#   Trans Mediated   trans significant, cis not
#   No Mediation     neither significant
#   Undirected       both significant
#
# but GMAC's statistic is the Wald test on the mediator coefficient in
# `outcome ~ mediator + treatment + confounders`, and that regression is symmetric in T1
# and T2 -- under V1 -> T1 -> T2 the reverse coefficient is nonzero too, so the trans test
# fires whenever the cis one does. GMAC has no mechanism for orienting the edge, and the
# run bears that out: across all 1500 trios "Cis Mediated" and "Trans Mediated" hold 2 and
# 3 trios, and 1356 are "Undirected". Reporting those cells as a direction call reports
# noise.
#
# So gmac.edge() (confusion_utils.R) pools the three edge-present calls and the table asks
# the question GMAC can answer: is there a T1-T2 edge? Truth maps to an edge status via
# EDGE.CORRECT -- M0 and M3 have no T1-T2 edge, M1, M2 and M4 do -- which is what makes
# precision and recall well defined for a two-row table. This is the framing the manuscript
# already used; Manuscript/other/tablescraps/MRGN.GMAC.class.inference.50conf reports GMAC
# the same way.
#
# Layout: generating model across the columns, edge call down the rows, precision in the
# last column and recall along the bottom row. All five sample sizes are stacked in one
# file.
#
# Writes no files of its own. Every cell it computes -- all five sample sizes, pooled and
# split by effect size -- is returned as long-format counts in gmac.confusion.long, which
# make_all_tables.R stacks into tables/confusion_counts_long.csv and renders into
# tables/confusion_matrices.md.

# TWO ARMS, as for MRGN and MRPC:
#
#   gmac    the confounders GMAC selected for itself, which is the attainable result
#   truth   that trio's true U block, the oracle -- the ceiling GMAC's own selection is
#           working toward
#
# The gap between them is the cost of GMAC's confounder selection rather than of the
# mediation test, which is the same reading METHODS.md section 5 takes for MRGN. Unlike
# MRGN's arms, these two come from different code paths: the selected arm is scored from the
# batch gmac() table, because selection and testing happen together there, while the oracle
# arm bypasses selection entirely and is scored from its own apply.gmac() call.

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

# column prefix -> arm label in the tables. The selected arm keeps the bare "gmac" prefix it
# has always had, so its rows stay comparable with the earlier run.
GMAC.ARMS <- c(gmac = "gmac", gmac.truth = "truth", gmac.CSi = "CSi")


# ---------------------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------------------

build.gmac.confusion <- function(results = load.method("gmac")) {
    long <- list()

    cat("=== GMAC confusion matrices (T1-T2 edge) ===\n")

    prefixes <- names(GMAC.ARMS)[paste0(names(GMAC.ARMS), ".model") %in% names(results)]
    if (length(prefixes) == 0) {
        stop("results have no gmac model column (looked for: ",
             paste0(names(GMAC.ARMS), ".model", collapse = ", "), ")")
    }
    cat("  arms present:", paste(GMAC.ARMS[prefixes], collapse = ", "), "\n")

    for (prefix in prefixes) {
        arm <- unname(GMAC.ARMS[prefix])
        model.col <- paste0(prefix, ".model")

        # guard the collapse: if gmac.model.call() ever grows a label, gmac.edge() would
        # file it under edge-present without saying so
        unknown <- setdiff(unique(as.character(results[[model.col]])), c(GMAC.LEVELS, NA))
        if (length(unknown)) {
            stop("unrecognised ", model.col, " value(s): ", paste(unknown, collapse = ", "),
                 ". gmac.edge() maps anything that is not \"No Mediation\" to edge-present; ",
                 "check that is right for these before tabulating.")
        }

        cat("  --", arm, "arm --\n")
        for (size in SAMPLE.SIZES) {
            rows <- results[results$sample.size == size, , drop = FALSE]
            if (nrow(rows) == 0) next

            # coarse.pred = FALSE: the edge labels are not M-labels and have no sub-type to
            # collapse. Passing them through coarse.model() would be a no-op, but saying so
            # explicitly keeps the two methods from looking interchangeable.
            m <- confusion(rows$truth.model, gmac.edge(rows[[model.col]]), EDGE.LEVELS,
                           coarse.pred = FALSE)
            long[[length(long) + 1]] <- confusion.long(m, "gmac", arm, size, "all", "edge")

            # the bottom-right cell of the scored table, re-derived rather than parsed back
            correct <- sum(vapply(TRUTH.LEVELS, function(t) m[t, EDGE.CORRECT[[t]]], numeric(1)))
            cat(sprintf("  n=%-4d  %4d trios | edge accuracy %5.1f%% | %4.1f%% called edge-present\n",
                        size, sum(m), 100 * correct / sum(m),
                        100 * sum(m[, EDGE.LEVELS[2]]) / sum(m)))

            for (eff in EFFECT.SIZES) {
                sub <- rows[rows$effect_size == eff, , drop = FALSE]
                me <- confusion(sub$truth.model, gmac.edge(sub[[model.col]]), EDGE.LEVELS,
                                coarse.pred = FALSE)
                long[[length(long) + 1]] <- confusion.long(me, "gmac", arm, size, eff, "edge")
            }
        }
    }

    do.call(rbind, long)
}

gmac.confusion.long <- build.gmac.confusion()
