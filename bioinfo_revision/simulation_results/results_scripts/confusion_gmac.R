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
# Writes, under simulation_results/tables/:
#   gmac/confusion_gmac.csv                       five stacked tables
#   by_effect_size/gmac/confusion_gmac_<eff>.csv  3, five stacked tables each
# and leaves the long-format counts in gmac.confusion.long for make_all_tables.R.

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")


# ---------------------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------------------

build.gmac.confusion <- function(results = load.method("gmac")) {
    long <- list()

    cat("=== GMAC confusion matrices (T1-T2 edge) ===\n")

    if (!"gmac.model" %in% names(results)) stop("results have no column 'gmac.model'")

    # guard the collapse: if gmac.model.call() ever grows a label, gmac.edge() would file
    # it under edge-present without saying so
    unknown <- setdiff(unique(as.character(results$gmac.model)), GMAC.LEVELS)
    if (length(unknown)) {
        stop("unrecognised gmac.model value(s): ", paste(unknown, collapse = ", "),
             ". gmac.edge() maps anything that is not \"No Mediation\" to edge-present; ",
             "check that is right for these before tabulating.")
    }

    pooled.blocks <- list()
    eff.blocks <- setNames(vector("list", length(EFFECT.SIZES)), EFFECT.SIZES)

    for (size in SAMPLE.SIZES) {
        rows <- results[results$sample.size == size, , drop = FALSE]

        # coarse.pred = FALSE: the edge labels are not M-labels and have no sub-type to
        # collapse. Passing them through coarse.model() would be a no-op, but saying so
        # explicitly keeps the two methods from looking interchangeable.
        m <- confusion(rows$truth.model, gmac.edge(rows$gmac.model), EDGE.LEVELS,
                       coarse.pred = FALSE)
        st <- scored.table(m, EDGE.CORRECT)

        pooled.blocks[[sprintf("n = %d", size)]] <- st
        long[[length(long) + 1]] <- confusion.long(m, "gmac", "gmac", size, "all", "edge")

        # the bottom-right cell of the scored table, re-derived rather than parsed back
        correct <- sum(vapply(TRUTH.LEVELS, function(t) m[t, EDGE.CORRECT[[t]]], numeric(1)))
        cat(sprintf("  n=%-4d  %4d trios | edge accuracy %5.1f%% | %4.1f%% called edge-present\n",
                    size, sum(m), 100 * correct / sum(m),
                    100 * sum(m[, EDGE.LEVELS[2]]) / sum(m)))

        for (eff in EFFECT.SIZES) {
            sub <- rows[rows$effect_size == eff, , drop = FALSE]
            me <- confusion(sub$truth.model, gmac.edge(sub$gmac.model), EDGE.LEVELS,
                            coarse.pred = FALSE)
            eff.blocks[[eff]][[sprintf("n = %d", size)]] <- scored.table(me, EDGE.CORRECT)
            long[[length(long) + 1]] <- confusion.long(me, "gmac", "gmac", size, eff, "edge")
        }
    }

    path <- file.path(tables.dir, "gmac", "confusion_gmac.csv")
    write.scored.csv(pooled.blocks, path)
    cat(sprintf("        -> %s (%d tables)\n", basename(path), length(pooled.blocks)))

    for (eff in EFFECT.SIZES) {
        write.scored.csv(eff.blocks[[eff]],
                         file.path(tables.dir, "by_effect_size", "gmac",
                                   sprintf("confusion_gmac_%s.csv", eff)))
    }

    do.call(rbind, long)
}

gmac.confusion.long <- build.gmac.confusion()

cat(sprintf("  wrote %d GMAC files (1 pooled, %d by effect size), %d tables each\n",
            1 + length(EFFECT.SIZES), length(EFFECT.SIZES), length(SAMPLE.SIZES)))
