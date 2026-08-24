# Confusion matrices for MRGN: true generating model vs inferred model, by sample size.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/confusion_mrgn.R
#
# MRGN was fitted three times per trio against three different confounder sets, and all
# three get their own matrices:
#
#   truth   the true confounders, i.e. the oracle arm -- the ceiling MRGN could reach if
#           selection were perfect. Not attainable in practice; it is the reference.
#   CSq     the CS-q selection
#   CSa     the CS-alpha selection
#
# Reading the three side by side is the point: the gap between the truth arm and CS-q/
# CS-alpha is the cost of confounder selection rather than of MRGN itself, and METHODS.md
# section 5 is about exactly that gap.
#
# Labels are coarse -- M0.1 and M0.2 are both reported as M0, using the same collapse that
# mrgn.*.correct.coarse is built on. The diagonal of each matrix is therefore the count of
# coarse-correct calls, and this script asserts that against the stored column; see the
# consistency check below.
#
# Layout: generating model across the columns, inferred model down the rows, precision in
# the last column and recall along the bottom row -- see scored.table() in
# confusion_utils.R. All five sample sizes are stacked in one file per arm.
#
# Each arm is tabulated twice, at two levels:
#
#   model   the six-way call, M0-M4 plus Other. What MRGN is actually for.
#   edge    the same trios re-scored on whether a T1-T2 edge was found, using exactly the
#           rows, columns and right answers of the GMAC table (mrgn.edge() / EDGE.CORRECT
#           in confusion_utils.R). GMAC cannot orient that edge and so cannot be scored at
#           the model level at all; this is the level the two methods share, and the only
#           one on which they can be compared head to head.
#
# Writes no files of its own. Every cell it computes -- both levels, all three arms, all
# five sample sizes, pooled and split by effect size -- is returned as long-format counts
# in mrgn.confusion.long, which make_all_tables.R stacks into
# tables/confusion_counts_long.csv and renders into tables/confusion_matrices.md.

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

MRGN.ARMS <- c("truth", "CSq", "CSa")


# ---------------------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------------------

build.mrgn.confusion <- function(results = load.method("mrgn")) {
    long <- list()

    cat("=== MRGN confusion matrices ===\n")

    for (arm in MRGN.ARMS) {
        model.col   <- paste0("mrgn.", arm, ".model")
        correct.col <- paste0("mrgn.", arm, ".correct.coarse")
        if (!model.col %in% names(results)) {
            stop("results have no column '", model.col, "'")
        }

        for (size in SAMPLE.SIZES) {
            rows <- results[results$sample.size == size, , drop = FALSE]
            m <- confusion(rows$truth.model, rows[[model.col]], MRGN.LEVELS)

            # ---- consistency check ----
            #
            # The diagonal counts trios whose coarse inferred label equals their coarse
            # truth label, which is the definition of mrgn.<arm>.correct.coarse
            # (inference_utils.R:409). The two are computed from the same two columns by
            # the same collapse, so this is an identity, not an approximation. If it fires,
            # the truth/inferred alignment is wrong and nothing else in this file should be
            # believed.
            expected <- sum(rows[[correct.col]], na.rm = TRUE)
            observed <- sum(diag(m[, TRUTH.LEVELS, drop = FALSE]))
            if (observed != expected) {
                stop(sprintf(paste0("MRGN %s n=%d: diagonal is %d but %s counts %d. ",
                                    "Truth and inferred labels are misaligned."),
                             arm, size, observed, correct.col, expected))
            }

            long[[length(long) + 1]] <- confusion.long(m, "mrgn", arm, size, "all", "model")

            e <- confusion(rows$truth.model, mrgn.edge(rows[[model.col]]),
                           MRGN.EDGE.LEVELS, coarse.pred = FALSE)
            long[[length(long) + 1]] <- confusion.long(e, "mrgn", arm, size, "all", "edge")

            es <- edge.scores(e)
            cat(sprintf("  %-5s n=%-4d  %4d trios | model accuracy %5.1f%% | edge accuracy %5.1f%% (%4.1f%% no call)\n",
                        arm, size, sum(m), 100 * observed / sum(m),
                        100 * es$accuracy, 100 * es$no.call))

            for (eff in EFFECT.SIZES) {
                sub <- rows[rows$effect_size == eff, , drop = FALSE]
                me <- confusion(sub$truth.model, sub[[model.col]], MRGN.LEVELS)
                long[[length(long) + 1]] <-
                    confusion.long(me, "mrgn", arm, size, eff, "model")

                ee <- confusion(sub$truth.model, mrgn.edge(sub[[model.col]]),
                                MRGN.EDGE.LEVELS, coarse.pred = FALSE)
                long[[length(long) + 1]] <-
                    confusion.long(ee, "mrgn", arm, size, eff, "edge")
            }
        }
    }

    do.call(rbind, long)
}

mrgn.confusion.long <- build.mrgn.confusion()
