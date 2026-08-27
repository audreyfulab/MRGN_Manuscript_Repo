# Confusion matrices for MRPC: true generating model vs inferred model, by sample size.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/confusion_mrpc.R
#
# MRPC infers a trio topology and is scored the same way MRGN is: the six-way model call,
# and the same trios re-scored on the T1-T2 edge so it can be read against GMAC and MR-GGI.
# classify.mrpc.adj() (inference_utils.R) matches the fitted 3x3 adjacency against the same
# eight topologies MRGN uses, so mrgn.edge() and EDGE.CORRECT apply unchanged.
#
# ARMS, set by mrpc.arms in inference_config.R:
#
#   truth   the true confounders -- the oracle arm, the ceiling MRPC could reach if
#           selection were perfect. This is the arm to read against MRGN's truth arm.
#   CSq     the CS-q selection
# CS-alpha is NOT an MRPC arm and is not listed here. 82-106 covariates is more than MRPC
# can run a PC algorithm over at any n in this study, so the fit is excluded outright rather
# than attempted and timed out -- see mrpc.arms in inference_config.R. The mrpc.CSa.* columns
# still exist and are always NA, with the reason in mrpc.CSa.error.
#
# An arm that was not attempted has no model column values at all and is skipped here with a
# note, rather than tabulated as 300 failures -- "not run" and "ran and did not finish" are
# different facts and only the second belongs in a table.
#
# ---------------------------------------------------------------------------------------
# THE TIMEOUT COLUMN IS THE POINT, NOT AN ANNOYANCE
# ---------------------------------------------------------------------------------------
#
# apply.mrpc() caps each fit at mrpc.timeout seconds and records model = NA when it expires.
# Unlike MRGN, where an NA model is a bug to chase, for MRPC it is a measurement: at 120 s
# it was 182 of 300 trios at n = 670 and 224 of 300 at n = 1000 in the CS-q arm. The cap is
# now 180 s and the new rate is itself a result.
#
# So "Failed" is a REAL LEVEL here (MRPC.LEVELS / MRPC.EDGE.LEVELS in confusion_utils.R) and
# it counts against accuracy, the same way MRGN's "Other" and MR-GGI's "Weak instrument" do.
# A trio whose fit never returned is a trio whose model was not identified. Folding those
# into any substantive cell -- or dropping them and scoring the remainder -- would report an
# accuracy MRPC did not achieve, and would do it most at the sample sizes where the timeouts
# concentrate, which is exactly where the comparison matters.
#
# Read the "Failed" row and the accuracy figure together. They are not independent: MRPC's
# accuracy on the trios it does finish is not its accuracy on the design.
#
# Layout: generating model across the columns, inferred model down the rows, precision in
# the last column and recall along the bottom row -- see scored.table() in
# confusion_utils.R. All five sample sizes are stacked in one file per arm.
#
# Writes no files of its own. Every cell it computes -- both levels, every attempted arm,
# all five sample sizes, pooled and split by effect size -- is returned as long-format
# counts in mrpc.confusion.long, which make_all_tables.R stacks into
# tables/confusion_counts_long.csv and renders into tables/confusion_matrices.md.

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

MRPC.ARMS <- c("truth", "CSq", "CSi")


# ---------------------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------------------

build.mrpc.confusion <- function(results = load.method("mrpc")) {
    long <- list()

    cat("=== MRPC confusion matrices ===\n")

    for (arm in MRPC.ARMS) {
        model.col <- paste0("mrpc.", arm, ".model")
        error.col <- paste0("mrpc.", arm, ".error")
        timeout.col <- paste0("mrpc.", arm, ".timed.out")

        if (!model.col %in% names(results)) {
            cat("  ", arm, ": no ", model.col, " column; skipping\n", sep = "")
            next
        }

        # An arm that was not attempted still gets its columns, with the reason in the error
        # column and model NA throughout. Tabulating that would report 300 failures for a
        # fit that was never asked to run -- and "did not finish in 180 s" and "was never
        # started" are different facts, only the first of which belongs in a table.
        #
        # THIS IS CHECKED PER SAMPLE SIZE, not once per arm. mrpc.truth.max.n attempts the
        # truth arm at n <= 300 and skips it above, so a single arm is legitimately half
        # real results and half not-attempted, and an arm-level test would mislabel one half
        # whichever way it went.
        not.attempted <- function(rows) {
            all(is.na(rows[[model.col]])) && error.col %in% names(rows) &&
                any(grepl("not attempted|not in mrpc.arms|disabled in mrpc.arms",
                          as.character(rows[[error.col]])))
        }

        for (size in SAMPLE.SIZES) {
            rows <- results[results$sample.size == size, , drop = FALSE]
            if (nrow(rows) == 0) next
            if (not.attempted(rows)) {
                why <- unique(stats::na.omit(as.character(rows[[error.col]])))[1]
                cat(sprintf("  %-5s n=%-4d  not attempted: %s\n", arm, size,
                            if (is.na(why)) "reason not recorded" else why))
                next
            }

            m <- confusion(rows$truth.model, rows[[model.col]], MRPC.LEVELS)
            long[[length(long) + 1]] <- confusion.long(m, "mrpc", arm, size, "all", "model")

            e <- confusion(rows$truth.model, mrgn.edge(rows[[model.col]]),
                           MRPC.EDGE.LEVELS, coarse.pred = FALSE)
            long[[length(long) + 1]] <- confusion.long(e, "mrpc", arm, size, "all", "edge")

            # both no-call columns count against accuracy; see edge.scores()
            es <- edge.scores(e, no.call.level = c("Other", "Failed"))
            correct <- sum(diag(m[, TRUTH.LEVELS, drop = FALSE]))
            timed.out <- if (timeout.col %in% names(rows))
                sum(rows[[timeout.col]], na.rm = TRUE) else NA_integer_

            cat(sprintf(paste("  %-5s n=%-4d  %4d trios | model accuracy %5.1f%% |",
                              "edge accuracy %5.1f%% | %4.1f%% no call | %3d timed out\n"),
                        arm, size, sum(m), 100 * correct / sum(m),
                        100 * es$accuracy, 100 * es$no.call, timed.out))

            for (eff in EFFECT.SIZES) {
                sub <- rows[rows$effect_size == eff, , drop = FALSE]
                me <- confusion(sub$truth.model, sub[[model.col]], MRPC.LEVELS)
                long[[length(long) + 1]] <-
                    confusion.long(me, "mrpc", arm, size, eff, "model")

                ee <- confusion(sub$truth.model, mrgn.edge(sub[[model.col]]),
                                MRPC.EDGE.LEVELS, coarse.pred = FALSE)
                long[[length(long) + 1]] <-
                    confusion.long(ee, "mrpc", arm, size, eff, "edge")
            }
        }
    }

    if (length(long) == 0) {
        stop("no MRPC arm produced any tabulated results. Check mrpc.arms in ",
             "inference_config.R and that the inference stage has been run.")
    }
    do.call(rbind, long)
}

mrpc.confusion.long <- build.mrpc.confusion()
