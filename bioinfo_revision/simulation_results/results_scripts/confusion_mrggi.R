# Confusion matrices for MR-GGI: true generating model vs T1-T2 edge call, by sample size.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/confusion_mrggi.R
#
# MR-GGI is scored on the T1-T2 edge, like GMAC and unlike MRGN, and for a reason internal
# to the method rather than a simplification of it. MR-GGI is a Mendelian randomisation
# estimator: it fits T2 on T1 instrumented by V1 and tests whether that causal effect is
# nonzero. Nothing in that test looks at whether V1 also acts on T2 directly, so model0 and
# model3 are the same trio to it, as are model1 and model4. It answers "is there a T1-T2
# edge", which is exactly the two-row question, and EDGE.CORRECT applies unchanged --
# M0 and M3 have no edge, M1, M2 and M4 do.
#
# ---------------------------------------------------------------------------------------
# TWO LEVELS, AND WHY THE ARMS ONLY APPEAR IN ONE OF THEM
# ---------------------------------------------------------------------------------------
#
# MR-GGI is now run with four covariate sets -- `none` (the bare trio), `truth`, `CSq` and
# `CSa` -- passed as extra columns of y. That does NOT make four different estimates.
# MRggi's estimator is strictly pairwise: .TSLS() reads only the two genes and their
# instruments, so the T1-T2 Wald ratio and its raw p-value are IDENTICAL in all four arms.
# Measured on one trio: Bg1g2 = 0.808, p = 0.000 with and without covariates in y.
#
# What the covariates change is MRggi's own multiplicity correction. It adjusts each gene's
# p-values across that gene's pairs, so with k covariates in y the T1-T2 p-value is
# corrected for T1's pairs against all of them. A wider covariate set is a harsher
# correction, and that is the entire difference between the arms.
#
# So the two levels report different things and both are needed:
#
#   edge      the raw-p call. ARM-INVARIANT by construction. This is the level that is
#             comparable with GMAC and MRGN, which get no multiplicity correction, and it
#             is tabulated ONCE rather than four identical times.
#
#   edge.fdr  the FDR-adjusted call. The only level on which the arms differ, and the one
#             that answers "what does carrying the confounders through MR-GGI cost".
#
# The arm-invariance of `edge` is ASSERTED below, not assumed. It is the strongest available
# check on the X/y alignment in mrggi.one.trio(): if X ever stopped lining up with the
# columns of y, the arms would start disagreeing on the raw call and this would stop the run.
#
# ---------------------------------------------------------------------------------------
#
# WHERE IT DIFFERS FROM GMAC: the third row. mrggi.fields() will not call an edge when the
# first-stage F for V1 -> T1 is below mrggi.min.F, and it records that trio as
# `edge = "none"`. Reading the edge column alone therefore scores an untested trio as a
# confident edge-absent call. It is not one, and there are far too many to absorb quietly:
# 570 of 1,500 trios, and 209 of 300 at n = 50, where a spurious edge-absent recall would do
# the most damage to the comparison. mrggi.edge() (confusion_utils.R) breaks them out into
# their own row, the same way MRGN's "Other" is broken out, and they count against accuracy
# -- a trio whose edge was never called is a trio whose edge was not identified.
#
# That row is also the finding, not an inconvenience. The weak-instrument rate is a property
# of the trio design rather than of the estimator: one variant, and MR-GGI needs it to be a
# strong instrument for the cis gene. It falls from 70% at n = 50 to 19% at n = 1000, so
# MR-GGI's small-sample behaviour on this simulation is governed by instrument strength
# rather than by the test. Neither MRGN nor GMAC has an analogous gate, which is what makes
# the row worth showing rather than dropping the trios.
#
# TWO COLUMNS TO READ BEFORE THE ACCURACY FIGURE, both of which are the method behaving as
# specified rather than failing. At n = 1000 every model recalls 0.65-0.78 except M2 (0.02)
# and M3 (0.03), and those two account for most of the gap between MR-GGI's accuracy and
# the others':
#
#   M3  V1 -> T1, V1 -> T2, no T1-T2 edge. Called edge-PRESENT 48 of 60 times. This is
#       exactly the exclusion-restriction violation MR cannot detect: V1 acts on T2 through
#       a path that does not go through T1, and two-stage least squares attributes that
#       association to a T1 -> T2 effect. Horizontal pleiotropy is the known failure mode
#       of single-instrument MR, and M3 is its textbook case. The column is the measurement
#       of it, not a defect in the tabulation.
#
#   M2  V1 -> T1 <- T2, edge present but oriented T2 -> T1. Called edge-ABSENT 47 of 60
#       times. Here MR-GGI is arguably right and the scoring is what is coarse: it tests
#       the cis -> trans effect, which under M2 genuinely is zero, and reports so. It is
#       penalised because EDGE.CORRECT asks only whether an edge exists, in either
#       direction.
#
# The M2 row is kept scored the same way regardless, because the alternative is worse. The
# whole purpose of the two-row framing is that MRGN, GMAC and MR-GGI are scored against
# identical rows, columns and right answers; giving MR-GGI its own definition of correct
# would make its accuracy incomparable with the other two, which is the one thing these
# tables exist to permit. GMAC's symmetric Wald test does find M2's edge, and that
# difference between the methods is a result worth being able to see. The note travels with
# the table instead -- in the report, in the header here, and it belongs in the manuscript
# text beside any MR-GGI accuracy figure.
#
# Only the cis -> trans direction is scored, and under the one-call design it is the only
# direction that exists: T2 is given no instrument, so Bg2g1 is NaN by construction. See
# MRGGI_METHODS.md and mrggi_feasibility.R sections 5-6.
#
# Layout: generating model across the columns, edge call down the rows, precision in the
# last column and recall along the bottom row. All five sample sizes are stacked in one
# file.
#
# Writes no files of its own. Every cell it computes -- both levels, all four arms, all five
# sample sizes, pooled and split by effect size -- is returned as long-format counts in
# mrggi.confusion.long, which make_all_tables.R stacks into
# tables/confusion_counts_long.csv and renders into tables/confusion_matrices.md.

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

MRGGI.ARMS <- c("none", "truth", "CSq", "CSa", "CSi")


# ---------------------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------------------

build.mrggi.confusion <- function(results = load.method("mrggi")) {
    long <- list()

    cat("=== MR-GGI confusion matrices (T1-T2 edge) ===\n")

    arms <- Filter(function(a) paste0("mrggi.", a, ".edge") %in% names(results), MRGGI.ARMS)
    if (length(arms) == 0) {
        stop("results have no mrggi.<arm>.edge columns. This script expects the four-arm ",
             "layout written by the current apply_mrggi.R; a checkpoint from the earlier ",
             "single-arm version has a bare 'mrggi.edge' column instead and needs rerunning.")
    }
    cat("  arms present:", paste(arms, collapse = ", "), "\n")

    # A trio that errored has no edge and no gate decision, and mrggi.edge() would hand
    # confusion() an NA that stops the run. There should be none -- the estimator is a
    # handful of lm() calls and nothing in it can fail on well-formed input -- so say so
    # here rather than letting the failure surface three frames down as an unrecognised
    # label.
    for (arm in arms) {
        failed <- sum(!is.na(results[[paste0("mrggi.", arm, ".error")]]))
        if (failed) {
            stop(failed, " trio(s) have a non-NA mrggi.", arm, ".error. Investigate them ",
                 "rather than tabulating them; the edge column for a failed fit is not an ",
                 "edge call.")
        }
    }

    # ---- the invariant: the raw-p call cannot depend on the arm ----
    # This is the check that the covariates are riding along in y without entering the
    # estimate, i.e. that X is still positionally aligned with the columns of y. If it ever
    # fails, mrggi.one.trio() is broken and nothing below is worth reading.
    if (length(arms) > 1) {
        base.arm <- arms[1]
        for (arm in arms[-1]) {
            disagree <- which(results[[paste0("mrggi.", base.arm, ".edge")]] !=
                              results[[paste0("mrggi.", arm, ".edge")]])
            if (length(disagree) > 0) {
                stop("the raw-p edge call differs between the '", base.arm, "' and '", arm,
                     "' arms on ", length(disagree), " trio(s) (first: dataset ",
                     results$dataset[disagree[1]], "). MRggi's estimator is pairwise, so ",
                     "covariates in y cannot move B or p -- this means X is no longer ",
                     "aligned with the columns of y in mrggi.one.trio().")
            }
        }
        cat("  raw-p edge call is identical across all", length(arms),
            "arms, as it must be\n")
    }

    # ---- level 1: the raw-p call, tabulated ONCE ----
    # Arm-invariant, so it is reported under the arm label "mrggi" -- the same label the
    # single-arm version used, which keeps this level's rows comparable with the earlier run.
    base.arm <- arms[1]
    edge.col <- paste0("mrggi.", base.arm, ".edge")
    weak.col <- paste0("mrggi.", base.arm, ".weak.instrument")

    for (size in SAMPLE.SIZES) {
        rows <- results[results$sample.size == size, , drop = FALSE]
        if (nrow(rows) == 0) next

        # coarse.pred = FALSE: the edge labels are not M-labels and have no sub-type to
        # collapse, the same as for GMAC.
        m <- confusion(rows$truth.model, mrggi.edge(rows[[edge.col]], rows[[weak.col]]),
                       MRGGI.EDGE.LEVELS, coarse.pred = FALSE)
        long[[length(long) + 1]] <- confusion.long(m, "mrggi", "mrggi", size, "all", "edge")

        es <- edge.scores(m, no.call.level = MRGGI.EDGE.LEVELS[3])
        cat(sprintf(paste("  n=%-4d  %4d trios | edge accuracy %5.1f%% |",
                          "%4.1f%% called edge-present | %4.1f%% weak instrument\n"),
                    size, sum(m), 100 * es$accuracy,
                    100 * sum(m[, EDGE.LEVELS[2]]) / sum(m), 100 * es$no.call))

        for (eff in EFFECT.SIZES) {
            sub <- rows[rows$effect_size == eff, , drop = FALSE]
            me <- confusion(sub$truth.model, mrggi.edge(sub[[edge.col]], sub[[weak.col]]),
                            MRGGI.EDGE.LEVELS, coarse.pred = FALSE)
            long[[length(long) + 1]] <- confusion.long(me, "mrggi", "mrggi", size, eff, "edge")
        }
    }

    # ---- level 2: the FDR-adjusted call, per arm ----
    # The only level on which the arms can differ. Reported at level "edge.fdr" so it never
    # gets mixed into a comparison with GMAC or MRGN, neither of which is corrected.
    cat("\n  --- FDR-adjusted call, by arm (this is where the arms differ) ---\n")
    for (arm in arms) {
        fdr.col  <- paste0("mrggi.", arm, ".edge.fdr")
        weak.col <- paste0("mrggi.", arm, ".weak.instrument")
        if (!fdr.col %in% names(results)) {
            cat("  ", arm, ": no edge.fdr column; skipping\n", sep = "")
            next
        }
        for (size in SAMPLE.SIZES) {
            rows <- results[results$sample.size == size, , drop = FALSE]
            if (nrow(rows) == 0) next

            m <- confusion(rows$truth.model, mrggi.edge(rows[[fdr.col]], rows[[weak.col]]),
                           MRGGI.EDGE.LEVELS, coarse.pred = FALSE)
            long[[length(long) + 1]] <- confusion.long(m, "mrggi", arm, size, "all", "edge.fdr")

            es <- edge.scores(m, no.call.level = MRGGI.EDGE.LEVELS[3])
            cat(sprintf(paste("  %-5s n=%-4d  %4d trios | edge accuracy %5.1f%% |",
                              "%4.1f%% called edge-present\n"),
                        arm, size, sum(m), 100 * es$accuracy,
                        100 * sum(m[, EDGE.LEVELS[2]]) / sum(m)))

            for (eff in EFFECT.SIZES) {
                sub <- rows[rows$effect_size == eff, , drop = FALSE]
                me <- confusion(sub$truth.model, mrggi.edge(sub[[fdr.col]], sub[[weak.col]]),
                                MRGGI.EDGE.LEVELS, coarse.pred = FALSE)
                long[[length(long) + 1]] <- confusion.long(me, "mrggi", arm, size, eff, "edge.fdr")
            }
        }
    }

    do.call(rbind, long)
}

mrggi.confusion.long <- build.mrggi.confusion()
