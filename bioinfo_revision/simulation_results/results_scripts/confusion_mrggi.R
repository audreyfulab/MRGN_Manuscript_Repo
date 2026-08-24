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
# WHERE IT DIFFERS FROM GMAC: the third row. mrggi.fields() (inference_utils.R:1253) will
# not call an edge when the first-stage F for V1 -> T1 is below mrggi.min.F, and it records
# that trio as `edge = "none"`. Reading the edge column alone therefore scores an untested
# trio as a confident edge-absent call. It is not one, and there are far too many to absorb
# quietly: 570 of 1,500 trios, and 209 of 300 at n = 50, where a spurious edge-absent
# recall would do the most damage to the comparison. mrggi.edge() (confusion_utils.R)
# breaks them out into their own row, the same way MRGN's "Other" is broken out, and they
# count against accuracy -- a trio whose edge was never called is a trio whose edge was not
# identified.
#
# That row is also the finding, not an inconvenience. The weak-instrument rate is a
# property of the trio design rather than of the estimator: one variant, and MR-GGI needs
# it to be a strong instrument for the cis gene. It falls from 70% at n = 50 to 19% at
# n = 1000, so MR-GGI's small-sample behaviour on this simulation is governed by instrument
# strength rather than by the test. Neither MRGN nor GMAC has an analogous gate, which is
# what makes the row worth showing rather than dropping the trios.
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
# Only the cis -> trans direction is scored. The reverse estimate is computed and stored in
# mrggi.B.T2T1/mrggi.p.T2T1 but is not an edge call: with a single instrument MRggi's test
# statistic reduces to the instrument -> outcome t-statistic and the exposure's first stage
# cancels out, so the reverse direction reports significance on trios where V1 is not an
# instrument for T2 at all. See MRGGI_METHODS.md and mrggi_feasibility.R sections 5-6.
#
# Layout: generating model across the columns, edge call down the rows, precision in the
# last column and recall along the bottom row. All five sample sizes are stacked in one
# file.
#
# Writes no files of its own. Every cell it computes -- all five sample sizes, pooled and
# split by effect size -- is returned as long-format counts in mrggi.confusion.long, which
# make_all_tables.R stacks into tables/confusion_counts_long.csv and renders into
# tables/confusion_matrices.md.

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")


# ---------------------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------------------

build.mrggi.confusion <- function(results = load.method("mrggi")) {
    long <- list()

    cat("=== MR-GGI confusion matrices (T1-T2 edge) ===\n")

    for (col in c("mrggi.edge", "mrggi.weak.instrument")) {
        if (!col %in% names(results)) stop("results have no column '", col, "'")
    }

    # A trio that errored has no edge and no gate decision, and mrggi.edge() would hand
    # confusion() an NA that stops the run. There are none -- the estimator is a handful of
    # lm() calls and nothing in it can fail on well-formed input -- so say so here rather
    # than letting the failure surface three frames down as an unrecognised label.
    failed <- sum(!is.na(results$mrggi.error))
    if (failed) {
        stop(failed, " trio(s) have a non-NA mrggi.error. Investigate them rather than ",
             "tabulating them; the edge column for a failed fit is not an edge call.")
    }

    for (size in SAMPLE.SIZES) {
        rows <- results[results$sample.size == size, , drop = FALSE]

        # coarse.pred = FALSE: the edge labels are not M-labels and have no sub-type to
        # collapse, the same as for GMAC.
        m <- confusion(rows$truth.model,
                       mrggi.edge(rows$mrggi.edge, rows$mrggi.weak.instrument),
                       MRGGI.EDGE.LEVELS, coarse.pred = FALSE)
        long[[length(long) + 1]] <- confusion.long(m, "mrggi", "mrggi", size, "all", "edge")

        es <- edge.scores(m, no.call.level = MRGGI.EDGE.LEVELS[3])
        cat(sprintf(paste("  n=%-4d  %4d trios | edge accuracy %5.1f%% |",
                          "%4.1f%% called edge-present | %4.1f%% weak instrument\n"),
                    size, sum(m), 100 * es$accuracy,
                    100 * sum(m[, EDGE.LEVELS[2]]) / sum(m), 100 * es$no.call))

        for (eff in EFFECT.SIZES) {
            sub <- rows[rows$effect_size == eff, , drop = FALSE]
            me <- confusion(sub$truth.model,
                            mrggi.edge(sub$mrggi.edge, sub$mrggi.weak.instrument),
                            MRGGI.EDGE.LEVELS, coarse.pred = FALSE)
            long[[length(long) + 1]] <- confusion.long(me, "mrggi", "mrggi", size, eff, "edge")
        }
    }

    do.call(rbind, long)
}

mrggi.confusion.long <- build.mrggi.confusion()
