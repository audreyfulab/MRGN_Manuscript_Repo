# Regression check on the MR-GGI rewrite.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/verify_mrggi_rewrite.R
#
# MR-GGI was changed from two MRggi() calls on the bare trio to ONE call per arm with the
# arm's covariates carried as extra columns of y. Two properties must hold afterwards, and
# both are checkable against data rather than by reading the code.
#
# 1. THE `none` ARM MUST REPRODUCE THE OLD IMPLEMENTATION EXACTLY.
#    Same estimator, same inputs, so B.T1T2 and p.T1T2 must be identical trio for trio
#    against the pre-rewrite results kept in legacy/mrggi_single_arm/. This is the
#    strongest available test that the rewrite did not change what MR-GGI computes.
#
#    (MRggi rounds its output to 3 decimals internally, so the comparison is exact on the
#    rounded values it returns -- not approximate.)
#
# 2. THE ARMS MUST AGREE ON B AND p, AND DIFFER ONLY ON FDR.
#    MRggi's estimator is pairwise: .TSLS() reads only the two genes and their instruments,
#    so covariates in y cannot move the T1-T2 estimate. If they do, X has stopped lining up
#    with the columns of y in mrggi.one.trio() and every arm is suspect. The FDR column is
#    then checked to actually vary, because if it does not, the arms are carrying no
#    information at all and the whole exercise is a no-op.
#
# Exits non-zero if any check fails, so it can be used as a gate.

suppressMessages(library(MRGN))
source("bioinfo_revision/simulation_results/inference_config.R")

LEGACY <- file.path(results.root, "legacy", "mrggi_single_arm", "inference_mrggi.RData")
ARMS   <- mrggi.arms

fail <- function(...) { cat("FAIL:", ..., "\n"); quit(status = 1) }
ok   <- function(...) cat("  ok:", ..., "\n")

new <- tryCatch(loadRData(file.path(out.dir, "inference_mrggi.RData")),
                error = function(e) fail("no current MR-GGI results:", conditionMessage(e)))
cat("current results:", nrow(new), "trios,", ncol(new), "columns\n")
cat("arms:", paste(ARMS, collapse = ", "), "\n\n")


# ---------------------------------------------------------------------------------------
# 1. the `none` arm against the pre-rewrite run
# ---------------------------------------------------------------------------------------

cat("1. `none` arm vs the pre-rewrite implementation\n")
if (!"none" %in% ARMS) {
    cat("  skipped: the `none` arm is not in mrggi.arms\n\n")
} else if (!file.exists(LEGACY)) {
    cat("  skipped:", LEGACY, "not found\n",
        "  (the pre-rewrite results are only there if the old run was moved aside)\n\n")
} else {
    old <- loadRData(LEGACY)
    if (!"mrggi.B.T1T2" %in% names(old)) {
        fail("legacy file has no mrggi.B.T1T2 column; it is not a pre-rewrite run")
    }
    m <- match(old$dataset, new$dataset)
    shared <- !is.na(m)
    cat("  ", sum(shared), "trios in both\n")
    if (sum(shared) == 0) fail("no overlapping datasets to compare")

    for (q in c("B.T1T2", "p.T1T2")) {
        a <- old[[paste0("mrggi.", q)]][shared]
        b <- new[[paste0("mrggi.none.", q)]][m[shared]]
        bad <- which(!( (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b) ))
        if (length(bad) > 0) {
            cat("    ", length(bad), "mismatches on", q, "-- first few:\n")
            print(utils::head(data.frame(dataset = old$dataset[shared][bad],
                                         old = a[bad], new = b[bad]), 5))
            fail("the `none` arm does not reproduce the pre-rewrite result for", q)
        }
        ok(q, "identical on all", sum(shared), "trios")
    }

    # the weak-instrument gate reads F.T1, which the rewrite also recomputes
    if ("mrggi.F.T1" %in% names(old)) {
        a <- old$mrggi.F.T1[shared]; b <- new$mrggi.none.F.T1[m[shared]]
        d <- max(abs(a - b), na.rm = TRUE)
        if (is.finite(d) && d > 1e-8) fail("F.T1 differs by up to", d)
        ok("F.T1 identical (max abs difference", format(d, digits = 3), ")")
    }
    cat("\n")
}


# ---------------------------------------------------------------------------------------
# 2. the arms agree on the estimate and differ on the correction
# ---------------------------------------------------------------------------------------

cat("2. arm invariance of the estimate\n")
present <- ARMS[paste0("mrggi.", ARMS, ".B.T1T2") %in% names(new)]
if (length(present) < 2) {
    cat("  skipped: fewer than two arms present\n\n")
} else {
    for (q in c("B.T1T2", "p.T1T2")) {
        cols <- new[, paste0("mrggi.", present, ".", q), drop = FALSE]
        varies <- apply(cols, 1, function(x) length(unique(x[!is.na(x)])) > 1)
        if (any(varies)) {
            cat("    ", sum(varies), "trios where", q, "differs between arms -- first:\n")
            print(cols[which(varies)[1], , drop = FALSE])
            fail("covariates in y moved the T1-T2 estimate; X is misaligned with y")
        }
        ok(q, "identical across", length(present), "arms on all", nrow(new), "trios")
    }

    fdr <- new[, paste0("mrggi.", present, ".FDR.T1T2"), drop = FALSE]
    varies <- apply(fdr, 1, function(x) length(unique(x[!is.na(x)])) > 1)
    if (!any(varies)) {
        fail("FDR.T1T2 is identical across every arm on every trio. The arms are then a ",
             "no-op -- check that the covariate frames are actually reaching MRggi()")
    }
    ok("FDR.T1T2 differs on", sum(varies), "of", nrow(new),
       "trios -- the arms are doing something")

    n.cov <- sapply(present, function(a) stats::median(new[[paste0("mrggi.", a, ".n.covars")]],
                                                      na.rm = TRUE))
    cat("  median covariates per arm:",
        paste(sprintf("%s=%g", present, n.cov), collapse = "  "), "\n")
    if ("none" %in% present && n.cov[["none"]] != 0) {
        fail("the `none` arm carries covariates; it must carry none")
    }
    cat("\n")
}


# ---------------------------------------------------------------------------------------
# 3. nothing errored
# ---------------------------------------------------------------------------------------

cat("3. failures\n")
for (a in present) {
    e <- new[[paste0("mrggi.", a, ".error")]]
    n.err <- sum(!is.na(e))
    if (n.err > 0) {
        cat("    ", a, "arm:", n.err, "errored trios. First message:\n      ",
            e[!is.na(e)][1], "\n")
        fail("MR-GGI errored on some trios; investigate rather than tabulating them")
    }
}
ok("no errors in any arm")

cat("\nall checks passed\n")
