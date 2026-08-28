# Feasibility checks behind the MR-GGI adaptation.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/mrggi_feasibility.R
#
# Kept as the permanent record for MRGGI_METHODS.md. Four questions:
#
#   1. Does MRggi()'s estimate really collapse to zero when both genes are
#      instrumented by the same SNP -- i.e. every trio?
#   2. What form must the outcome gene's "no instrument" entry take so MRggi()
#      accepts it? scale() is applied to every element of X, so NULL may not work.
#   3. Does the adapted call recover a known effect?
#   4. Does FineMapping() actually hand V1 to the trans gene, which is what
#      triggers (1)?

suppressMessages({ library(MRGN); library(MRggi); library(susieR) })

TRUE.EFFECT <- 0.7

# One synthetic trio: V1 -> T1 -> T2, with confounding on both genes.
make.trio <- function(n = 1000, seed = 1, direct.V.to.T2 = 0) {
    set.seed(seed)
    V1 <- rbinom(n, 2, 0.3)
    U  <- rnorm(n)
    T1 <- 0.8 * V1 + 0.6 * U + rnorm(n)
    T2 <- TRUE.EFFECT * T1 + direct.V.to.T2 * V1 + 0.6 * U + rnorm(n)
    list(V1 = as.matrix(V1), y = cbind(T1 = T1, T2 = T2))
}

report <- function(label, res) {
    if (inherits(res, "error")) {
        cat(sprintf("  %-42s ERROR: %s\n", label, substr(conditionMessage(res), 1, 46)))
    } else {
        cat(sprintf("  %-42s Bg1g2 = %8.4f (p %6.4f) | Bg2g1 = %8.4f (p %6.4f)\n",
                    label, res$Bg1g2[1], res$pval_Bg1g2[1],
                    res$Bg2g1[1], res$pval_Bg2g1[1]))
    }
}
try.mrggi <- function(y, X) tryCatch(MRggi::MRggi(y = y, X = X, cor.thr = 0),
                                    error = function(e) e)

tr <- make.trio()
n  <- nrow(tr$y)

cat("=== 1-2. what can the outcome gene's instrument entry be? ===\n")
cat("    true T1 -> T2 effect =", TRUE.EFFECT, "\n\n")

# the failure mode: both genes given the same instrument
report("both genes given V1  (the collapse)",
       try.mrggi(tr$y, list(T1 = tr$V1, T2 = tr$V1)))

# candidate forms for "the outcome has no instrument"
report("outcome entry = NULL",
       try.mrggi(tr$y, list(T1 = tr$V1, T2 = NULL)))
report("outcome entry = 0-column matrix",
       try.mrggi(tr$y, list(T1 = tr$V1, T2 = matrix(numeric(0), nrow = n, ncol = 0))))
report("outcome entry = column of zeros",
       try.mrggi(tr$y, list(T1 = tr$V1, T2 = matrix(0, nrow = n, ncol = 1))))
report("outcome entry = independent noise",
       try.mrggi(tr$y, list(T1 = tr$V1, T2 = matrix(rnorm(n), ncol = 1))))

cat("\n=== 3. does FineMapping hand V1 to the trans gene? ===\n")
fm <- tryCatch(MRggi::FineMapping(X = tr$V1, y = tr$y), error = function(e) e)
if (inherits(fm, "error")) {
    cat("  FineMapping errored:", conditionMessage(fm), "\n")
} else {
    for (g in names(fm)) {
        cat(sprintf("  instruments assigned to %-4s : %s\n", g,
                    if (is.null(fm[[g]]) || length(fm[[g]]) == 0) "none"
                    else paste(ncol(as.matrix(fm[[g]])), "column(s)")))
    }
    report("FineMapping's own assignment", try.mrggi(tr$y, fm))
}

cat("\n=== 4. the p.adjust bug ===\n")
cat("  MRggi() calls p.adjust(pval, method = p.adjust.methods) -- the base R\n")
cat("  CONSTANT VECTOR, not its own p.adjust.method argument. match.arg takes\n")
cat("  the first element, so the correction is always:", p.adjust.methods[1], "\n")
cat("  With one gene pair per trio this adjusts nothing, so FDR_* == pval_*.\n")


# ---------------------------------------------------------------------------------------
# 5. both directions, by swapping which gene holds the instrument
# ---------------------------------------------------------------------------------------
#
# The zeros trick gives one direction per call: whichever gene holds V1 is the exposure,
# and the other direction comes back NaN. So each trio needs two calls.

cat("\n=== 5. both directions on known truths ===\n")
zeros <- matrix(0, n, 1)
mk <- function(mod, n = 1000, seed = 1) {
    set.seed(seed)
    V1 <- rbinom(n, 2, 0.3); U <- rnorm(n)
    T1 <- 0.8*V1 + 0.6*U + rnorm(n)
    T2 <- switch(mod,
        model0 = 0.6*U + rnorm(n),                      # no edge, V1 not -> T2
        model1 = TRUE.EFFECT*T1 + 0.6*U + rnorm(n),     # T1 -> T2
        model3 = 0.8*V1 + 0.6*U + rnorm(n))             # V1 -> T2 direct, NO T1-T2 edge
    if (mod == "model2") { T2 <- 0.6*U + rnorm(n)
                           T1 <- 0.8*V1 + TRUE.EFFECT*T2 + 0.6*U + rnorm(n) }
    list(V1 = as.matrix(V1), y = cbind(T1 = T1, T2 = T2))
}
Fstat <- function(x, v) { s <- summary(lm(x ~ v)); unname(s$fstatistic[1]) }

cat(sprintf("%-8s | %-24s | %-24s\n", "model", "T1->T2  (V1 arms T1)", "T2->T1  (V1 arms T2)"))
for (mod in c("model0","model1","model2","model3")) {
    tr <- mk(mod); V <- tr$V1
    a <- suppressMessages(MRggi::MRggi(tr$y, list(T1 = V,     T2 = zeros), cor.thr = 0))
    b <- suppressMessages(MRggi::MRggi(tr$y, list(T1 = zeros, T2 = V),     cor.thr = 0))
    cat(sprintf("%-8s | B=%7.3f p=%.4f F=%5.1f | B=%8.3f p=%.4f F=%5.1f\n", mod,
                a$Bg1g2[1], a$pval_Bg1g2[1], Fstat(tr$y[,"T1"], V),
                b$Bg2g1[1], b$pval_Bg2g1[1], Fstat(tr$y[,"T2"], V)))
}
cat("truth: model0 none | model1 T1->T2 = 0.7 | model2 T2->T1 = 0.7 | model3 none\n")

# ---------------------------------------------------------------------------------------
# 6. what the p-value actually tests when there is ONE instrument
# ---------------------------------------------------------------------------------------
#
#   Bxy    = sum(Bzx*Bzy*IV) / sum(Bzx^2*IV),   se_Bxy = sqrt(1/sum(Bzx^2*IV))
#   =>  Bxy/se_Bxy = Bzy/se_Bzy * sign(Bzx)     for a single instrument
#
# The exposure's first stage cancels out of the test statistic entirely. The p-value is
# therefore the INSTRUMENT -> OUTCOME p-value and is significant whenever V1 is associated
# with the outcome gene, even when the exposure has no instrument and the reported ratio
# is a division by ~zero. This is why the reverse direction must be gated on the
# exposure's first-stage F, not on the reported p-value.
cat("\n=== 6. the p-value tests instrument -> outcome, not direction ===\n")
tr <- mk("model0")
b <- suppressMessages(MRggi::MRggi(tr$y, list(T1 = zeros, T2 = tr$V1), cor.thr = 0))
cat(sprintf("  model0, T2 -> T1 :  reported B = %.2f, p = %.3e\n", b$Bg2g1[1], b$pval_Bg2g1[1]))
cat(sprintf("  plain p for V1 -> T1 (the OUTCOME)      : %.3e   <- what the p-value tracks\n",
            coef(summary(lm(tr$y[,"T1"] ~ tr$V1)))[2,4]))
cat(sprintf("  first stage V1 -> T2 (the EXPOSURE)     : p = %.3f, F = %.1f  <- no instrument\n",
            coef(summary(lm(tr$y[,"T2"] ~ tr$V1)))[2,4], Fstat(tr$y[,"T2"], tr$V1)))
cat("\n  => report T1 -> T2 (V1 is the cis gene's eQTL by construction) and gate any\n")
cat("     reverse call on first-stage F, conventionally F > 10.\n")
