# Calibration check for the revised trio simulation.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
# Answers one question: do the simulated datasets carry the confounding that real Whole
# Blood trios carry? Reference values come from pc_distribution_invest/, which regresses
# each of 3,248 GTEx trios' genes on the PCs selected for that trio.
#
#   Rscript bioinfo_revision/verify_simulation.R
#   Rscript bioinfo_revision/verify_simulation.R path/to/simulated_trios.RData
#
# Nothing here feeds back into the simulation. It only reports measured against real.

library(MRGN)
library(ggplot2)
source("bioinfo_revision/simulation_utils.R")

args <- commandArgs(trailingOnly = TRUE)
SIM.FILE  <- if (length(args) > 0) args[1] else "bioinfo_revision/simulated_data/simulated_trios.RData"
POOL.FILE <- "bioinfo_revision/real_pc_effect_pools.RData"
PLOT.FILE <- "bioinfo_revision/simulated_vs_real_conf_effects.png"

# real Whole Blood reference values
REAL <- list(R2.cis = 0.412, R2.trans = 0.307,     # median R2 of a gene on its own PCs
             sd.cis = 0.117, sd.trans = 0.1035,    # sd of the per-PC correlations
             snp.cor = 0.13)                       # SNP partial correlation, modal value


# ---------------------------------------------------------------------------------------
# per-dataset measurements, all on the standardized scale
# ---------------------------------------------------------------------------------------

measure <- function(dat, params) {
    nms <- colnames(dat)
    K_n <- params$K_n; U_n <- params$U_n; n <- nrow(dat)
    U <- as.matrix(dat[, 3 + K_n + seq_len(U_n), drop = FALSE])
    V1 <- dat[, 1]; T1 <- dat[, 2]; T2 <- dat[, 3]
    w.idx <- grep("^W[0-9]+([.]|$)", nms)
    z.idx <- grep("^Z[0-9]+([.]|$)", nms)

    # adjusted R2: the raw value is inflated by U_n/n, which at n = 50 with U_n = 26 adds
    # about 0.35 whatever the population value is. That is the trap conf.r.squared() falls
    # into, so the adjusted value is the one to compare against the real 0.41 / 0.31.
    fit <- function(y) if (U_n > n - 3) NA_real_ else summary(lm(y ~ U))$adj.r.squared

    list(n = n, U_n = U_n, model = params$model, effect_size = params$effect_size,
         R2.cis = fit(T1), R2.trans = fit(T2),
         r.cis = as.vector(cor(U, T1)), r.trans = as.vector(cor(U, T2)),
         cor.V1.T1 = cor(V1, T1),
         cor.V1.W = if (length(w.idx)) cor(V1, dat[, w.idx[1]]) else NA_real_,
         cor.V1.Z = if (length(z.idx)) cor(V1, dat[, z.idx[1]]) else NA_real_,
         n.cov = ncol(dat) - 3)
}


plot.effect.comparison <- function(sim.cis, sim.trans, pools, filename = PLOT.FILE) {
    # Simulated against real confounder effects, both on the standardized scale, which is
    # the only scale the two share: sd(gene) is arbitrary in the simulation and spans 8
    # orders of magnitude across real trios.
    #
    # Densities rather than counts: the real pool has 95,564 effects against a few
    # thousand simulated ones, so raw histograms are not comparable. The two are drawn on
    # a common x range so the tails line up visually.
    df <- rbind(
        data.frame(effect = sim.cis,     source = "simulated", target = "cis gene"),
        data.frame(effect = pools$cis,   source = "real (GTEx Whole Blood)", target = "cis gene"),
        data.frame(effect = sim.trans,   source = "simulated", target = "trans gene"),
        data.frame(effect = pools$trans, source = "real (GTEx Whole Blood)", target = "trans gene"))
    df$source <- factor(df$source, levels = c("simulated", "real (GTEx Whole Blood)"))

    # per-panel sd, printed into the panel: the single number that says whether the
    # simulated confounding is the right size
    lab <- do.call(rbind, lapply(split(df, list(df$source, df$target)), function(d) {
        data.frame(target = d$target[1],
                   label = sprintf("%s: sd = %.3f, n = %s", d$source[1], sd(d$effect),
                                   format(nrow(d), big.mark = ",")),
                   source = d$source[1])
    }))
    lab$y <- ifelse(lab$source == "simulated", Inf, Inf)
    lab$vjust <- ifelse(lab$source == "simulated", 1.6, 3.0)

    p <- ggplot(df, aes(x = effect, fill = source, colour = source)) +
        geom_density(alpha = 0.35, linewidth = 0.6, adjust = 1.2) +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40", linetype = "dashed") +
        facet_wrap(~ target, ncol = 1, scales = "free_y") +
        geom_text(data = lab, aes(x = -Inf, y = y, label = label, vjust = vjust),
                  hjust = -0.05, size = 3, show.legend = FALSE, inherit.aes = FALSE,
                  colour = "grey20") +
        scale_fill_manual(values = c("simulated" = "#D95F02",
                                     "real (GTEx Whole Blood)" = "#1B9E77")) +
        scale_colour_manual(values = c("simulated" = "#D95F02",
                                       "real (GTEx Whole Blood)" = "#1B9E77")) +
        labs(title = "Confounder effects: simulated vs real",
             subtitle = paste("standardized effect b * sd(U) / sd(gene), which equals",
                              "cor(U, gene) for orthogonal confounders"),
             x = "Standardized effect (correlation)", y = "Density",
             fill = NULL, colour = NULL) +
        theme_minimal(base_size = 11) +
        theme(legend.position = "top",
              panel.grid.minor = element_blank(),
              strip.text = element_text(face = "bold"))

    ggsave(filename, plot = p, width = 7, height = 7, dpi = 150)
    cat("\nwrote", filename, "\n")
    return(invisible(p))
}


# ---------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------

cat("simulated data:", SIM.FILE, "\n")
pools <- loadRData(POOL.FILE)
sim <- loadRData(SIM.FILE)
cat("loaded", length(sim), "datasets\n")

m <- lapply(sim, function(x) measure(x$data, x$params))
sizes <- sapply(m, function(x) x$n)

# ---- 1. per-confounder effect distribution ----
#
# n = 670 only. A sample correlation carries measurement noise of about 1/sqrt(n) on top
# of the population effect, so the observed sd is sqrt(true^2 + 1/n): at n = 50 that
# inflates 0.117 to 0.184 with nothing wrong. The real pool is measured on 670 donors, so
# only the n = 670 datasets are a like-for-like comparison.
cat("\n=== per-confounder standardized effect, cor(U_i, gene), n = 670 only ===\n")
ref <- which(sizes == 670)
if (length(ref) == 0) ref <- seq_along(m)
r.cis <- unlist(lapply(m[ref], function(x) x$r.cis))
r.trans <- unlist(lapply(m[ref], function(x) x$r.trans))
cat(sprintf("  %-16s %8s %8s %8s %8s\n", "", "n", "sd", "med|r|", "max|r|"))
show.dist <- function(lbl, v) cat(sprintf("  %-16s %8d %8.4f %8.4f %8.3f\n",
                                          lbl, length(v), sd(v), median(abs(v)), max(abs(v))))
show.dist("simulated cis", r.cis)
show.dist("real cis", pools$cis)
show.dist("simulated trans", r.trans)
show.dist("real trans", pools$trans)

plot.effect.comparison(r.cis, r.trans, pools)

# ---- 2. adjusted R2 by sample size ----
cat("\n=== adjusted R2 of each gene on its own U block ===\n")
cat(sprintf("  %6s %9s %8s %9s %9s\n", "n", "datasets", "med U_n", "R2 cis", "R2 trans"))
for (s in sort(unique(sizes))) {
    i <- which(sizes == s)
    cat(sprintf("  %6d %9d %8.0f %9.3f %9.3f\n", s, length(i),
                median(sapply(m[i], function(x) x$U_n)),
                median(sapply(m[i], function(x) x$R2.cis), na.rm = TRUE),
                median(sapply(m[i], function(x) x$R2.trans), na.rm = TRUE)))
}
cat(sprintf("\n  overall median: cis %.3f (real %.3f), trans %.3f (real %.3f)\n",
            median(sapply(m, function(x) x$R2.cis), na.rm = TRUE), REAL$R2.cis,
            median(sapply(m, function(x) x$R2.trans), na.rm = TRUE), REAL$R2.trans))
cat("  Both genes are fed by one conf.coef.ranges$U, so they share a target and the\n")
cat("  real cis/trans asymmetry (0.41 vs 0.31) is not reproduced.\n")

# ---- 3. R2 against confounder count ----
cat("\n=== adjusted R2 against confounder count, n >= 300 ===\n")
cat("  real: 0.084 at 9.5 PCs, 0.279 at 18, 0.378 at 26, 0.451 at 34, 0.521 at 43\n")
big <- which(sizes >= 300)
bin <- cut(sapply(m[big], function(x) x$U_n), c(0, 10, 20, 30, 40, 50))
agg <- tapply(sapply(m[big], function(x) x$R2.cis), bin, median, na.rm = TRUE)
for (b in names(agg)) cat(sprintf("  U_n in %-9s  R2 = %.3f\n", b, agg[[b]]))

# ---- 4. SNP effect ----
cat("\n=== SNP effect ===\n")
snp <- abs(sapply(m, function(x) x$cor.V1.T1))
cat(sprintf("  median |cor(V1, T1)| = %.3f   (real partial correlation peaks at %.2f,\n",
            median(snp), REAL$snp.cor))
cat("                                  central 99% within +/-0.45)\n")
by.eff <- tapply(snp, sapply(m, function(x) x$effect_size), median)
for (e in c("small", "medium", "large")) {
    if (!is.na(by.eff[e])) cat(sprintf("    %-7s %.3f\n", e, by.eff[e]))
}

# ---- 5. can the filter see W and Z? ----
cat("\n=== intermediate / common child detectability ===\n")
cat("  get.conf.trios(filter_int_child = TRUE) correlates each trio's variant against the\n")
cat("  pooled covariates and q-values them all together, so a covariate must clear roughly\n")
cat("  filter_fdr * (true positives) / (total tests) to survive.\n\n")
cat("  Two separate obstacles, split out below. 'alpha' is the pass rate at an\n")
cat("  uncorrected alpha = 0.05, which isolates raw power; 'pooled' is the pass rate at\n")
cat("  the actual q-value threshold, so the gap between them is the multiplicity cost.\n")
cat("  'true r' de-noises the observed correlation, sqrt(mean(r^2) - 1/n): a sample\n")
cat("  correlation carries about 1/sqrt(n) of noise, which is why the observed medians\n")
cat("  FALL with n while the underlying effect is constant.\n\n")
cat(sprintf("  %6s | %8s %8s | %8s %8s | %10s %7s %7s | %7s %7s\n",
            "n", "med|rVW|", "true rW", "med|rVZ|", "true rZ",
            "threshold", "W pool", "Z pool", "W alpha", "Z alpha"))
for (s in sort(unique(sizes))) {
    i <- which(sizes == s)
    thr <- 0.1 * (2 * length(i)) / (length(i) * sum(sapply(m[i], function(x) x$n.cov)))
    pv <- function(r) { tv <- r * sqrt((s - 2) / (1 - r^2)); 2 * pt(-abs(tv), s - 2) }
    # E[r_obs^2] ~ r_true^2 + 1/n, so subtracting the noise floor recovers the effect
    denoise <- function(r) sqrt(max(0, mean(r^2, na.rm = TRUE) - 1 / s))
    rw <- sapply(m[i], function(x) x$cor.V1.W)
    rz <- sapply(m[i], function(x) x$cor.V1.Z)
    cat(sprintf("  %6d | %8.3f %8.3f | %8.3f %8.3f | %10.2e %7.3f %7.3f | %7.3f %7.3f\n",
                s, median(abs(rw), na.rm = TRUE), denoise(rw),
                median(abs(rz), na.rm = TRUE), denoise(rz), thr,
                mean(pv(rw) < thr, na.rm = TRUE), mean(pv(rz) < thr, na.rm = TRUE),
                mean(pv(rw) < 0.05, na.rm = TRUE), mean(pv(rz) < 0.05, na.rm = TRUE)))
}
# sample size needed to clear each hurdle for the median trio, from n = (z/r)^2
cat("\n  n needed for the median trio (n = (z/r)^2, using the de-noised r at n = 1000):\n")
big <- which(sizes == max(sizes))
s.big <- max(sizes)
for (blk in c("W", "Z")) {
    r <- sapply(m[big], function(x) if (blk == "W") x$cor.V1.W else x$cor.V1.Z)
    r.true <- sqrt(max(0, mean(r^2, na.rm = TRUE) - 1 / s.big))
    thr <- 0.1 * (2 * length(big)) / (length(big) * sum(sapply(m[big], function(x) x$n.cov)))
    cat(sprintf("    %s: true r = %.3f  ->  n ~ %s at alpha = 0.05,  n ~ %s pooled\n",
                blk, r.true,
                format(round((qnorm(0.975) / r.true)^2), big.mark = ","),
                format(round((qnorm(1 - thr / 2) / r.true)^2), big.mark = ",")))
}
cat("\n  |cor(V1, W)| and |cor(V1, Z)| are bounded above by |cor(V1, T1)|, since W and Z\n")
cat("  are downstream of T1. That ceiling is set by the real cis-eQTL effect (~0.13), so\n")
cat("  detecting these blocks through their correlation with the variant is intrinsically\n")
cat("  hard and cannot be fixed by recalibrating the confounders.\n")
