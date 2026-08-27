# Edge precision and recall for all four methods, in a sample-size x effect-treatment grid.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/make_edge_pr_figure.R
#
# Standalone. Reads ../tables/confusion_counts_long.csv -- so run make_all_tables.R first --
# and writes:
#
#   ../../reports/figures/fig_edge_pr_grid.png    the figure
#   ../tables/edge_pr_grid.csv                    the numbers behind it
#
# The figure is referenced from ../../reports/INFERENCE_PERFORMANCE.md, so rerun
# make_inference_report.R after this if the caption numbers need to move with it.
#
# ---------------------------------------------------------------------------------------
# WHY THE EDGE, AND WHY ONE ARM PER METHOD
# ---------------------------------------------------------------------------------------
#
# The four methods do not answer the same question. MRGN and MRPC name a model; GMAC
# returns a mediation call and MR-GGI a pairwise 2SLS edge test, and neither can name one.
# The only question all four answer is whether the trio has a T1-T2 edge, which is what
# `level = "edge"` in confusion_counts_long.csv scores. Comparing them at the model level
# is not possible, and this figure does not try.
#
# One arm each, the one the method would actually be run under:
#
#   MRGN    CSq     the CS-q selection
#   MRPC    CSq     ditto. Its oracle arm is a handicap rather than a ceiling -- see
#                   INFERENCE_PERFORMANCE.md section 3.1 -- so CS-q is also its best arm.
#   MR-GGI  mrggi   the raw-p edge call, which is ARM-INVARIANT by construction (the
#                   estimator is strictly pairwise; the arms differ only in MRggi()'s own
#                   multiplicity correction, i.e. at the edge.fdr level). So there is
#                   exactly one MR-GGI edge series, not four.
#   GMAC    gmac    its own selection; it has no other arm in this run.
#
# Precision and recall are for the EDGE-PRESENT class, matching edge.scores()
# (confusion_utils.R:380) and the `Edge precision`/`Edge recall` columns of the section-2
# scorecards. Trios that got no edge call at all -- MRGN's "Other", MRPC's "Other"/"Failed",
# MR-GGI's "Weak instrument" -- are in neither the numerator nor the denominator of
# precision, but they ARE in the denominator of recall. That asymmetry is the point: a
# method that declines to call an edge is not punished on precision and is punished on
# recall, which is exactly how MR-GGI ends up at 63% precision and 49% recall at n = 1000.
#
# ---------------------------------------------------------------------------------------
# WHAT THE BOX IS -- READ THIS BEFORE INTERPRETING THE SPREAD
# ---------------------------------------------------------------------------------------
#
# Precision and recall are SET-level quantities. One trio does not have a precision, so
# there is no per-trio distribution to draw a box over, and the boxes here are NOT spread
# across trios.
#
# Each box is a STRATIFIED NONPARAMETRIC BOOTSTRAP of the 100 trios in that cell: it shows
# how precisely the rate is pinned down by 100 trios, i.e. a confidence interval drawn as a
# box. The observed rate for the cell is marked separately as a white point, so the box can
# never be mistaken for the estimate itself.
#
# NOTE that a panel here is a FINER split than either scorecard view. The section-2 tables
# pool one axis -- the sample-size view sums the three effect treatments, the effect-size
# view sums the five sample sizes -- so these cells aggregate to those rows and do not
# match them one for one. Verified: the three per-effect cells at n = 1000 for GMAC sum to
# exactly the 67.2% / 93.3% that the sample-size view reports.
#
# The bootstrap is STRATIFIED BY TRUTH MODEL, resampling the 20 trios within each of the
# five generating models rather than 100 trios at large. The design fixes 20 replicates per
# model per cell, so the class balance (60 edge-present, 40 edge-absent) is a property of
# the design and not something that varied between runs; an unstratified bootstrap would
# jitter the recall denominator and report design constants as sampling noise.
#
# The obvious alternative -- one point per replicate -- is degenerate here. A replicate is
# five trios, one per model, three of them edge-present, so its recall can only be 0, 1/3,
# 2/3 or 1. Twenty such values are not a distribution worth drawing.

library(ggplot2)

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

FIG.DIR <- "bioinfo_revision/reports/figures"
dir.create(FIG.DIR, recursive = TRUE, showWarnings = FALSE)

N.BOOT <- 2000
set.seed(20260827)   # fixed: the figure must not move between runs of the same data

# method, arm and the label to print. Order is the x-axis order, and it is the same
# fastest-to-slowest ordering used by compute_time.R so the two figures line up.
SERIES <- list(
    list(method = "mrgn",  arm = "CSq",   label = "MRGN"),
    list(method = "mrpc",  arm = "CSq",   label = "MRPC"),
    list(method = "mrggi", arm = "mrggi", label = "MR-GGI"),
    list(method = "gmac",  arm = "gmac",  label = "GMAC"))

METHOD.ORDER <- vapply(SERIES, `[[`, "", "label")

# The house precision/recall pair from make_figures.R, with the yellow darkened. The
# original #F0C808 sits at L 0.84 against a light surface and clears only 1.6:1 contrast,
# which is fine behind the points and lines of Figures 3-4 but too weak as a box fill.
# #B8860B is the same hue at a passing lightness:
#   node scripts/validate_palette.js "#2C7FB8,#B8860B" --mode light   -> all checks pass
METRIC.COLS <- c(Precision = "#2C7FB8", Recall = "#B8860B")

PRESENT <- EDGE.LEVELS[2]
TRUTH.PRESENT <- names(EDGE.CORRECT)[EDGE.CORRECT == PRESENT]


# ---------------------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------------------

# Same arithmetic as edge.scores(), on a truth x predicted matrix. Kept separate because
# edge.scores() returns all six rates and this is called 2000 times per cell; and because
# the bootstrap needs it to accept a resampled matrix rather than a data frame.
pr <- function(m) {
    present.col <- if (PRESENT %in% colnames(m)) m[, PRESENT] else
        stats::setNames(rep(0, nrow(m)), rownames(m))
    tp     <- sum(present.col[intersect(TRUTH.PRESENT, rownames(m))])
    called <- sum(present.col)
    true   <- sum(m[intersect(TRUTH.PRESENT, rownames(m)), , drop = FALSE])
    c(Precision = if (called > 0) tp / called else NA_real_,
      Recall    = if (true   > 0) tp / true   else NA_real_)
}

# One cell of the grid, as a truth x predicted matrix of counts.
cell.matrix <- function(counts, method, arm, n, effect) {
    d <- counts[counts$method == method & counts$arm == arm & counts$level == "edge" &
                    counts$sample_size == n & counts$effect_size == effect, ]
    if (nrow(d) == 0) return(NULL)
    m <- tapply(d$n, list(d$truth, d$predicted), sum)
    m[is.na(m)] <- 0
    m
}

# Stratified bootstrap: each truth row is resampled to its own total, independently, so the
# row margins -- the 20 trios per generating model that the design fixes -- are held.
# rmultinom is the resample: drawing `row total` labels with replacement from a row of
# counts is exactly a multinomial draw on that row's proportions.
boot.pr <- function(m, B = N.BOOT) {
    rows <- rownames(m)
    draws <- lapply(rows, function(r) {
        total <- sum(m[r, ])
        if (total == 0) return(matrix(0, nrow = ncol(m), ncol = B))
        stats::rmultinom(B, total, m[r, ])
    })
    out <- vapply(seq_len(B), function(b) {
        mb <- t(vapply(draws, function(d) d[, b], numeric(ncol(m))))
        dimnames(mb) <- dimnames(m)
        pr(mb)
    }, numeric(2))
    data.frame(replicate = seq_len(B), Precision = out["Precision", ],
               Recall = out["Recall", ])
}


# ---------------------------------------------------------------------------------------
# building the grid
# ---------------------------------------------------------------------------------------

counts <- utils::read.csv(file.path(tables.dir, "confusion_counts_long.csv"),
                          stringsAsFactors = FALSE)

cat("building", length(SERIES), "series x", length(SAMPLE.SIZES), "sizes x",
    length(EFFECT.SIZES), "effects...\n")

boot.rows <- list()
obs.rows  <- list()
missing   <- list()

for (s in SERIES) {
    for (n in SAMPLE.SIZES) {
        for (eff in EFFECT.SIZES) {
            m <- cell.matrix(counts, s$method, s$arm, n, eff)
            if (is.null(m)) {
                # MRPC has no n = 670 or n = 1000 group. That is an absent run, not a zero
                # score, and it is recorded so the panel can say so rather than quietly
                # dropping a box and leaving three where there should be four.
                missing[[length(missing) + 1]] <-
                    data.frame(method = s$label, sample.size = n, effect_size = eff)
                next
            }
            stopifnot(sum(m) == 100)   # 5 models x 20 replicates; a short cell is a bug

            b <- boot.pr(m)
            boot.rows[[length(boot.rows) + 1]] <- data.frame(
                method = s$label, sample.size = n, effect_size = eff, b)

            # n.called is the DENOMINATOR of precision, and it is carried through to the
            # CSV because without it a precision cannot be read. MRPC at n = 150, small
            # effect scores 1.00 -- on two calls. That is not the same claim as GMAC's
            # 0.63 on ninety-odd, and a box alone cannot tell them apart: the bootstrap
            # of a two-call cell collapses onto 1.00 and looks like certainty.
            o <- pr(m)
            obs.rows[[length(obs.rows) + 1]] <- data.frame(
                method = s$label, sample.size = n, effect_size = eff,
                Precision = o[["Precision"]], Recall = o[["Recall"]],
                n.called.present = if (PRESENT %in% colnames(m)) sum(m[, PRESENT]) else 0L,
                n.true.present = sum(m[intersect(TRUTH.PRESENT, rownames(m)), ]))
        }
    }
    cat("  ", s$label, "\n")
}

boot <- do.call(rbind, boot.rows)
obs  <- do.call(rbind, obs.rows)
gaps <- if (length(missing)) do.call(rbind, missing) else NULL

# Wide -> long, once, for both frames.
to.long <- function(d) {
    keys <- c("method", "sample.size", "effect_size")
    rbind(data.frame(d[, keys], metric = "Precision", value = d$Precision),
          data.frame(d[, keys], metric = "Recall",    value = d$Recall))
}
boot.long <- to.long(boot)
obs.long  <- to.long(obs)

# Precision is undefined when a resample called no edge present at all. That is a real
# property of the method in that cell -- it happens where a method declines to call almost
# everything -- so it is reported rather than dropped silently, and the box is drawn from
# the resamples where it IS defined.
undef <- boot.long[boot.long$metric == "Precision" & is.na(boot.long$value), ]
if (nrow(undef)) {
    tb <- table(paste(undef$method, undef$sample.size, undef$effect_size))
    cat(sprintf("\nprecision undefined in %d of %d resamples (no edge called present):\n",
                nrow(undef), N.BOOT * length(unique(paste(boot$method, boot$sample.size,
                                                          boot$effect_size)))))
    for (k in names(tb)) cat(sprintf("    %-28s %4d / %d\n", k, tb[[k]], N.BOOT))
}
boot.long <- boot.long[!is.na(boot.long$value), ]

as.factors <- function(d) {
    d$method      <- factor(d$method, levels = METHOD.ORDER)
    d$metric      <- factor(d$metric, levels = names(METRIC.COLS))
    d$size.lab    <- factor(paste0("n = ", d$sample.size),
                            levels = paste0("n = ", SAMPLE.SIZES))
    d$effect.lab  <- factor(paste0(d$effect_size, " effect"),
                            levels = paste0(EFFECT.SIZES, " effect"))
    d
}
boot.long <- as.factors(boot.long)
# Same treatment for the observed points: a 0/0 precision has no point to draw. Dropped
# here rather than left for ggplot to discard with a warning, and annotated on the panel
# below so the blank slot reads as "the question does not apply" and not as an oversight.
obs.long  <- as.factors(obs.long[!is.na(obs.long$value), ])


# ---------------------------------------------------------------------------------------
# the figure
# ---------------------------------------------------------------------------------------
#
# Rows are the effect treatment, columns the sample size, so reading ACROSS a row is "what
# does more data buy" and reading DOWN a column is "what does a stronger effect buy". Both
# axes are ordered, so the grid is a surface, not a set of unrelated panels.
#
# Precision and recall are dodged side by side within each method rather than given their
# own facet row. The comparison that matters is between the two -- a method with high
# precision and low recall is doing something quite different from one with the reverse,
# and MR-GGI vs GMAC is exactly that contrast -- and splitting them across rows puts the
# two halves of it a panel apart.
#
# The observed value is drawn as a white point on top of each box. It is what the section-2
# scorecards report; the box around it is the bootstrap. Having both in the figure is what
# stops the box being read as spread across trios.

DODGE <- position_dodge(width = 0.78)

p <- ggplot(boot.long, aes(method, value, fill = metric)) +
    geom_boxplot(width = 0.7, position = DODGE, outlier.shape = NA,
                 linewidth = 0.35, colour = "grey25", alpha = 0.85) +
    geom_point(data = obs.long, position = DODGE, shape = 21, size = 1.5,
               fill = "white", colour = "grey15", stroke = 0.45) +
    scale_fill_manual(values = METRIC.COLS, name = "Metric", drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    # No hard limits. Both rates live in [0, 1] by construction, so a limit adds nothing,
    # and limits = c(0, 1) clips the top of a box sitting exactly at 1.00 -- which several
    # do, since a method that calls three edges and gets all three right has precision 1.
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = c(0.04, 0.04))) +
    facet_grid(effect.lab ~ size.lab) +
    labs(x = NULL, y = "Edge precision / recall",
         title = "Edge precision and recall, four methods, by sample size and effect treatment",
         subtitle = paste0(
             "Each box is a stratified bootstrap of that cell's 100 trios -- how precisely ",
             "100 trios pin the rate down, NOT spread across trios.\nThe white point is the ",
             "observed value for that cell -- section 2's tables pool one axis, these ",
             "panels are the cross, so cells aggregate to those rows.\nScored on ",
             "the T1-T2 edge, the only question all four answer; MRGN / MRPC on CS-q, ",
             "GMAC on its own selection, MR-GGI arm-invariant at this level.  ",
             "† = precision from fewer than 10 edge calls; read it with the count, ",
             "not on its own.")) +
    theme_bw(base_size = 11) +
    theme(legend.position    = "top",
          legend.title       = element_text(face = "bold"),
          strip.background   = element_rect(fill = "grey85", colour = "grey30"),
          strip.text         = element_text(face = "bold"),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.subtitle      = element_text(colour = "grey30", size = 8.5),
          axis.text.x        = element_text(angle = 45, hjust = 1))

# Two kinds of empty slot, and neither is allowed to read as "scored zero".
#
#   not run       MRPC has no n = 670 / n = 1000 group at all.
#   no edge call  the method ran and called ZERO trios edge-present in that cell, so its
#                 precision is 0/0. That is not a precision of zero -- the question does
#                 not apply -- and it is the same blank-not-zero rule the confusion tables
#                 use for a label nothing was assigned to. Its RECALL is still defined and
#                 still drawn: it is 0, and correctly so.
annotate.gap <- function(p, d, label) {
    if (is.null(d) || nrow(d) == 0) return(p)
    d <- as.factors(cbind(d, metric = "Precision", value = 0.5))
    p + geom_text(data = d, aes(method, 0.5, label = label), inherit.aes = FALSE,
                  angle = 90, size = 2.6, colour = "grey45", fontface = "italic")
}

no.edge <- unique(obs[is.na(obs$Precision), c("method", "sample.size", "effect_size")])
p <- annotate.gap(p, gaps, "not run")
p <- annotate.gap(p, no.edge, "no edge called")

# A precision computed from a handful of calls is not a precise precision, and the box
# cannot say so: bootstrap a cell where 2 of 2 calls were right and almost every resample
# returns 1.00, so it collapses onto a flat line at the top of the panel and reads as
# certainty. Those cells are daggered, and the count is in edge_pr_grid.csv.
THIN.CALLS <- 10
thin <- obs[!is.na(obs$Precision) & obs$n.called.present < THIN.CALLS, ]
if (nrow(thin)) {
    # Above the box, not above the observed value: several of these cells have a wide
    # bootstrap, so an offset from the point lands inside the box it is meant to flag.
    top <- vapply(seq_len(nrow(thin)), function(i) {
        v <- boot.long$value[boot.long$metric == "Precision" &
                                 boot.long$method == thin$method[i] &
                                 boot.long$sample.size == thin$sample.size[i] &
                                 boot.long$effect_size == thin$effect_size[i]]
        min(1.06, max(v) + 0.05)
    }, numeric(1))
    thin <- as.factors(cbind(thin[, c("method", "sample.size", "effect_size")],
                             metric = "Precision", value = top))
    p <- p + geom_text(data = thin, aes(method, value, label = "†"),
                       position = DODGE, size = 3.4, colour = "grey20")
    cat(sprintf("\n%d cells have precision from fewer than %d edge calls (daggered)\n",
                nrow(thin), THIN.CALLS))
}
if (nrow(no.edge)) {
    cat("\ncells where the method called no edge present (precision undefined):\n")
    for (i in seq_len(nrow(no.edge))) {
        cat(sprintf("    %-8s n = %-5d %s effect\n", no.edge$method[i],
                    no.edge$sample.size[i], no.edge$effect_size[i]))
    }
}

path <- file.path(FIG.DIR, "fig_edge_pr_grid.png")
ggsave(path, p, width = 15, height = 9.5, dpi = 150, bg = "white")
cat(sprintf("\n  wrote %s\n", path))

# The numbers behind the figure: the observed rate and the bootstrap interval for each box.
summ <- do.call(rbind, lapply(split(boot.long, list(
    boot.long$method, boot.long$sample.size, boot.long$effect_size, boot.long$metric),
    drop = TRUE), function(d) {
        row <- obs[obs$method == as.character(d$method[1]) &
                       obs$sample.size == d$sample.size[1] &
                       obs$effect_size == as.character(d$effect_size[1]), ]
        metric <- as.character(d$metric[1])
        q <- stats::quantile(d$value, c(0.025, 0.25, 0.5, 0.75, 0.975), names = FALSE)
        data.frame(method = as.character(d$method[1]), sample_size = d$sample.size[1],
                   effect_size = as.character(d$effect_size[1]), metric = metric,
                   observed = row[[metric]],
                   # the denominator, without which the rate above cannot be read:
                   # trios called edge-present for precision, trios truly edge-present
                   # (60 by design in every cell) for recall
                   denominator = if (metric == "Precision") row$n.called.present
                                 else row$n.true.present,
                   boot_lo95 = q[1], boot_q1 = q[2], boot_median = q[3],
                   boot_q3 = q[4], boot_hi95 = q[5], n_resamples = nrow(d),
                   stringsAsFactors = FALSE)
    }))
rownames(summ) <- NULL
summ <- summ[order(match(summ$method, METHOD.ORDER), summ$sample_size,
                   match(summ$effect_size, EFFECT.SIZES), summ$metric), ]
write.csv(summ, file.path(tables.dir, "edge_pr_grid.csv"), row.names = FALSE)
cat(sprintf("  wrote %s (%d rows)\n", file.path(tables.dir, "edge_pr_grid.csv"),
            nrow(summ)))

cat("\ndone.\n")
