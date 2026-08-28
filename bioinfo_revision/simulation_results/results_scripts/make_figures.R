# Figures 3 and 4, reproduced on the revised simulation.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/make_figures.R
#
# Writes PNGs to reports/figures/. Two figures, each in a pooled form that matches the
# published layout for direct comparison and a by-effect form that splits the three effect
# treatments:
#
#   fig3_pooled.png      selection precision and recall against a trio's true-confounder
#   fig3_by_effect.png   count, one panel per selection method. Reproduces
#                        Manuscript/figures/Figure_Confounder_Selection_Performance.pdf
#
#   fig4_pooled.png      MRGN edge-based precision and recall against the number of false
#   fig4_by_effect.png   POSITIVE confounders the selection handed it (panel a) and the
#                        number of false NEGATIVES it missed (panel b). Reproduces
#                        Manuscript/figures/SF7_Confounder_Selection_Impact_On_Inference.pdf
#                        Scored over ALL THREE trio edges (V1-T1, V1-T2, T1-T2) at the
#                        skeleton level, not on T1-T2 alone -- see the note above
#                        edge.pr.by().
#
# THIS IS THE ONE SCRIPT IN THE STAGE THAT NEEDS ggplot2. confusion_utils.R is deliberately
# base R so the tables can be rebuilt without it; the figures cannot be, so the dependency
# is isolated here rather than pushed into the shared helpers.

library(ggplot2)
library(grid)   # unit(), for the undefined-metric rug in fig4.base()

source("bioinfo_revision/simulation_results/results_scripts/selection_metrics.R")

FIG.DIR <- "bioinfo_revision/reports/figures"
dir.create(FIG.DIR, recursive = TRUE, showWarnings = FALSE)

METRIC.COLS <- c(Precision = "#2C7FB8", Recall = "#F0C808")

fig.theme <- theme_bw(base_size = 12) +
    theme(legend.position   = "top",
          legend.title      = element_text(face = "bold"),
          strip.background  = element_rect(fill = "grey85", colour = "grey30"),
          strip.text        = element_text(face = "bold"),
          panel.grid.minor  = element_blank())

save.fig <- function(p, name, width, height) {
    path <- file.path(FIG.DIR, name)
    ggsave(path, p, width = width, height = height, dpi = 150, bg = "white")
    cat(sprintf("  wrote %s (%.0f x %.0f in)\n", path, width, height))
}


# ---------------------------------------------------------------------------------------
# Figure 3: selection performance against true-confounder count
# ---------------------------------------------------------------------------------------
#
# One point per trio per metric. Trios that selected nothing have undefined precision and
# are dropped from the precision series only -- na.rm in the smooth would do this silently,
# so it is done explicitly and the count is reported to stdout. At n = 50 that is over half
# the CS-q trios, which is itself a result and is stated in the report rather than left to
# an invisible NA.
#
# The smooth is LOESS on the per-trio points, matching the published figure. Points are
# drawn under it at low alpha so the density at each confounder count stays visible.

per.trio <- selection.per.trio(load.selection())

KEEP <- c("selection.lab", "effect_size", "sample.size", "n.true.confs")
fig3.long <- rbind(
    data.frame(per.trio[, KEEP], metric = "Precision", value = per.trio$precision),
    data.frame(per.trio[, KEEP], metric = "Recall",    value = per.trio$recall))
fig3.long$selection.lab <- factor(fig3.long$selection.lab,
                                  levels = c("CS-alpha", "CS-q"),
                                  labels = c("CS-α  (α < 0.01)", "CS-q  (q < 0.05)"))
fig3.long$effect_size <- factor(fig3.long$effect_size, levels = EFFECT.SIZES)

cat("Figure 3: dropping", sum(is.na(fig3.long$value)), "undefined metric values",
    "(trios that selected nothing)\n")
fig3.long <- fig3.long[!is.na(fig3.long$value), ]

fig3.base <- function(d) {
    ggplot(d, aes(n.true.confs, value, colour = metric, fill = metric)) +
        geom_point(shape = 21, size = 1.6, alpha = 0.25, stroke = 0.2) +
        geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 1.1) +
        scale_colour_manual(values = METRIC.COLS, name = "Metric") +
        scale_fill_manual(values = METRIC.COLS, name = "Metric") +
        scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1.02)) +
        labs(x = "# of True Confounders", y = "Selection Performance") +
        fig.theme
}

# DO NOT POOL THE SAMPLE SIZES INTO ONE PANEL. The published figure is a single sample
# size. Pooling all five here mixes power regimes that differ by two orders of magnitude --
# CS-q selects nothing at all for 156 of 300 trios at n = 50 and a median of 23 covariates
# at n = 1000 -- so a pooled curve sits between two populations and describes neither, and
# the point cloud goes visibly bimodal with a band of zeros along the axis. The three
# figures below are the honest decomposition.

# (i) one sample size, the layout of the published figure, for direct comparison.
# n = 670 is the group METHODS.md section 5 compares against the published runs.
save.fig(fig3.base(fig3.long[fig3.long$sample.size == 670, ]) +
             facet_wrap(~ selection.lab, ncol = 1) +
             labs(subtitle = "n = 670"),
         "fig3_n670.png", 10, 8)

# (ii) the sample-size view: the same two panels at every n.
save.fig(fig3.base(fig3.long) +
             facet_grid(selection.lab ~ sample.size,
                        labeller = labeller(sample.size = function(x) paste0("n = ", x))),
         "fig3_by_size.png", 16, 8)

# (iii) the effect-size view. Each panel still pools all five sample sizes, so the LEVELS
# here carry the same mixing as a pooled figure would; what is readable is the comparison
# BETWEEN the three panels, since each pools the identical set of sizes.
save.fig(fig3.base(fig3.long) + facet_grid(selection.lab ~ effect_size),
         "fig3_by_effect.png", 13, 8)


# ---------------------------------------------------------------------------------------
# Figure 4: what selection errors cost MRGN
# ---------------------------------------------------------------------------------------
#
# Precision and recall are SET-level quantities -- they cannot be computed for a single
# trio -- so the x-axis is binned and the two metrics are computed within each bin, then
# smoothed. This is the only way to put them against a per-trio count, and it is what the
# published figure does.
#
# Scored on the EDGES, not the model. TRIO.EDGES (confusion_utils.R) gives the SKELETON of
# each of the eight topologies -- which of V1-T1, V1-T2 and T1-T2 are adjacent, orientation
# ignored -- and the counts are summed over the trios in a bin, so
#
#   precision = edges called that are truly there / edges called
#   recall    = edges called that are truly there / edges truly there
#
# THIS IS ALL THREE TRIO EDGES, NOT JUST T1-T2. The two-row confusion tables score T1-T2
# alone because GMAC and MR-GGI resolve nothing else and the tables exist to be read
# against them; Figure 4 is MRGN-only and has no such constraint. A trio therefore
# contributes 1-3 to each denominator rather than 0-1, and a partially-right call -- M3
# returned for an M4 trio, say -- scores 2/2 precision and 2/3 recall instead of being a
# flat miss. Orientation is not scored here (M1.1 and M2.1 share a skeleton); that is what
# the model rows of Table 2 measure.
#
# ---------------------------------------------------------------------------------------
# HOW THE X-AXIS IS BINNED, AND THE THREE BUGS THIS REPLACES
# ---------------------------------------------------------------------------------------
#
# An earlier version cut 20 equal-count bins with
# unique(quantile(x, probs = seq(0, 1, length.out = 21))) and dropped any bin holding
# fewer than 15 trios. That produced lines with holes in them, isolated points and
# truncated ends, for three separate reasons:
#
# 1. QUANTILE BINS CANNOT SPLIT A DISCRETE COUNT. Both x-axes are small integers with
#    heavy ties, and unique() then collapses the breaks. CS-q's false-positive count takes
#    3 distinct values at n = 50 (0 -> 157 trios, 1 -> 141, 2 -> 2), so 21 requested breaks
#    became 3, giving 2 bins, one of which held 2 trios and was dropped -- one point, no
#    line. Distinct values by group ran 3, 4, 5, 7, 8 against the 20 bins asked for.
#
# 2. THE FLOOR SAT EXACTLY AT THE MEAN BIN SIZE. A sample-size group is 300 trios and the
#    code asked for 20 bins: 15 per bin on average, against a MIN.BIN of 15. Occupancies
#    for CS-q false negatives at n = 150 were
#        17 13 20 11 14 16 19 10 26 8 13 16 16 13 14 17 12 20 15 10
#    so which bins survived was a coin flip, the survivors were scattered through the
#    range, and geom_line() then drew straight across every hole. The pooled figure was
#    unaffected only because 1,500 trios put 75 in a bin, comfortably clear of the floor --
#    which is why this went unnoticed.
#
# 3. THE GEOMETRY WAS NOT THE PUBLISHED ONE. SF7 is a smooth with a confidence ribbon;
#    this drew geom_line() + geom_point() over the raw bin values, so each bin's sampling
#    noise was rendered as a spike.
#
# The replacement: group by VALUE where the variable is discrete enough to allow it, derive
# the bin count from the data rather than fixing it at 20, keep every bin that is not
# degenerate, and let a WEIGHTED smooth downweight the sparse ones instead of deleting them.
# Gaps are then impossible by construction rather than by a threshold being set low enough.

# A degeneracy floor, not a power threshold. A bin of 2 trios can only take a handful of
# distinct precision values; a bin of 5 still informs a weighted smooth. Anything dropped
# here is counted and reported rather than silently removed -- see BINS.DROPPED.
MIN.BIN <- 5

# Bins are sized to hold about this many trios, and the COUNT follows from the data. This
# is the fix for bug 2: with 300 trios a group gets ~10 bins rather than 20, so occupancy
# is roughly 30 and no bin is near the floor.
TARGET.PER.BIN <- 30

# At or below this many distinct values, x is grouped BY VALUE and not binned at all. This
# is the fix for bug 1: there is nothing to bin when the variable takes three values, and
# an exact-value grouping is both simpler and lossless.
MAX.DISCRETE.LEVELS <- 15

BINS.DROPPED <- 0L   # accumulated across every panel, reported at the end

edge.pr.by <- function(df, xcol, arm) {
    model.col <- paste0("mrgn.", arm, ".model")
    truth.skel <- trio.skeleton(df$truth.model)
    pred.skel  <- trio.skeleton(df[[model.col]])

    x <- df[[xcol]]
    ok <- which(!is.na(x))
    if (length(ok) == 0) return(NULL)

    levels.x <- sort(unique(x[ok]))
    groups <- if (length(levels.x) <= MAX.DISCRETE.LEVELS) {
        lapply(levels.x, function(v) ok[x[ok] == v])
    } else {
        n.bins <- max(2L, min(20L, as.integer(floor(length(ok) / TARGET.PER.BIN))))
        br <- unique(stats::quantile(x[ok], probs = seq(0, 1, length.out = n.bins + 1),
                                     na.rm = TRUE))
        if (length(br) < 3) return(NULL)
        bin <- cut(x, breaks = br, include.lowest = TRUE)
        lapply(levels(bin), function(b) which(bin == b))
    }

    out <- lapply(groups, function(i) {
        if (length(i) < MIN.BIN) {
            BINS.DROPPED <<- BINS.DROPPED + 1L
            return(NULL)
        }
        # Edges, not trios: each trio contributes its own skeleton to the three counts,
        # so a bin of 30 trios carries roughly 60 true edges rather than 30 binary calls.
        tp <- sum(vapply(i, function(k) length(intersect(truth.skel[[k]], pred.skel[[k]])),
                         numeric(1)))
        cp <- sum(vapply(i, function(k) length(pred.skel[[k]]),  numeric(1)))
        tr <- sum(vapply(i, function(k) length(truth.skel[[k]]), numeric(1)))
        data.frame(x = mean(x[i]), trios = length(i),
                   Precision = if (cp > 0) tp / cp else NA_real_,
                   Recall    = if (tr > 0) tp / tr else NA_real_)
    })
    out <- do.call(rbind, out)
    if (is.null(out) || nrow(out) == 0) return(NULL)
    rbind(data.frame(out[, c("x", "trios")], metric = "Precision", value = out$Precision),
          data.frame(out[, c("x", "trios")], metric = "Recall",    value = out$Recall))
}

mrgn <- load.method("mrgn")
sel  <- load.selection()
stopifnot(identical(mrgn$dataset, sel$dataset))

# by = "none"   every trio in one series per arm (all five sample sizes together)
#      "effect"  split by effect treatment
#      "size"    split by sample size
#
# The same caution as Figure 3 applies: the number of false positives a selection makes is
# itself strongly a function of n, so a series pooled over sample sizes traverses several
# power regimes as it moves along its own x-axis. The "size" split is the one to read for
# mechanism; the pooled version is kept only because it is the published layout.
fig4.build <- function(by = c("none", "effect", "size")) {
    by   <- match.arg(by)
    rows <- list()
    for (arm in c("CSq", "CSa")) {
        d <- cbind(mrgn[, c("dataset", "truth.model", paste0("mrgn.", arm, ".model"))],
                   fp = sel[[paste0(arm, ".n.fp")]],
                   fn = sel[[paste0(arm, ".n.fn")]],
                   effect_size = sel$effect_size,
                   sample.size = sel$sample.size)
        groups <- switch(by, none = "all", effect = EFFECT.SIZES, size = SAMPLE.SIZES)
        col    <- switch(by, none = NA,    effect = "effect_size", size = "sample.size")
        for (g in groups) {
            dg <- if (by == "none") d else d[d[[col]] == g, , drop = FALSE]
            for (panel in c("fp", "fn")) {
                r <- edge.pr.by(dg, panel, arm)
                if (is.null(r)) next
                r$arm   <- unname(SELECTION.METHODS[arm])
                r$panel <- if (panel == "fp") "a  Number of False Positive Confounders"
                           else               "b  Number of False Negative Confounders"
                r$group <- as.character(g)
                rows[[length(rows) + 1]] <- r
            }
        }
    }
    do.call(rbind, rows)
}

# ---------------------------------------------------------------------------------------
# The smooth, fitted here rather than by geom_smooth()
# ---------------------------------------------------------------------------------------
#
# SF7 is a smooth with a confidence ribbon, so this reproduces that rather than the
# geom_line() + geom_point() it used to draw. The fit is done explicitly, for two reasons
# geom_smooth() cannot handle inside a facet:
#
#   WEIGHTS. Each point is a bin, and bins hold different numbers of trios -- a bin of 8
#   should not pull the curve as hard as one of 60. loess() takes the trio count as a
#   weight, which is also what makes the MIN.BIN floor a formality rather than the load
#   bearing decision it used to be.
#
#   PANELS THAT CANNOT SUPPORT A CURVE. loess needs several distinct x values, and CS-q's
#   false-positive panels have 3-4 at the small sample sizes. geom_smooth() would error or
#   silently drop those facets; here they fall back to points alone, which is the honest
#   rendering -- there is no curve to draw through three points and pretending otherwise is
#   what the old figure did.
MIN.X.FOR.SMOOTH <- 5

fit.smooth <- function(d) {
    xs <- sort(unique(d$x))
    if (length(xs) < MIN.X.FOR.SMOOTH) return(NULL)
    grid <- seq(min(xs), max(xs), length.out = 100)
    # degree = 1 and a wide span: these series are short and noisy, and a quadratic local
    # fit on ~10 points produces excursions outside [0, 1] that are an artefact of the
    # smoother rather than anything in the data.
    fit <- try(stats::loess(value ~ x, data = d, weights = d$trios, span = 1,
                            degree = 1, family = "gaussian",
                            control = stats::loess.control(surface = "direct")),
               silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    p <- try(stats::predict(fit, newdata = data.frame(x = grid), se = TRUE), silent = TRUE)
    if (inherits(p, "try-error") || all(is.na(p$fit))) return(NULL)
    # Clamped to [0, 1]. Precision and recall are bounded by definition, so a local fit
    # predicting 1.04 is an artefact of the smoother, not a value in the data. Clamping is
    # done here rather than left to the y-scale limits, which would DROP those vertices and
    # break the line and ribbon at exactly the points where they run to the boundary.
    data.frame(x = grid,
               value = pmin(1, pmax(0, p$fit)),
               lo = pmax(0, p$fit - 1.96 * p$se.fit),
               hi = pmin(1, p$fit + 1.96 * p$se.fit))
}

# Fits one curve per (metric x facet) cell. The facet columns are whatever the caller
# facets on, passed in by name so one helper serves all three figures.
smooth.by <- function(d, keys) {
    parts <- split(d, d[c("metric", keys)], drop = TRUE)
    out <- lapply(parts, function(p) {
        p <- p[!is.na(p$value), , drop = FALSE]
        f <- fit.smooth(p)
        if (is.null(f)) return(NULL)
        cbind(p[1, c("metric", keys), drop = FALSE], f, row.names = NULL)
    })
    out <- do.call(rbind, out)
    out
}

fig4.base <- function(d, keys) {
    sm <- smooth.by(d, keys)
    p <- ggplot(d, aes(x, value, colour = metric)) +
        scale_colour_manual(values = METRIC.COLS, name = "Metric") +
        scale_fill_manual(values = METRIC.COLS, name = "Metric") +
        scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
        labs(x = NULL, y = "MRGN Edge-Based Inference") +
        fig.theme
    if (!is.null(sm)) {
        p <- p + geom_ribbon(data = sm, aes(x, ymin = lo, ymax = hi, fill = metric),
                             inherit.aes = FALSE, alpha = 0.18) +
            geom_line(data = sm, aes(x, value, colour = metric), linewidth = 1.1)
    }
    # ---- bins where the metric is UNDEFINED, not zero ----
    #
    # Precision is 0/0 in any bin where MRGN named no edges at all -- every trio in it came
    # back "Other", whose skeleton is empty -- and that is not a precision of zero: the
    # question does not apply. Those bins carry NA, so the curve legitimately stops at the
    # last bin that had a value while recall, whose denominator is the edges the generating
    # topologies actually contain and so is never zero (every model has at least the V1-T1
    # or V1-T2 edge), carries on to the end of the panel. In the CS-alpha false-negative
    # panels that is a real result: past roughly 15 missed confounders MRGN returns "Other"
    # for everything and stops calling edges at all.
    #
    # Without a mark, that reads as a truncated line -- indistinguishable from the binning
    # bug this figure used to have. A rug at the top of the panel says where the metric
    # existed but could not be computed. It goes at the TOP deliberately: a rug along the
    # bottom would sit where y = 0 and invite exactly the "precision fell to zero" reading
    # it is there to prevent.
    undef <- d[is.na(d$value), , drop = FALSE]
    if (nrow(undef) > 0) {
        p <- p + geom_rug(data = undef, aes(x = x, colour = metric), inherit.aes = FALSE,
                          sides = "t", linewidth = 0.5, alpha = 0.8, length = unit(3, "mm"))
    }

    # The bin estimates stay visible under the curve. SF7 hides them, and hiding them here
    # is exactly what let a figure built from 2 surviving bins look like a result.
    p + geom_point(aes(size = trios), alpha = 0.55, stroke = 0) +
        scale_size_continuous(range = c(0.6, 2.6), name = "Trios in bin")
}

# facet_wrap, not facet_grid, and free x on EVERY panel. facet_grid(arm ~ panel) frees the
# x scale per column, which puts CS-alpha and CS-q on one axis -- and their false-positive
# counts differ by a factor of fifteen (54-124 against 0-7), so CS-q's whole series
# collapsed into the left margin of a 0-110 axis and read as a truncated line. The two
# arms are not on a comparable x here and should not share one.
save.fig(fig4.base(fig4.build("none"), c("arm", "panel")) +
             facet_wrap(~ panel + arm, scales = "free_x", ncol = 2),
         "fig4_pooled.png", 12, 8)

# These two keep facet_grid: their columns ARE (panel, arm) pairs, so each column already
# carries its own x range and the five rows within it are meant to share one -- that is the
# comparison the split exists to make.
f4s <- fig4.build("size")
f4s$group <- factor(f4s$group, levels = as.character(SAMPLE.SIZES),
                    labels = paste0("n = ", SAMPLE.SIZES))
save.fig(fig4.base(f4s, c("group", "arm", "panel")) +
             facet_grid(group ~ panel + arm, scales = "free_x") +
             theme(strip.text = element_text(size = 8)),
         "fig4_by_size.png", 14, 12)

f4e <- fig4.build("effect")
f4e$group <- factor(f4e$group, levels = EFFECT.SIZES)
save.fig(fig4.base(f4e, c("group", "arm", "panel")) +
             facet_grid(group ~ panel + arm, scales = "free_x") +
             theme(strip.text = element_text(size = 8)),
         "fig4_by_effect.png", 14, 9)

cat(sprintf("\nFigure 4: %d bin(s) dropped for holding fewer than %d trios\n",
            BINS.DROPPED, MIN.BIN))
cat("\ndone.\n")
