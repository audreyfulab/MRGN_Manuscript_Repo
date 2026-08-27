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
#
# THIS IS THE ONE SCRIPT IN THE STAGE THAT NEEDS ggplot2. confusion_utils.R is deliberately
# base R so the tables can be rebuilt without it; the figures cannot be, so the dependency
# is isolated here rather than pushed into the shared helpers.

library(ggplot2)

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
# Scored on the EDGE, not the model: EDGE.CORRECT maps M0 and M3 to "edge absent" and
# M1/M2/M4 to "edge present", so
#
#   precision = called present and truly present / called present
#   recall    = called present and truly present / truly present
#
# Bins holding fewer than MIN.BIN trios are dropped -- a bin of three trios produces a
# precision of 0, 1/3, 2/3 or 1 and would dominate a smooth with noise that is not signal.
MIN.BIN <- 15

edge.pr.by <- function(df, xcol, arm) {
    model.col <- paste0("mrgn.", arm, ".model")
    truth.edge <- unname(EDGE.CORRECT[coarse.model(df$truth.model)])
    pred.edge  <- mrgn.edge(df[[model.col]])

    present <- EDGE.LEVELS[2]
    x <- df[[xcol]]
    # equal-count binning: quantile breaks keep every bin populated even where the
    # distribution of x is heavily skewed, which it is for CS-alpha false positives.
    br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
    if (length(br) < 3) return(NULL)
    bin <- cut(x, breaks = br, include.lowest = TRUE)

    out <- lapply(levels(bin), function(b) {
        i <- which(bin == b)
        if (length(i) < MIN.BIN) return(NULL)
        tp <- sum(pred.edge[i] == present & truth.edge[i] == present)
        cp <- sum(pred.edge[i] == present)
        tr <- sum(truth.edge[i] == present)
        data.frame(x = mean(x[i]), trios = length(i),
                   Precision = if (cp > 0) tp / cp else NA_real_,
                   Recall    = if (tr > 0) tp / tr else NA_real_)
    })
    out <- do.call(rbind, out)
    if (is.null(out)) return(NULL)
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

fig4.base <- function(d) {
    ggplot(d, aes(x, value, colour = metric)) +
        geom_line(linewidth = 1.1) +
        geom_point(size = 1.4) +
        scale_colour_manual(values = METRIC.COLS, name = "Metric") +
        scale_y_continuous(breaks = seq(0, 1, 0.1)) +
        labs(x = NULL, y = "MRGN Edge-Based Inference") +
        fig.theme
}

save.fig(fig4.base(fig4.build("none")) + facet_grid(arm ~ panel, scales = "free_x"),
         "fig4_pooled.png", 12, 8)

f4s <- fig4.build("size")
f4s$group <- factor(f4s$group, levels = as.character(SAMPLE.SIZES),
                    labels = paste0("n = ", SAMPLE.SIZES))
save.fig(fig4.base(f4s) + facet_grid(group ~ panel + arm, scales = "free_x") +
             theme(strip.text = element_text(size = 8)),
         "fig4_by_size.png", 14, 12)

f4e <- fig4.build("effect")
f4e$group <- factor(f4e$group, levels = EFFECT.SIZES)
save.fig(fig4.base(f4e) + facet_grid(group ~ panel + arm, scales = "free_x") +
             theme(strip.text = element_text(size = 8)),
         "fig4_by_effect.png", 14, 9)

cat("\ndone.\n")
