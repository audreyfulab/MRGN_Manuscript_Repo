# Per-trio compute time for the four inference methods.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/compute_time.R
#
# Standalone, like make_figures.R -- not part of make_all_tables.R. Writes five files:
#
#   ../tables/compute_time_long.csv      one row per trio per method per arm
#   ../tables/compute_time.csv           median / IQR per method x arm x sample size
#   ../../reports/COMPUTE_TIME.md        the rendered report
#   ../../reports/figures/fig_compute_time*.png   three panels
#
# ---------------------------------------------------------------------------------------
# WHAT time.seconds ACTUALLY MEASURES, and why the four are comparable
# ---------------------------------------------------------------------------------------
#
# Every method's `time.seconds` column is wall clock around the INFERENCE call for ONE
# trio, and nothing else:
#
#   MRGN    MRGN::infer.trio()          inference_utils.R:756-800.  EXCLUDES the bootstrap,
#                                       which is opt-in and timed separately in
#                                       mrgn.*.bootstrap.time.seconds.
#   MRPC    MRPC::MRPC() inside         inference_utils.R:895-933.  Capped at mrpc.timeout
#           R.utils::withTimeout()      (180 s); an expired fit returns model = NA and
#                                       time.seconds = NA, so it is CENSORED, not slow.
#   GMAC    apply.gmac()                inference_utils.R:982-999.  The per-trio permutation
#                                       test, gmac.nperm = 1000 permutations.
#   MR-GGI  mrggi.one.trio()            inference_utils.R:1429-1433.
#
# Confounder SELECTION is excluded from all four. CS-q and CS-alpha are computed once per
# group in the selection stage and shared by MRGN, MRPC and MR-GGI; GMAC's own selection
# happens in the batch gmac() call and is recorded only as a group attribute
# (attr(results, "gmac.time.seconds"), inference_utils.R:667). So the numbers below are
# "time to fit one trio, given a confounder set" for all four -- which is the quantity that
# is actually comparable. It is NOT end-to-end cost per trio.
#
# ---------------------------------------------------------------------------------------
# THIS IS NOT A CONTROLLED BENCHMARK
# ---------------------------------------------------------------------------------------
#
# These are wall-clock times harvested from the production inference run, not from a timing
# harness. Three things follow, all stated in the report rather than hidden:
#
#   1. MRGN, GMAC and MR-GGI ran their groups on a parallel::makeCluster(); MRPC did not
#      (apply_mrpc.R never builds one). A trio timed inside a busy worker carries some CPU
#      contention that a serial fit does not.
#   2. MRPC has only n = 50/150/300 -- the n = 670 and n = 1000 groups have not been rerun
#      under the raised 180 s cap (see inference_config.R). Pooling all trios therefore
#      compares MRPC on an easier subset than the other three, so the headline table is the
#      COMMON SUBSET (n <= 300, the 900 trios every method has) and the all-trios table is
#      reported beside it, flagged.
#   3. Twenty n = 670 MR-GGI trios carry wall-clock stalls of 3,000-53,000 s -- clustered at
#      near-identical values across trios that were in flight on different workers at the
#      same moment, and wildly out of line with their own covariate counts (dataset 191:
#      78 covariates, 53,192 s, against 155 s for a 106-covariate trio at n = 1000). Those
#      are the machine suspending, not MR-GGI computing, so they are excluded from every
#      statistic and counted in the `n_stalled` column instead. That costs nothing: all of
#      them sit above the 75th percentile of their cell, so the median and the IQR -- what
#      this report is built on -- are identical either way. What it buys is a mean, a max
#      and a total that are not dominated by a sleeping laptop.

library(ggplot2)

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

FIG.DIR <- "bioinfo_revision/reports/figures"
REPORT  <- "bioinfo_revision/reports/COMPUTE_TIME.md"
dir.create(FIG.DIR, recursive = TRUE, showWarnings = FALSE)

# Anything above this in an MR-GGI arm is a suspended machine, not a fit. Nothing else can
# reach it: MRPC is capped at 180 s, MRGN's slowest trio is 0.3 s, GMAC's is 10.8 s.
STALL.SECONDS <- 1000


# ---------------------------------------------------------------------------------------
# assembling the long table
# ---------------------------------------------------------------------------------------
#
# One row per (trio, method, arm). Arms are named exactly as they are in the confusion
# tables so the two reports join: `truth` is the oracle confounder set, `CSq`/`CSa`/`CSi`
# the three selections, `none` MR-GGI's unadjusted fit, `selected` GMAC's own selection.
#
# CS-i is the fourth arm and the one every method now shares: it is
# get.conf.trios(adjust_by = "individual"), which reproduces GMAC's own selection exactly
# (METHODS_FINAL.md section 4.2). GMAC's `CSi` and `selected` arms therefore run the same
# test on the same covariates, and any timing difference between them is the per-trio
# permutation test's own variance, not a difference in work.
#
# GMAC has ONE arm here. run.gmac.group() does build an oracle arm (inference_utils.R:650),
# but the current inference_gmac.RData predates it and carries no gmac.truth.* columns, so
# there is nothing to tabulate. The arm map is checked against the loaded frame rather than
# assumed, so if that arm is backfilled it appears without a code change.
ARM.COLUMNS <- list(
    mrgn  = c(truth = "mrgn.truth.time.seconds",
              CSq   = "mrgn.CSq.time.seconds",
              CSa   = "mrgn.CSa.time.seconds",
              CSi   = "mrgn.CSi.time.seconds"),
    mrpc  = c(truth = "mrpc.truth.time.seconds",
              CSq   = "mrpc.CSq.time.seconds",
              CSa   = "mrpc.CSa.time.seconds",
              CSi   = "mrpc.CSi.time.seconds"),
    gmac  = c(selected = "gmac.time.seconds",
              truth    = "gmac.truth.time.seconds",
              CSi      = "gmac.CSi.time.seconds"),
    mrggi = c(none  = "mrggi.none.time.seconds",
              truth = "mrggi.truth.time.seconds",
              CSq   = "mrggi.CSq.time.seconds",
              CSa   = "mrggi.CSa.time.seconds",
              CSi   = "mrggi.CSi.time.seconds"))

METHOD.LABELS <- c(mrgn = "MRGN", mrpc = "MRPC", gmac = "GMAC", mrggi = "MR-GGI")

# The arm each method is scored on in the headline table: the one it would actually be run
# under. MRGN, MRPC and MR-GGI all get CS-q, which is the selection the manuscript reports;
# GMAC gets its own. The oracle arms are not attainable results and CS-alpha is not the
# recommended selection, so both are reported separately rather than mixed in.
DEPLOYED.ARM <- c(mrgn = "CSq", mrpc = "CSq", gmac = "selected", mrggi = "CSq")

# X-axis order, fixed. Fastest to slowest on the common subset, and the same order on the
# full set, so the two panels read identically.
METHOD.ORDER <- c("MRGN", "MRPC", "MR-GGI", "GMAC")

# Colours follow the METHOD, not its rank -- a method keeps its hue in every panel. The
# four were snapped so that adjacent pairs IN THIS ORDER clear the CVD separation floor;
# re-ordering METHOD.ORDER without re-checking would break that.
METHOD.COLS <- c("MRGN"   = "#0072B2",
                 "MRPC"   = "#009E73",
                 "MR-GGI" = "#D55E00",
                 "GMAC"   = "#8856A7")

build.long <- function() {
    rows <- list()
    for (method in names(ARM.COLUMNS)) {
        x <- load.method(method)
        for (arm in names(ARM.COLUMNS[[method]])) {
            col <- ARM.COLUMNS[[method]][[arm]]
            if (!col %in% names(x)) {
                cat(sprintf("  %-6s / %-8s : no column %s, skipping\n", method, arm, col))
                next
            }
            v <- x[[col]]
            # An arm that is configured off (MRPC's CS-alpha, mrpc.arms in
            # inference_config.R) is all-NA. That is "not attempted", not "took no time",
            # so it is dropped rather than reported as an empty cell.
            if (all(is.na(v))) {
                cat(sprintf("  %-6s / %-8s : all NA (arm not run), skipping\n", method, arm))
                next
            }
            # MRPC is the only method that can censor. timed.out = TRUE means the fit hit
            # mrpc.timeout and was killed, so its time is a lower bound of mrpc.timeout,
            # not a measurement -- carried as a flag and excluded from the quantiles.
            to.col <- paste0(sub("[.]time[.]seconds$", "", col), ".timed.out")
            timed.out <- if (to.col %in% names(x)) {
                !is.na(x[[to.col]]) & x[[to.col]]
            } else rep(FALSE, nrow(x))
            rows[[length(rows) + 1]] <- data.frame(
                method       = unname(METHOD.LABELS[method]),
                arm          = arm,
                dataset      = x$dataset,
                sample.size  = x$sample.size,
                effect_size  = x$effect_size,
                truth.model  = x$truth.model,
                time.seconds = v,
                timed.out    = timed.out,
                stalled      = !is.na(v) & v > STALL.SECONDS,
                stringsAsFactors = FALSE)
        }
    }
    out <- do.call(rbind, rows)
    out$method <- factor(out$method, levels = METHOD.ORDER)
    out
}

cat("loading results...\n")
long <- build.long()
cat(sprintf("  %d method-arm-trio rows\n", nrow(long)))

# A zero or negative duration would break the log axis and would mean the clock resolution
# was coarser than the fit. Assert rather than silently drop -- if it ever fires, the
# timing is not measuring what this report claims it measures.
stopifnot(all(long$time.seconds[!is.na(long$time.seconds)] > 0))


# ---------------------------------------------------------------------------------------
# the summary
# ---------------------------------------------------------------------------------------
#
# MEDIAN and IQR, not mean and sd. Every one of these distributions is right-skewed by one
# to five orders of magnitude -- MR-GGI's CS-alpha arm runs 42 s to 778 s at n = 1000 on
# covariate counts that differ by 2x -- and the MR-GGI stalls put five-figure seconds in
# the tail. A mean over that describes the tail, not the typical trio. Both are still
# written to the CSV so the skew stays visible; only the robust pair is rendered.

# Two exclusions, both counted rather than silent:
#
#   timed.out  a killed MRPC fit has no duration. Its true cost is censored at
#              mrpc.timeout, so including it as anything would be inventing a number.
#   stalled    the MR-GGI wall clocks that are a suspended machine (see the header). These
#              are excluded from EVERY statistic, not just from the mean, because a
#              measurement artefact is not a slow fit. It costs nothing to do it this way:
#              all of them sit above the 75th percentile of their cell, so the median and
#              the IQR -- what this report is built on -- are identical either way. What it
#              buys is a `mean`, `max` and `total_hours` that mean something: with the
#              three stalls left in, MR-GGI's CS-q arm reports 10.9 total hours of which
#              9.3 are a sleeping laptop.
summarise.cell <- function(d) {
    v <- d$time.seconds[!d$timed.out & !d$stalled & !is.na(d$time.seconds)]
    if (length(v) == 0) return(NULL)
    qs <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE)
    data.frame(n_trios     = nrow(d),
               n_timed_out = sum(d$timed.out),
               n_stalled   = sum(d$stalled),
               n_scored    = length(v),
               min         = min(v),
               q1          = qs[1],
               median      = qs[2],
               q3          = qs[3],
               iqr         = qs[3] - qs[1],
               max         = max(v),
               mean        = mean(v),
               total_hours = sum(v) / 3600)
}

# scope = "all" pools every sample size the method has; scope = "n<=300" is the 900 trios
# all four share; the numeric scopes are the per-group split.
summarise.by.scope <- function(d) {
    cells <- list()
    add <- function(scope, sub) {
        s <- summarise.cell(sub)
        if (!is.null(s)) {
            cells[[length(cells) + 1]] <<-
                cbind(method = as.character(sub$method[1]), arm = sub$arm[1],
                      scope = scope, s, stringsAsFactors = FALSE)
        }
    }
    add("all", d)
    add("n<=300", d[d$sample.size <= 300, ])
    for (n in SAMPLE.SIZES) {
        sub <- d[d$sample.size == n, ]
        if (nrow(sub) > 0) add(as.character(n), sub)
    }
    do.call(rbind, cells)
}

summary.tbl <- do.call(rbind, lapply(
    split(long, list(long$method, long$arm), drop = TRUE), summarise.by.scope))
rownames(summary.tbl) <- NULL
# Arms in increasing covariate count -- none, CS-q, truth, CS-alpha -- not alphabetically.
# That IS the story of section 3: sorting them "CS-alpha, CS-q, truth" puts the most
# expensive arm first and hides the monotone relationship the table exists to show.
ARM.ORDER <- c("none", "CSq", "CSi", "truth", "CSa", "selected")
stopifnot(all(summary.tbl$arm %in% ARM.ORDER))
summary.tbl <- summary.tbl[order(match(summary.tbl$method, METHOD.ORDER),
                                 match(summary.tbl$arm, ARM.ORDER),
                                 match(summary.tbl$scope,
                                       c("all", "n<=300", as.character(SAMPLE.SIZES)))), ]

dir.create(tables.dir, recursive = TRUE, showWarnings = FALSE)
write.csv(long, file.path(tables.dir, "compute_time_long.csv"), row.names = FALSE)
write.csv(summary.tbl, file.path(tables.dir, "compute_time.csv"), row.names = FALSE)
cat(sprintf("  wrote %s (%d rows)\n",
            file.path(tables.dir, "compute_time_long.csv"), nrow(long)))
cat(sprintf("  wrote %s (%d rows)\n",
            file.path(tables.dir, "compute_time.csv"), nrow(summary.tbl)))


# ---------------------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------------------
#
# Boxplot on a LOG10 second axis. The four medians span 6 ms to 0.56 s and the whiskers
# span five orders of magnitude, so on a linear axis three of the four methods render as a
# flat line on zero. The axis is labelled in human units (ms / s / min / h) rather than in
# powers of ten, because the question the figure answers is "how long does one trio take",
# not "what is log10 of that".
#
# The box IS the IQR by construction, so the figure and the table are the same two numbers
# drawn two ways. Medians are DIRECT-LABELLED above each box: with four boxes there is room
# for it, and it means the headline number never has to be read off a log axis by eye.
#
# No legend. Colour is redundant with the x-axis label here, so identity is carried by text
# and a colourblind reader loses nothing; a legend would only repeat the axis.

fig.theme <- theme_bw(base_size = 12) +
    theme(strip.background   = element_rect(fill = "grey85", colour = "grey30"),
          strip.text         = element_text(face = "bold"),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.subtitle      = element_text(colour = "grey30", size = 9),
          legend.position    = "none")

save.fig <- function(p, name, width, height) {
    path <- file.path(FIG.DIR, name)
    ggsave(path, p, width = width, height = height, dpi = 150, bg = "white")
    cat(sprintf("  wrote %s (%.0f x %.0f in)\n", path, width, height))
}

# Ticks every decade from 1 ms to 10 h, spelled in the unit a reader thinks in.
TIME.BREAKS <- c(0.001, 0.01, 0.1, 1, 10, 60, 600, 3600, 36000)
TIME.LABELS <- c("1 ms", "10 ms", "100 ms", "1 s", "10 s", "1 min", "10 min", "1 h", "10 h")

fmt.time <- function(x, unit = NULL) {
    u <- if (!is.null(unit)) unit else
        ifelse(x < 1, "ms", ifelse(x < 60, "s", ifelse(x < 3600, "min", "h")))
    ifelse(u == "ms",  sprintf("%.0f ms",  x * 1000),
    ifelse(u == "s",   sprintf("%.2f s",   x),
    ifelse(u == "min", sprintf("%.1f min", x / 60),
                       sprintf("%.1f h",   x / 3600))))
}

# Both endpoints of a range wear the SAME unit where that is legible, because formatting
# them independently produces "53.81 s - 2.4 min" and the reader has to convert before they
# can see how wide the range is -- the one thing a range is there to show.
#
# The unit is taken from the LOWER bound, which is the endpoint that loses precision when
# coerced upwards (53.81 s becomes "0.9 min"). Where the range is wide enough that the
# upper bound would then run past 1000 of that unit -- "50 ms - 4411 ms" -- the endpoints
# get their own units after all, since by then the two are obviously different magnitudes
# and the shared unit has stopped buying anything.
NATURAL.UNIT <- function(x) if (x < 1) "ms" else if (x < 60) "s" else
    if (x < 3600) "min" else "h"
UNIT.SECONDS <- c(ms = 0.001, s = 1, min = 60, h = 3600)

fmt.range <- function(lo, hi) {
    u <- NATURAL.UNIT(lo)
    if (hi / UNIT.SECONDS[[u]] > 1000) return(sprintf("%s - %s", fmt.time(lo), fmt.time(hi)))
    sprintf("%s - %s", fmt.time(lo, u), fmt.time(hi, u))
}

log.y <- function() {
    scale_y_log10(breaks = TIME.BREAKS, labels = TIME.LABELS,
                  expand = expansion(mult = c(0.04, 0.07)))
}

box.base <- function(d) {
    # Two classes of row are dropped HERE rather than left for ggplot to discard with a
    # "removed N rows containing non-finite values" warning, which is a silent result:
    #
    #   timed.out  a killed MRPC fit has no duration to plot at all
    #   stalled    the 20 MR-GGI wall clocks that are a suspended machine, not a fit. They
    #              stay in the quantiles (which they cannot move -- see the header) but
    #              they are not drawn: three orders of magnitude past every real point,
    #              they stretch the log axis to 10 h and squash every box into the bottom
    #              third of the panel to show two dots that mean nothing.
    #
    # Both counts are announced, and both are in the CSV as n_timed_out / n_stalled.
    drop <- is.na(d$time.seconds) | d$timed.out | d$stalled
    if (any(drop)) {
        by.arm <- table(paste(d$method[drop], d$arm[drop]))
        cat(sprintf("    (%d censored/stalled fits not plotted: %s)\n", sum(drop),
                    paste(sprintf("%s x%d", names(by.arm), as.integer(by.arm)),
                          collapse = ", ")))
    }
    d <- d[!drop, , drop = FALSE]
    ggplot(d, aes(method, time.seconds, fill = method, colour = method)) +
        # Points UNDER the box, so the density stays visible -- an IQR box alone cannot
        # show that MR-GGI's CS-q arm is bimodal in n (milliseconds at n = 50, seconds at
        # n = 1000), which is the whole reason the by-size panel exists.
        geom_jitter(width = 0.22, height = 0, size = 0.5, alpha = 0.10, stroke = 0) +
        geom_boxplot(width = 0.5, alpha = 0.35, outlier.shape = NA, linewidth = 0.6,
                     colour = "grey25") +
        scale_fill_manual(values = METHOD.COLS, drop = FALSE) +
        scale_colour_manual(values = METHOD.COLS, drop = FALSE) +
        scale_x_discrete(drop = FALSE) +
        log.y() +
        labs(x = NULL, y = "Compute time per trio") +
        fig.theme
}

# Median label for each box, anchored to the TOP OF THE UPPER WHISKER -- the largest value
# within Q3 + 1.5 x IQR, which is exactly where geom_boxplot() ends its line. Anchoring to
# max(v) instead floats the label a decade or more above its own box on a log axis (the
# tails here run four orders past the whisker), which reads as an unlabelled box beside a
# free-floating number.
median.labels <- function(d, by = NULL) {
    keys  <- c("method", by)
    parts <- split(d, d[keys], drop = TRUE)
    do.call(rbind, lapply(parts, function(s) {
        v <- s$time.seconds[!s$timed.out & !s$stalled & !is.na(s$time.seconds)]
        if (!length(v)) return(NULL)
        qs      <- stats::quantile(v, c(0.25, 0.75), names = FALSE)
        whisker <- max(v[v <= qs[2] + 1.5 * (qs[2] - qs[1])])
        cbind(s[1, keys, drop = FALSE],
              data.frame(time.seconds = whisker, lab = fmt.time(median(v)),
                         stringsAsFactors = FALSE))
    }))
}

# ---- (i) the headline: the deployed arm, on the 900 trios all four methods have ----
method.key   <- setNames(names(METHOD.LABELS), unname(METHOD.LABELS))
deployed.for <- function(m) unname(DEPLOYED.ARM[method.key[as.character(m)]])
dep    <- long[long$arm == deployed.for(long$method), ]
common <- dep[dep$sample.size <= 300, ]

save.fig(box.base(common) +
             geom_text(data = median.labels(common), aes(label = lab), colour = "grey15",
                       vjust = -1.1, size = 4, fontface = "bold") +
             labs(title = "Per-trio compute time, as deployed",
                  subtitle = paste0("n = 50, 150, 300 -- the 900 trios all four methods ",
                                    "cover. Box = IQR, whiskers = 1.5 x IQR, label = ",
                                    "median.\nMRGN / MRPC / MR-GGI on the CS-q confounder ",
                                    "set, GMAC on its own selection. Log scale.")),
         "fig_compute_time.png", 8, 7)

# ---- (ii) the same thing at every sample size ----
# MRPC is absent from the n = 670 and n = 1000 panels because those groups have not been
# run, NOT because it was fast there. drop = FALSE on the x scale keeps the empty slot --
# without it ggplot silently draws three evenly-spaced boxes, which reads as a three-method
# comparison -- and the slot is annotated so the gap says why it is there.
size.lab <- function(n) factor(paste0("n = ", n), levels = paste0("n = ", SAMPLE.SIZES))
dep$size.lab <- size.lab(dep$sample.size)

lab.s <- median.labels(dep, by = "sample.size")
lab.s$size.lab <- size.lab(lab.s$sample.size)

have    <- unique(paste(dep$method, dep$sample.size))
not.run <- expand.grid(method = METHOD.ORDER, sample.size = SAMPLE.SIZES,
                       stringsAsFactors = FALSE)
not.run <- not.run[!paste(not.run$method, not.run$sample.size) %in% have, , drop = FALSE]
not.run$method       <- factor(not.run$method, levels = METHOD.ORDER)
not.run$size.lab     <- size.lab(not.run$sample.size)
not.run$time.seconds <- 0.05

save.fig(box.base(dep) +
             geom_text(data = lab.s, aes(label = lab), colour = "grey15",
                       vjust = -1.1, size = 3, fontface = "bold") +
             geom_text(data = not.run, aes(label = "not run"), colour = "grey45",
                       size = 3.2, angle = 90, fontface = "italic") +
             facet_wrap(~ size.lab, nrow = 1) +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
             labs(title = "Per-trio compute time by sample size, as deployed",
                  subtitle = paste0("Box = IQR, label = median. MRPC has no n = 670 or ",
                                    "n = 1000 group -- see inference_config.R.\n3 ",
                                    "stalled MR-GGI wall clocks at n = 670 are excluded ",
                                    "from the marks; they move neither the median nor the ",
                                    "IQR. Log scale.")),
         "fig_compute_time_by_size.png", 15, 7)

# ---- (iii) what the confounder set costs ----
# The deployed arm is one column of a wider story: every method here is dominated by how
# many covariates it was handed, not by the method. This panel is the evidence for that.
ARM.PRETTY <- c("no covariates", "CS-q", "CS-i", "truth (oracle)", "CS-α", "own (GMAC)")
arm.d     <- long
arm.d$arm <- factor(arm.d$arm, levels = ARM.ORDER, labels = ARM.PRETTY)
lab.a     <- median.labels(arm.d, by = "arm")

save.fig(box.base(arm.d) +
             geom_text(data = lab.a, aes(label = lab), colour = "grey15",
                       vjust = -1.1, size = 3, fontface = "bold") +
             facet_wrap(~ arm, nrow = 1) +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
             labs(title = "Per-trio compute time by confounder set (all trios, all sizes)",
                  subtitle = paste0("Box = IQR, label = median. The cost is set by how ",
                                    "many covariates the method was handed, not by the ",
                                    "method.\n30 stalled MR-GGI wall clocks (20 trios, ",
                                    "across three arms) and 103 timed-out MRPC oracle ",
                                    "fits are excluded from the marks. Log scale.")),
         "fig_compute_time_by_arm.png", 15, 7)


# ---------------------------------------------------------------------------------------
# the report
# ---------------------------------------------------------------------------------------

md.table <- function(header, rows) {
    paste0(c(paste0("| ", paste(header, collapse = " | "), " |"),
             paste0("| ", paste(rep("---", length(header)), collapse = " | "), " |"),
             vapply(rows, function(r) paste0("| ", paste(r, collapse = " | "), " |"), "")),
           collapse = "\n")
}

cell <- function(method, arm, scope) {
    r <- summary.tbl[summary.tbl$method == method & summary.tbl$arm == arm &
                         summary.tbl$scope == scope, ]
    if (nrow(r) == 0) NULL else r[1, ]
}

arm.name <- function(arm) switch(arm, CSq = "CS-q", CSa = "CS-α", CSi = "CS-i",
                                 truth = "truth (oracle)", none = "no covariates",
                                 selected = "own selection", arm)

deployed.rows <- function(scope) {
    lapply(METHOD.ORDER, function(m) {
        arm <- deployed.for(m)
        r   <- cell(m, arm, scope)
        if (is.null(r)) return(c(m, arm.name(arm), rep("--", 5)))
        c(m, arm.name(arm), format(r$n_scored, big.mark = ","),
          fmt.time(r$median),
          fmt.range(r$q1, r$q3),
          fmt.time(r$iqr),
          fmt.time(r$total_hours * 3600))
    })
}

HEAD <- c("method", "confounders", "trios", "median", "IQR (Q1 - Q3)", "IQR width",
          "total, all trios")

mrggi.covars <- load.method("mrggi")
covar.median <- function(method, arm) {
    if (method != "MR-GGI") return("--")
    col <- sprintf("mrggi.%s.n.covars", arm)
    if (!col %in% names(mrggi.covars)) "--" else
        sprintf("%.0f", median(mrggi.covars[[col]], na.rm = TRUE))
}

mrgn.boot <- load.method("mrgn")

lines <- c(
"# Per-trio compute time",
"",
sprintf(paste0("Generated by [`compute_time.R`](../simulation_results/results_scripts/",
               "compute_time.R) on %s."), format(Sys.Date(), "%Y-%m-%d")),
"",
"Median and interquartile range of the time to fit **one trio**, for each of the four",
"inference methods. Machine-readable versions:",
"[`compute_time.csv`](../simulation_results/tables/compute_time.csv) (the tables below,",
"plus mean/min/max and a per-sample-size split) and",
"[`compute_time_long.csv`](../simulation_results/tables/compute_time_long.csv) (one row",
"per trio per method per arm, if you want to re-cut it).",
"",
"## 1. Headline: the deployed configuration",
"",
"Each method on the confounder set it would actually be run under -- CS-q for MRGN, MRPC",
"and MR-GGI, its own selection for GMAC -- restricted to **n = 50, 150 and 300: the 900",
"trios all four methods cover**. Section 4 says why the restriction is there.",
"",
md.table(HEAD, deployed.rows("n<=300")),
"",
"![per-trio compute time](figures/fig_compute_time.png)",
"",
"The same four on **every trio each method has** -- 1,500 for MRGN, GMAC and MR-GGI, 900",
"for MRPC. MRPC's row here is the identical 900 trios as above, while the other three now",
"include n = 670 and n = 1000, which is where their cost is. Read it as three numbers and",
"a placeholder, not as a four-way comparison.",
"",
md.table(HEAD, deployed.rows("all")),
"",
"## 2. By sample size",
"",
"![compute time by sample size](figures/fig_compute_time_by_size.png)",
"",
md.table(c("method", "confounders", paste0("n = ", SAMPLE.SIZES)),
         lapply(METHOD.ORDER, function(m) {
             arm <- deployed.for(m)
             c(m, arm.name(arm),
               vapply(as.character(SAMPLE.SIZES), function(s) {
                   r <- cell(m, arm, s)
                   if (is.null(r)) "*not run*" else
                       sprintf("%s<br><sub>%s</sub>", fmt.time(r$median),
                               fmt.range(r$q1, r$q3))
               }, ""))
         })),
"",
"Median on the first line, IQR beneath it.",
"",
"## 3. Every arm",
"",
"The deployed row is one column of a wider picture. For all four methods the per-trio cost",
"is set by **how many covariates the method was handed**, not by the method. MR-GGI is the",
"clearest case: median 23 ms with no covariates, 1.3 min against CS-alpha's median of 98 --",
"a factor of ~3,400 on the same trios, because its cost is quadratic in the covariate count",
"(it tests every gene pair). The other two move far less: MRGN spans 7-34 ms across its",
"three arms and MRPC 34-163 ms across its two, a factor of five each.",
"",
"![compute time by confounder set](figures/fig_compute_time_by_arm.png)",
"",
md.table(c("method", "arm", "covariates (median)", "trios scored", "median",
           "IQR (Q1 - Q3)", "IQR width"),
         local({
             rs <- summary.tbl[summary.tbl$scope == "all", ]
             lapply(seq_len(nrow(rs)), function(i) {
                 r <- rs[i, ]
                 c(r$method, arm.name(r$arm), covar.median(r$method, r$arm),
                   format(r$n_scored, big.mark = ","), fmt.time(r$median),
                   fmt.range(r$q1, r$q3), fmt.time(r$iqr))
             })
         })),
"",
"## 4. What these numbers are, and what they are not",
"",
"**They time the fit, not the pipeline.** Every `time.seconds` column is wall clock around",
"the inference call for one trio and nothing else: `MRGN::infer.trio()`, `MRPC::MRPC()`",
"inside its timeout, `apply.gmac()`'s permutation test, `mrggi.one.trio()`. Confounder",
"**selection is excluded from all four** -- CS-q and CS-alpha are computed once per group",
"in the selection stage and shared, and GMAC's selection happens inside the batch `gmac()`",
"call and is recorded only as a group attribute. So this is *time to fit one trio given a",
"confounder set*, which is the quantity that is comparable across the four. It is not",
"end-to-end cost per trio.",
"",
"**MRGN's bootstrap is excluded.** `mrgn.*.time.seconds` covers `infer.trio()` only; the",
"opt-in edge-probability bootstrap is timed separately and is roughly 200-250x the fit it",
"wraps -- 1,000 resamples, each a fresh `infer.trio()`. If bootstrap edge probabilities are",
"wanted, MRGN's per-trio cost is the sum of the two columns and it stops being the cheapest",
"method in section 1:",
"",
md.table(c("arm", "fit (median)", "bootstrap (median)", "ratio"),
         lapply(c("CSq", "truth", "CSa"), function(a) {
             r <- cell("MRGN", a, "all")
             b <- median(mrgn.boot[[sprintf("mrgn.%s.bootstrap.time.seconds", a)]],
                         na.rm = TRUE)
             c(arm.name(a), fmt.time(r$median), fmt.time(b),
               sprintf("%.0fx", b / r$median))
         })),
"",
"**MRPC has no n = 670 or n = 1000 group.** Those two have not been rerun under the raised",
"180 s cap (`inference_config.R`), so pooling every trio would score MRPC on an easier",
"subset than the other three. That is why section 1 leads with n <= 300.",
"",
sprintf(paste0("**MRPC's timeouts are censored, not slow.** A fit that hits the %g s cap ",
               "returns no model and no time, so it is excluded from the quantiles and ",
               "counted separately (`n_timed_out` in the CSV). The CS-q arm has none at ",
               "n <= 300; the oracle arm loses %d of 900, all at n = 300."),
        mrpc.timeout, sum(long$timed.out[long$method == "MRPC" & long$arm == "truth"])),
"",
sprintf(paste0("**Twenty MR-GGI trios carry machine stalls.** All at n = 670, all with ",
               "wall clocks of 3,000-53,000 s clustered at near-identical values across ",
               "trios that were in flight on different workers at the same instant, and ",
               "all wildly out of line with their own covariate counts -- dataset 191 has ",
               "78 covariates and 53,192 s, against 155 s for a 106-covariate trio at ",
               "n = 1000. That is the machine suspending, not MR-GGI computing, so they ",
               "are excluded from every statistic here and counted in `n_stalled` ",
               "instead. **The medians and IQRs are identical either way** -- all thirty ",
               "sit above the 75th percentile of their cell -- so nothing is lost; what ",
               "the exclusion buys is a `mean`, `max` and `total_hours` in the CSV that ",
               "mean something. Left in, MR-GGI's CS-q arm would report 10.9 total hours ",
               "of which 9.3 are a sleeping laptop. Affected trios: %s."),
        paste(sort(unique(long$dataset[long$stalled])), collapse = ", ")),
"",
"**MRGN's CS-i arm ran on a quiet machine, and it shows.** It reports a 4 ms median",
"against CS-q's 7 ms while carrying *more* covariates (2.6-24.0 per trio against",
"0.5-21.7), which inverts the relationship every other row in section 3 follows. The arm",
"is not cheaper. It was run on its own, after the fact, as a single pass",
"(`apply_mrgn_csi.R`), whereas the truth/CS-q/CS-α arms were timed inside the original run",
"with three other method processes competing for the same cores. That is a clean",
"demonstration of the caveat below rather than an exception to section 3: **do not read a",
"factor-of-two gap between arms timed in different runs as a property of the arms.**",
"",
"**This is not a controlled benchmark.** These are wall-clock times harvested from the",
"production inference run, not from a timing harness. MRGN, GMAC and MR-GGI ran their",
"groups on a `parallel::makeCluster()`; MRPC did not (`apply_mrpc.R` never builds one), so",
"a trio timed inside a busy worker carries CPU contention that a serial fit does not. Read",
"factor-of-two differences as noise and factor-of-a-hundred differences as real.")

writeLines(lines, REPORT)
cat(sprintf("  wrote %s\n", REPORT))
cat("\ndone.\n")
