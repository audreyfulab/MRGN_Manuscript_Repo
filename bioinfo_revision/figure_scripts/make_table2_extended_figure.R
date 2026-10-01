# Table 2 extended as a two-panel main figure: MRGN with CS-q across the four covariate
# structures.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/figure_scripts/make_table2_extended_figure.R
#
# The "Table2 Extended" sheet of Final_tables/Final_rec_prec_tables.xlsx holds four blocks
# -- one per covariate structure at n = 670 -- each with recall and precision for M0-M4 and
# the T1-T2 edge, under three MRGN arms. That is 144 cells, and the result buried in them
# is hard to see: the common child is what damages MRGN, the intermediate is harmless.
#
# This figure shows the CS-q arm only, as two panels:
#
#   A  recall    x = the four structures, y = recall,    six lines
#   B  precision x = the four structures, y = precision, six lines
#
# Six lines = the five generating models plus the T1-T2 edge. The edge is not a sixth model
# -- it scores the edge-present class rather than a model label -- so it is drawn black,
# dashed and thicker, to keep it from reading as one.
#
# X ORDER IS NOT THE SHEET ORDER. The sheet runs U, U+Z, U+W, U+W+Z; the figure runs
# U, U+W, U+Z, U+W+Z so the two hazards are introduced one at a time in increasing severity
# and the drop at U+Z reads left to right. Blocks are matched by title, not by position, so
# reordering the sheet cannot silently mis-label the figure.
#
# READS THE WORKBOOK, not confusion_structures.R. The workbook is what ships with the
# manuscript, so sourcing the figure from it means the figure and the table cannot
# disagree. The cost is a string parse; the assertions below cover it.
#
# Writes fig_table2_extended_csq.{pdf,png} to bioinfo_revision/reports/figures/.

library('ggpubr')
library('ggthemes')
library('readxl')

FIG.DIR <- "bioinfo_revision/reports/figures"
dir.create(FIG.DIR, recursive = TRUE, showWarnings = FALSE)

XLSX <- "bioinfo_revision/Final_tables/Final_rec_prec_tables.xlsx"
SHEET <- "Table2 Extended"

# columns of the sheet: A row label, B/C MRGN (Truth), D/E MRGN (CS-q), F/G MRGN (CS-alpha)
COL.LABEL <- 1
COL.RECALL <- 4
COL.PRECISION <- 5

ROW.LABELS <- c("M0", "M1", "M2", "M3", "M4", "T1-T2 Edge")

# substring -> short scenario name, in the order they are plotted. Matched on a substring
# because the fourth block title carries an em-dash suffix naming it the main simulation.
SCENARIOS <- list(
    list(key = "only",                        label = "U only"),
    list(key = "+ intermediate (U + W)",      label = "U + W"),
    list(key = "+ common child (U + Z)",      label = "U + Z"),
    list(key = "intermediate + common child", label = "U + W + Z"))

SERIES.LABELS <- c(M0 = "Null (M0)",
                   M1 = "Mediation (M1)",
                   M2 = "V-structure (M2)",
                   M3 = "Cond. indep. (M3)",
                   M4 = "Fully connected (M4)",
                   `T1-T2 Edge` = "T1 - T2 edge")

# Okabe-Ito for the five models; the edge is black, and set apart again by linetype and
# width below, so it does not have to be told apart by hue alone.
SERIES.COLORS <- c("Null (M0)"            = "#0072B2",
                   "Mediation (M1)"       = "#D55E00",
                   "V-structure (M2)"     = "#009E73",
                   "Cond. indep. (M3)"    = "#CC79A7",
                   "Fully connected (M4)" = "#E69F00",
                   "T1 - T2 edge"         = "#000000")

# ---------------------------------------------------------------------------------------
# read and parse
# ---------------------------------------------------------------------------------------

raw <- as.data.frame(readxl::read_excel(XLSX, sheet = SHEET, col_names = FALSE,
                                        .name_repair = "minimal"))

cat("=== Table 2 extended figure (MRGN with CS-q) ===\n")
cat("  read ", XLSX, " sheet '", SHEET, "' | ", nrow(raw), " x ", ncol(raw), "\n", sep = "")

# "0.683 (0.060)" -> c(estimate = 0.683, se = 0.060). A cell that does not match this shape
# is an error rather than an NA: it means the sheet was reformatted and every number below
# is suspect.
parse.cell <- function(x, where) {
    x <- trimws(as.character(x))
    m <- regmatches(x, regexec("^([0-9.]+)\\s*\\(([0-9.]+)\\)$", x))[[1]]
    if (length(m) != 3)
        stop("cannot parse '", x, "' at ", where, " -- expected 'estimate (SE)'")
    c(estimate = as.numeric(m[2]), se = as.numeric(m[3]))
}

# Blocks are found by their M0 row, and the title is three rows above it. Position in the
# sheet is never assumed.
m0.rows <- which(trimws(as.character(raw[[COL.LABEL]])) == "M0")
if (length(m0.rows) != length(SCENARIOS))
    stop("found ", length(m0.rows), " blocks in '", SHEET, "', expected ",
         length(SCENARIOS))

blocks <- lapply(m0.rows, function(r0) {
    rows <- r0:(r0 + length(ROW.LABELS) - 1L)
    labs <- trimws(as.character(raw[[COL.LABEL]][rows]))
    if (!identical(labs, ROW.LABELS))
        stop("block at row ", r0, " has labels ", paste(labs, collapse = "/"),
             " -- expected ", paste(ROW.LABELS, collapse = "/"))
    list(title = trimws(as.character(raw[[2]][r0 - 3L])), rows = rows)
})

titles <- vapply(blocks, function(b) b$title, character(1))
cat("  blocks found (sheet order):\n")
for (t in titles) cat("    ", t, "\n", sep = "")

# ---------------------------------------------------------------------------------------
# long frame, in plotting order
# ---------------------------------------------------------------------------------------

long <- list()
for (sc in SCENARIOS) {
    hit <- which(grepl(sc$key, titles, fixed = TRUE))
    if (length(hit) != 1)
        stop("scenario key '", sc$key, "' matched ", length(hit), " block titles")
    blk <- blocks[[hit]]

    for (k in seq_along(ROW.LABELS)) {
        r <- blk$rows[k]
        rec <- parse.cell(raw[[COL.RECALL]][r],
                          sprintf("row %d, recall (%s)", r, ROW.LABELS[k]))
        prc <- parse.cell(raw[[COL.PRECISION]][r],
                          sprintf("row %d, precision (%s)", r, ROW.LABELS[k]))
        long[[length(long) + 1L]] <- data.frame(
            scenario = sc$label,
            series   = unname(SERIES.LABELS[ROW.LABELS[k]]),
            metric   = c("Recall", "Precision"),
            estimate = c(rec[["estimate"]], prc[["estimate"]]),
            se       = c(rec[["se"]], prc[["se"]]),
            stringsAsFactors = FALSE)
    }
}
dat <- do.call(rbind, long)

dat$scenario <- factor(dat$scenario, levels = vapply(SCENARIOS, function(s) s$label,
                                                     character(1)))
dat$series <- factor(dat$series, levels = unname(SERIES.LABELS))

stopifnot(nrow(dat) == length(SCENARIOS) * length(ROW.LABELS) * 2L)
stopifnot(!any(is.na(dat$estimate)), !any(is.na(dat$se)))
stopifnot(all(dat$estimate >= 0 & dat$estimate <= 1))

# The two numbers the figure exists to show, echoed so a mis-parse is visible in the log.
show <- function(sc, sr, mt) {
    v <- dat$estimate[dat$scenario == sc & dat$series == sr & dat$metric == mt]
    sprintf("%.3f", v)
}
cat("  CS-q recall, Null (M0)         : ",
    paste(vapply(levels(dat$scenario), show, character(1),
                 sr = "Null (M0)", mt = "Recall"), collapse = " -> "), "\n", sep = "")
cat("  CS-q recall, Cond. indep. (M3) : ",
    paste(vapply(levels(dat$scenario), show, character(1),
                 sr = "Cond. indep. (M3)", mt = "Recall"), collapse = " -> "), "\n", sep = "")
cat("  CS-q precision, T1 - T2 edge   : ",
    paste(vapply(levels(dat$scenario), show, character(1),
                 sr = "T1 - T2 edge", mt = "Precision"), collapse = " -> "), "\n", sep = "")

# ---------------------------------------------------------------------------------------
# the panels
# ---------------------------------------------------------------------------------------

# get_legend() and ggarrange() below both call ggplot_build(), which measures text and so
# needs an open device. Under Rscript there is none, so R opens the default pdf() and
# leaves a stray Rplots.pdf in the working directory. Send that to the null device for the
# duration of the build; closed again just before the real output devices are opened.
grDevices::pdf(NULL)

# six series at four x positions, each with an error bar; without a dodge the bars overlap
# wherever two series sit close, which is most of the precision panel
pd <- position_dodge(width = 0.4)

panel <- function(metric, ylab) {
    d <- dat[dat$metric == metric, , drop = FALSE]
    ggplot(d, aes(x = scenario, y = estimate, colour = series, group = series,
                  linetype = series, linewidth = series)) +
        geom_errorbar(aes(ymin = pmax(0, estimate - se), ymax = pmin(1, estimate + se)),
                      width = 0.18, linewidth = 0.5, position = pd) +
        geom_line(position = pd) +
        geom_point(size = 2.6, position = pd) +
        scale_colour_manual(values = SERIES.COLORS) +
        scale_linetype_manual(values = c(rep("solid", 5), "22")) +
        scale_linewidth_manual(values = c(rep(0.9, 5), 1.6)) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        theme_hc() +
        theme(legend.position = 'top', legend.title = element_blank(),
              legend.text = element_text(size = 16),
              axis.text = element_text(size = 16),
              axis.text.x = element_text(size = 15),
              axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 18),
              axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 18),
              # theme_hc draws no vertical gridlines, and a four-category line plot needs
              # them to read a series across the x axis
              panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.4)) +
        xlab("Covariate structure") +
        ylab(ylab)
}

A <- panel("Recall", "Recall")
B <- panel("Precision", "Precision")

# ggpubr::get_legend broke on ggplot2 3.5; this digs the grob out directly. Verified
# against ggplot2 4.0.3, where the guide grob is still named "guide-box".
get_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend_obj <- get_legend(A)

E <- ggarrange(A + theme(legend.position = 'none'),
               B + theme(legend.position = 'none'),
               labels = c("A", "B"),
               nrow = 1, ncol = 2,
               common.legend = F,
               legend = 'top',
               legend.grob = legend_obj,
               font.label = list(size = 16, face = "bold", color = "black")) +
     theme(plot.margin = unit(c(0, 0, 0.4, 0), "cm"))

E <- annotate_figure(E,
    bottom = grid::textGrob(
        paste("MRGN with CS-q, n = 670, 300 trios per structure.",
              "Error bars are +/- 1 Wald SE."),
        gp = grid::gpar(fontsize = 13, col = "grey25")))

# ---------------------------------------------------------------------------------------
# write
# ---------------------------------------------------------------------------------------

invisible(grDevices::dev.off())   # the null device opened above

pdf.path <- file.path(FIG.DIR, "fig_table2_extended_csq.pdf")
pdf(pdf.path, height = 7, width = 14)
plot(E)
invisible(dev.off())
cat(sprintf("  wrote %s (14 x 7 in)\n", pdf.path))

png.path <- file.path(FIG.DIR, "fig_table2_extended_csq.png")
png(png.path, height = 7, width = 14, units = 'in', res = 300)
plot(E)
invisible(dev.off())
cat(sprintf("  wrote %s (14 x 7 in, 300 dpi)\n", png.path))

cat("\ndone.\n")
