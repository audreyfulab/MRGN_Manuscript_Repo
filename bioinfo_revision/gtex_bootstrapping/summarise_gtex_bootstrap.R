# A summary of the GTEx edge-probability bootstrap, as one workbook.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/gtex_bootstrapping/summarise_gtex_bootstrap.R
#
# Writes bioinfo_revision/Final_tables/GTEx_bootstrap_summary.xlsx. Reads
# reports/gtex_bootstrapping/gtex_mrgn_bootstrap.csv, the per-trio output of
# apply_mrgn_gtex.R, so run that first.
#
# WHAT IS BEING SUMMARISED. MRGN on the 3,248 GTEx Whole Blood trios, each refitted on 1,000
# resamples. THERE IS NO GROUND TRUTH HERE -- this is real data, so nothing can be scored
# correct and there is no confusion matrix to build. The bootstrap is the only available
# statement about confidence, and what it supports is two things: whether the resampled
# majority names the same model as the point fit (`boot.agrees`), and how much of the
# resample mass sits behind each individual edge (`boot.p.*`). Both are summarised here.
#
# THE CENTRAL TABLE is point_vs_bootstrap: the point model call cross-tabulated against the
# model a majority vote of the resamples supports. It is NOT a confusion matrix -- neither
# axis is truth -- so it carries a stability rate on the diagonal rather than recall and
# precision. Read it as: of the trios MRGN called M4, how many does resampling still call M4.
#
# ---------------------------------------------------------------------------------------
# TWO PROPERTIES OF THE UNDERLYING COLUMNS THAT THE TABLES HAVE TO ACCOUNT FOR
# ---------------------------------------------------------------------------------------
#
# 1. boot.p.T1T2 AND boot.p.T2T1 ARE THE SAME NUMBER. infer.trio() returns six indicators,
#    b11, b12, b21, b22, V1:T1, V1:T2, and b12 and b22 are one indicator for the T1-T2
#    ADJACENCY rather than two orientation-specific ones -- they are identical in every one
#    of the 3,248 trios here, and in every resample behind them. Orientation of T1-T2 is
#    decided by the two interaction terms, not by these. So the edge_support sheet reports
#    the T1-T2 edge ONCE, labelled unoriented, instead of listing one number twice as if it
#    were evidence about direction. The script re-checks the identity and stops if it ever
#    fails, since a future MRGN could split them and the collapse would then be wrong.
#
# 2. boot.min.edge.prob IS NOT ALWAYS AN EDGE PROBABILITY. It is built in
#    boostrap_edge_probabilities() as min(indicator.means[indicator.means >= 0.5]) over ALL
#    SIX indicators, and the last two are the V1:T1 and V1:T2 INTERACTION terms, which are
#    not edges. In 891 of the 3,248 trios the minimum is one of those two rather than any
#    edge -- 870 where it sits below the weakest supported edge, and 21 where no edge is
#    supported at all yet the column still carries a value. So the column is a floor on the
#    whole supported indicator vector, not on the skeleton.
#
#    Both readings are reported. `min.supported.indicator` is the column as recorded, which
#    is what apply_mrgn_gtex.R wrote and what any earlier reading of this file used;
#    `min.supported.edge` is recomputed here over the four edge indicators only, and is NA
#    for the trios whose majority vote supports no edge. Neither is dropped -- the first is
#    the provenance, the second is the quantity the phrase "weakest edge" actually means.

library(writexl)

in.file  <- "bioinfo_revision/reports/gtex_bootstrapping/gtex_mrgn_bootstrap.csv"
out.dir  <- "bioinfo_revision/Final_tables"
out.file <- file.path(out.dir, "GTEx_bootstrap_summary.xlsx")

d <- utils::read.csv(in.file, stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------------------
# label sets
# ---------------------------------------------------------------------------------------
#
# Fixed levels, so a model that never occurs still gets a row of zeros. The tables are meant
# to be read against the simulation tables and against each other, and a silently dropped
# row would not line up. M2.1 is the case in point: it does not occur in this run at all.
MODEL.LEVELS <- c("M0.1", "M0.2", "M1.1", "M1.2", "M2.1", "M2.2", "M3", "M4", "Other")

unknown <- setdiff(c(d$model, d$boot.model), c(MODEL.LEVELS, NA))
if (length(unknown)) {
    stop("unrecognised model label(s): ", paste(unknown, collapse = ", "))
}

EDGE.PROB.COLS <- c(`V1 -> T1` = "boot.p.V1T1",
                    `T1 - T2 (unoriented)` = "boot.p.T1T2",
                    `V1 -> T2` = "boot.p.V1T2")

# Property 1 above. Checked rather than assumed.
if (!identical(d$boot.p.T1T2, d$boot.p.T2T1)) {
    stop("boot.p.T1T2 and boot.p.T2T1 differ in ", sum(d$boot.p.T1T2 != d$boot.p.T2T1),
         " trios. They have always been one unoriented indicator reported twice; if MRGN",
         " now splits them, edge_support must report both rather than collapsing them.")
}

# Property 2 above: the weakest SUPPORTED EDGE, over the four edge indicators only. NA when
# the majority vote supports no edge -- a trio with an empty skeleton has no weakest edge,
# and 0 would be a different and false claim.
edge.p <- as.matrix(d[, c("boot.p.V1T1", "boot.p.T1T2", "boot.p.V1T2", "boot.p.T2T1")])
min.supported.edge <- apply(edge.p, 1, function(r) {
    r <- r[r >= 0.5]
    if (length(r)) min(r) else NA_real_
})
d$min.supported.edge <- min.supported.edge

# How often the recorded column is set by an interaction term rather than an edge. Quoted in
# run_summary, so the number in the README and the number in the sheet cannot drift apart.
n.no.edge.supported <- sum(is.na(min.supported.edge))
n.interaction.floor <- sum(!is.na(min.supported.edge) &
                           min.supported.edge > d$boot.min.edge.prob + 1e-9)

# ---------------------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------------------

# Vectorised over the counts: several sheets take a whole column of them, against either one
# shared denominator or a column of them. A zero denominator gives NA rather than NaN -- "no
# trios to take a share of" is not the same statement as 0%.
#
# Written out rather than as ifelse(), which takes its length from the TEST: a scalar `n`
# there silently collapses a six-element result to one.
pct <- function(x, n) {
    n <- rep_len(n, length(x))
    out <- rep(NA_real_, length(x))
    ok <- !is.na(n) & n != 0
    out[ok] <- 100 * x[ok] / n[ok]
    out
}

# Median and quartiles as three columns, so a sheet can carry a distribution without a plot.
quart <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(c(q1 = NA_real_, median = NA_real_, q3 = NA_real_))
    q <- stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
    c(q1 = q[1], median = q[2], q3 = q[3])
}

n.trios <- nrow(d)

# ---------------------------------------------------------------------------------------
# sheet: run_summary
# ---------------------------------------------------------------------------------------

t.point <- quart(d$time.seconds)
t.boot  <- quart(d$bootstrap.time.seconds)
cov.q   <- quart(d$n.covariates)
pc.q    <- quart(d$n.pcs)
ind.q   <- quart(d$boot.min.edge.prob)
edge.q  <- quart(d$min.supported.edge)

run.summary <- data.frame(
    quantity = c(
        "Tissue",
        "Trios",
        "Samples per trio",
        "Covariates per trio (Q1 / median / Q3)",
        "  of which selected PCs (Q1 / median / Q3)",
        "Resamples requested per trio",
        "Resamples used (min)",
        "Resamples dropped: trios with any",
        "Resamples dropped: largest for one trio",
        "Point inference failures",
        "Bootstrap failures",
        "Bootstrap majority agrees with point call",
        "Bootstrap majority differs from point call",
        "Trios whose supported vector contains no edge",
        "Trios whose min.supported.indicator is an interaction term, not an edge",
        "Weakest supported indicator (Q1 / median / Q3)",
        "Weakest supported edge (Q1 / median / Q3)",
        "Point inference seconds per trio (Q1 / median / Q3)",
        "Bootstrap seconds per trio (Q1 / median / Q3)",
        "Total bootstrap CPU time (hours)"),
    value = c(
        paste(unique(d$tissue), collapse = ", "),
        format(n.trios, big.mark = ","),
        paste(unique(d$n.samples), collapse = ", "),
        sprintf("%g / %g / %g", cov.q[1], cov.q[2], cov.q[3]),
        sprintf("%g / %g / %g", pc.q[1], pc.q[2], pc.q[3]),
        format(unique(d$boot.n.requested), big.mark = ","),
        format(min(d$boot.n.used), big.mark = ","),
        sprintf("%d (%.1f%%)", sum(d$boot.n.dropped > 0),
                pct(sum(d$boot.n.dropped > 0), n.trios)),
        as.character(max(d$boot.n.dropped)),
        as.character(sum(!is.na(d$error))),
        as.character(sum(!is.na(d$bootstrap.error))),
        sprintf("%d (%.1f%%)", sum(d$boot.agrees), pct(sum(d$boot.agrees), n.trios)),
        sprintf("%d (%.1f%%)", sum(!d$boot.agrees), pct(sum(!d$boot.agrees), n.trios)),
        sprintf("%d (%.1f%%)", n.no.edge.supported, pct(n.no.edge.supported, n.trios)),
        sprintf("%d (%.1f%%)", n.interaction.floor + n.no.edge.supported,
                pct(n.interaction.floor + n.no.edge.supported, n.trios)),
        sprintf("%.3f / %.3f / %.3f", ind.q[1], ind.q[2], ind.q[3]),
        sprintf("%.3f / %.3f / %.3f", edge.q[1], edge.q[2], edge.q[3]),
        sprintf("%.3f / %.3f / %.3f", t.point[1], t.point[2], t.point[3]),
        sprintf("%.1f / %.1f / %.1f", t.boot[1], t.boot[2], t.boot[3]),
        sprintf("%.1f", sum(d$bootstrap.time.seconds, na.rm = TRUE) / 3600)),
    stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------------------
# sheet: model_calls
# ---------------------------------------------------------------------------------------
#
# One row per model label. `point` and `bootstrap` are the two marginals of the cross-tab
# below; `stability` is the share of the POINT trios the bootstrap still calls that model,
# so it is a row rate of that table and reads per label rather than overall.

point.n <- as.integer(table(factor(d$model, levels = MODEL.LEVELS)))
boot.n  <- as.integer(table(factor(d$boot.model, levels = MODEL.LEVELS)))

model.calls <- do.call(rbind, lapply(seq_along(MODEL.LEVELS), function(i) {
    m <- MODEL.LEVELS[i]
    rows <- d[d$model == m, , drop = FALSE]
    q <- quart(rows$min.supported.edge)
    data.frame(
        model               = m,
        point.n             = point.n[i],
        point.pct           = pct(point.n[i], n.trios),
        bootstrap.n         = boot.n[i],
        bootstrap.pct       = pct(boot.n[i], n.trios),
        stability.n         = sum(rows$boot.agrees),
        stability.pct       = pct(sum(rows$boot.agrees), nrow(rows)),
        min.supported.edge.q1         = q[["q1"]],
        min.supported.edge.median     = q[["median"]],
        min.supported.edge.q3         = q[["q3"]],
        min.supported.edge.ge.0.8.pct = pct(sum(rows$min.supported.edge >= 0.8, na.rm = TRUE),
                                  nrow(rows)),
        stringsAsFactors = FALSE)
}))
model.calls <- rbind(model.calls, data.frame(
    model               = "All trios",
    point.n             = n.trios,
    point.pct           = 100,
    bootstrap.n         = n.trios,
    bootstrap.pct       = 100,
    stability.n         = sum(d$boot.agrees),
    stability.pct       = pct(sum(d$boot.agrees), n.trios),
    min.supported.edge.q1         = quart(d$min.supported.edge)[["q1"]],
    min.supported.edge.median     = quart(d$min.supported.edge)[["median"]],
    min.supported.edge.q3         = quart(d$min.supported.edge)[["q3"]],
    min.supported.edge.ge.0.8.pct = pct(sum(d$min.supported.edge >= 0.8, na.rm = TRUE), n.trios),
    stringsAsFactors = FALSE))

# ---------------------------------------------------------------------------------------
# sheet: point_vs_bootstrap
# ---------------------------------------------------------------------------------------
#
# Rows are the point call, columns the bootstrap-supported call, matching the row-is-what-
# the-method-said convention of the simulation matrices. The diagonal is agreement. There is
# no Precision column and no Recall row, because neither margin is truth: what stands in
# their place is `stable.pct`, the diagonal over the row total.

tab <- table(factor(d$model,      levels = MODEL.LEVELS),
             factor(d$boot.model, levels = MODEL.LEVELS))
tab <- matrix(as.integer(tab), nrow = length(MODEL.LEVELS),
              dimnames = list(MODEL.LEVELS, MODEL.LEVELS))

row.total <- rowSums(tab)
point.vs.bootstrap <- data.frame(
    point.call = c(MODEL.LEVELS, "Total"),
    rbind(tab, Total = colSums(tab)),
    Total      = c(row.total, sum(tab)),
    stable.pct = c(ifelse(row.total == 0, NA_real_, 100 * diag(tab) / row.total),
                   pct(sum(diag(tab)), sum(tab))),
    check.names = FALSE, stringsAsFactors = FALSE)
rownames(point.vs.bootstrap) <- NULL

# ---------------------------------------------------------------------------------------
# sheet: disagreements
# ---------------------------------------------------------------------------------------
#
# The off-diagonal of the table above, as a ranked list. Same numbers, but the cross-tab is
# mostly zeros and the question "where does resampling actually move a call" is easier to
# answer from 20-odd rows than from an 81-cell grid.

off <- which(tab > 0 & row(tab) != col(tab), arr.ind = TRUE)
disagreements <- data.frame(
    point.call     = MODEL.LEVELS[off[, "row"]],
    bootstrap.call = MODEL.LEVELS[off[, "col"]],
    n              = tab[off],
    pct.of.point.call = pct(tab[off], row.total[off[, "row"]]),
    pct.of.all.disagreements = pct(tab[off], sum(!d$boot.agrees)),
    stringsAsFactors = FALSE)
disagreements <- disagreements[order(-disagreements$n), ]
rownames(disagreements) <- NULL

# ---------------------------------------------------------------------------------------
# sheet: edge_support
# ---------------------------------------------------------------------------------------
#
# The resample mass behind each of the three distinct edges, over all trios. The thresholds
# are read two ways on purpose: `pct.ge.*` is how often an edge is supported, `pct.le.0.05`
# how often it is confidently ABSENT. A bootstrap that is uninformative would pile up in the
# middle and show low numbers in both.

edge.support <- do.call(rbind, lapply(seq_along(EDGE.PROB.COLS), function(i) {
    x <- d[[EDGE.PROB.COLS[i]]]
    q <- quart(x)
    data.frame(
        edge        = names(EDGE.PROB.COLS)[i],
        mean        = mean(x, na.rm = TRUE),
        q1          = q[["q1"]],
        median      = q[["median"]],
        q3          = q[["q3"]],
        pct.le.0.05 = pct(sum(x <= 0.05, na.rm = TRUE), n.trios),
        pct.ge.0.50 = pct(sum(x >= 0.50, na.rm = TRUE), n.trios),
        pct.ge.0.80 = pct(sum(x >= 0.80, na.rm = TRUE), n.trios),
        pct.ge.0.95 = pct(sum(x >= 0.95, na.rm = TRUE), n.trios),
        stringsAsFactors = FALSE)
}))
rownames(edge.support) <- NULL

# ---------------------------------------------------------------------------------------
# sheet: confidence_bands
# ---------------------------------------------------------------------------------------
#
# How the per-trio confidence floor is distributed. Bands start at 0.50 because both columns
# are a minimum over indicators that cleared 0.50 by construction -- nothing below that can
# occur, and a band there would read as an empty region of a range rather than as an
# impossible one. The NA row is the 21 trios with no supported edge; they have a
# min.supported.indicator but no min.supported.edge.

BREAKS <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0000001)
BAND.LABELS <- c("0.50 - 0.60", "0.60 - 0.70", "0.70 - 0.80", "0.80 - 0.90",
                 "0.90 - 0.95", "0.95 - 1.00")

band.of <- function(x) cut(x, breaks = BREAKS, labels = BAND.LABELS, right = FALSE)

band.counts <- function(x) as.integer(table(band.of(x)))
ind.n  <- band.counts(d$boot.min.edge.prob)
edge.n <- band.counts(d$min.supported.edge)

confidence.bands <- data.frame(
    band = c(BAND.LABELS, "no supported edge", "Total"),
    min.supported.indicator.n   = c(ind.n, 0L, n.trios),
    min.supported.indicator.pct = c(pct(ind.n, n.trios), 0, 100),
    min.supported.edge.n        = c(edge.n, n.no.edge.supported, n.trios),
    min.supported.edge.pct      = c(pct(edge.n, n.trios),
                                    pct(n.no.edge.supported, n.trios), 100),
    stringsAsFactors = FALSE)

stopifnot(sum(ind.n) == n.trios,
          sum(edge.n) + n.no.edge.supported == n.trios)

# ---------------------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------------------

stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

readme <- data.frame(
    field = c(
        "What this workbook is",
        "No ground truth",
        "Sheet: run_summary",
        "Sheet: model_calls",
        "Sheet: point_vs_bootstrap",
        "Sheet: disagreements",
        "Sheet: edge_support",
        "Sheet: confidence_bands",
        "Point call vs bootstrap call",
        "stability.pct",
        "T1-T2 is unoriented",
        "min.supported.indicator",
        "min.supported.edge",
        "Model labels",
        "Dropped resamples",
        "Generated",
        "Source",
        "Script"),
    value = c(
        paste("A summary of MRGN's edge-probability bootstrap over the", n.trios,
              "GTEx Whole Blood trios: what MRGN called, what resampling supports, and how",
              "much resample mass sits behind each edge."),
        paste("This is real data. Nothing here is scored against truth and there is no",
              "confusion matrix: the bootstrap measures STABILITY under resampling, not",
              "correctness. A call that is perfectly stable can still be wrong."),
        "Run shape and headline totals: trios, covariates, resamples, failures, timing.",
        paste("One row per model label: how many trios the point fit gave it, how many the",
              "bootstrap majority gives it, and how often the two agree for that label."),
        paste("The central table. Point call (rows) x bootstrap-supported call (columns).",
              "The diagonal is agreement."),
        "The off-diagonal of that table as a ranked list -- where resampling moves a call.",
        "Per-edge resample support across all trios, with absence and support thresholds.",
        "How the per-trio confidence floor is distributed.",
        paste("The point call is infer.trio() on the full 670 samples. The bootstrap call is",
              "MRGN::class.vec() applied to the majority vote of the six indicators across",
              "the usable resamples -- a second, independent label for the same trio."),
        paste("Of the trios given a model by the point fit, the share the bootstrap majority",
              "still gives that model. It replaces recall/precision, which need a truth axis",
              "this data does not have."),
        paste("infer.trio() returns b12 and b22 as one indicator for T1-T2 ADJACENCY, not",
              "two orientation-specific ones; they are identical in every trio here.",
              "boot.p.T1T2 and boot.p.T2T1 are therefore one number written twice, and",
              "edge_support reports the edge once. Orientation is decided by the two",
              "interaction terms, not by these."),
        paste("The boot.min.edge.prob column as recorded: the smallest mean among ALL SIX",
              "indicators the majority vote supports -- including the V1:T1 and V1:T2",
              "INTERACTION terms, which are not edges. In", n.interaction.floor +
              n.no.edge.supported, "of the", n.trios, "trios the minimum is one of those",
              "two rather than any edge, so despite the name it is a floor on the whole",
              "supported indicator vector, not on the skeleton."),
        paste("Recomputed here over the four edge indicators only, which is what 'weakest",
              "supported edge' actually means. NA for the", n.no.edge.supported,
              "trios whose majority vote supports no edge at all -- those have no weakest",
              "edge, and 0 would be a different and false claim."),
        paste("MRGN's eight topologies plus Other. The .1 / .2 suffix is which gene the",
              "variant acts on directly: .1 through the cis gene T1, .2 through the trans",
              "gene T2. M2.1 does not occur in this run and keeps a zero row."),
        paste("A resample that loses every copy of the minor allele leaves V1 constant and",
              "says nothing about the V1 edges, so it is dropped rather than averaged in.",
              "At n = 670 this is rare; the counts are in run_summary."),
        stamp,
        in.file,
        "bioinfo_revision/gtex_bootstrapping/summarise_gtex_bootstrap.R"),
    stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------------------
# write
# ---------------------------------------------------------------------------------------

if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

writexl::write_xlsx(
    list(README            = readme,
         run_summary       = run.summary,
         model_calls       = model.calls,
         point_vs_bootstrap = point.vs.bootstrap,
         disagreements     = disagreements,
         edge_support      = edge.support,
         confidence_bands  = confidence.bands),
    out.file)

cat(sprintf("%d trios | %.1f%% bootstrap agreement | median weakest supported edge %.3f\n",
            n.trios, pct(sum(d$boot.agrees), n.trios),
            stats::median(d$min.supported.edge, na.rm = TRUE)))
cat(sprintf("wrote %s (7 sheets)\n", out.file))
