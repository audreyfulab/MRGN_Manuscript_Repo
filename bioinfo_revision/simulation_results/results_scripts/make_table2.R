# Table 2 rebuilt on the revised simulation: four methods side by side.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/make_table2.R
#
# Reproduces the layout of Table 2 of MRGN_v8.pdf -- rows are the five generating models
# plus the T1-T2 edge, columns are method-arm blocks each carrying Recall and Precision --
# extended to the four methods of the revision. Reads tables/confusion_counts_long.csv, so
# run make_all_tables.R first.
#
# GMAC CANNOT FILL THE MODEL ROWS. It returns a mediation call (Cis Mediated / Trans
# Mediated / No Mediation / Undirected), which is not one of the five models, so it
# contributes the T1-T2 row alone. The original Table 2 handles it the same way -- named in
# the caption and marked with an asterisk on the T1-T2 row, with no column of its own.
#
# MR-GGI now DOES fill them, from mrggi.model. That column is built by MRGN::class.vec()
# from MR-GGI's own pairwise estimates plus the instrument-gene and marginal tests, and it
# is reported because it FAILS: M0 and M2 are structurally unreachable, because with one
# instrument the Wald-ratio p reduces to the instrument->outcome t-test and the causal tests
# are therefore not independent of the marginals. See MRGGI_METHODS.md section 5.2.
#
# MRPC COVERS THREE OF THE FIVE SAMPLE SIZES. Its n = 670 and n = 1000 groups were dropped
# (see inference_config.R), so two tables are written: one over all 1,500 trios with MRPC
# flagged, and one restricting every method to the 900 trios at n <= 300 where all four are
# directly comparable. Reading them together is what separates a real MRPC gap from a
# coverage artefact.

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

REPORT.DIR <- "bioinfo_revision/reports"
dir.create(REPORT.DIR, recursive = TRUE, showWarnings = FALSE)
OUT.MD  <- file.path(REPORT.DIR, "table2_inference_summary.md")
OUT.CSV <- file.path(REPORT.DIR, "table2_inference_summary.csv")

counts <- utils::read.csv(file.path(tables.dir, "confusion_counts_long.csv"),
                          stringsAsFactors = FALSE)

# The columns, in the order Table 2 uses: the oracle first, then the realistic arms, then
# the two edge-only methods.
COLUMNS <- list(
    list(method = "mrgn",  arm = "truth", model.level = "model", edge.level = "edge",
         label = "MRGN (True Confs.)", oracle = TRUE),
    list(method = "mrgn",  arm = "CSq",   model.level = "model", edge.level = "edge",
         label = "MRGN (CS-q)", oracle = FALSE),
    list(method = "mrgn",  arm = "CSa",   model.level = "model", edge.level = "edge",
         label = "MRGN (CS-α)", oracle = FALSE),
    list(method = "mrpc",  arm = "CSq",   model.level = "model", edge.level = "edge",
         label = "MRPC (CS-q)", oracle = FALSE),
    list(method = "gmac",  arm = "gmac",  model.level = NA,      edge.level = "edge",
         label = "GMAC", oracle = FALSE),
    # MR-GGI is shown on its ARM-INVARIANT columns, both of which are reported under the arm
    # label "mrggi": the raw-p edge call, which is the one comparable with GMAC and MRGN, and
    # the model call, which class.vec() builds from the same pairwise estimates. Its four
    # covariate arms move only the FDR-adjusted edge call, so an arm-specific column here
    # would differ from its neighbours in the correction applied rather than in the method.
    # The per-arm edge.fdr figures are in INFERENCE_PERFORMANCE.md.
    list(method = "mrggi", arm = "CSq",   model.level = "model", edge.level = "edge.fdr",
         label = "MR-GGI (CS-q)", oracle = FALSE))

MODEL.ROWS <- c(M0 = "Null Model (**M0**)",
                M1 = "Mediation (**M1**)",
                M2 = "V-structure (**M2**)",
                M3 = "Cond. Indep. (**M3**)",
                M4 = "Fully Connect. (**M4**)")

# Pools the per-sample-size matrices into one. Every cell is a raw count, so summing is
# exact; `sizes` restricts which groups are pooled, which is the whole mechanism behind the
# matched table.
pool.matrix <- function(method, arm, level, sizes) {
    x <- counts[counts$method == method & counts$arm == arm & counts$level == level &
                counts$effect_size == "all" & counts$sample_size %in% sizes, , drop = FALSE]
    if (nrow(x) == 0) return(NULL)
    matrix.from.long(stats::aggregate(n ~ truth + predicted, data = x, FUN = sum),
                     pred.levels.for(method, level))
}

# confusion() puts truth in the rows and the inferred label in the columns, so for model M:
#   recall    = m[M, M] / sum(m[M, ])    of the trios generated under M, the share called M
#   precision = m[M, M] / sum(m[, M])    of the trios called M, the share generated under M
model.scores <- function(m) {
    lapply(TRUTH.LEVELS, function(t) {
        if (is.null(m) || !t %in% rownames(m) || !t %in% colnames(m))
            return(c(recall = NA_real_, prec = NA_real_))
        hit <- m[t, t]
        c(recall = if (sum(m[t, ]) > 0) hit / sum(m[t, ]) else NA_real_,
          prec   = if (sum(m[, t]) > 0) hit / sum(m[, t]) else NA_real_)
    })
}

# The T1-T2 row scores the edge-PRESENT class, matching Table 2 -- where MRGN's oracle
# column reads recall 1.000 against precision 0.767, i.e. it finds nearly every real edge
# and pays for it in false positives.
edge.row <- function(method, arm, level, sizes) {
    m <- pool.matrix(method, arm, level, sizes)
    if (is.null(m)) return(c(recall = NA_real_, prec = NA_real_))
    nc <- if (method == "mrpc") c("Other", "Failed")
          else if (method == "mrggi") c("Weak instrument", "Screened out")
          else if (method == "gmac") "" else "Other"
    s <- edge.scores(m, no.call.level = nc)
    c(recall = s$present.recall, prec = s$present.prec)
}

# ---------------------------------------------------------------------------------------
# assemble one table
# ---------------------------------------------------------------------------------------
build <- function(sizes) {
    cells <- lapply(COLUMNS, function(cc) {
        mm <- if (is.na(cc$model.level)) NULL
              else pool.matrix(cc$method, cc$arm, cc$model.level, sizes)
        ms <- model.scores(mm)
        names(ms) <- TRUTH.LEVELS
        list(model = ms,
             edge  = edge.row(cc$method, cc$arm, cc$edge.level, sizes),
             n     = if (is.null(mm)) {
                         em <- pool.matrix(cc$method, cc$arm, cc$edge.level, sizes)
                         if (is.null(em)) NA_integer_ else sum(em)
                     } else sum(mm))
    })
    names(cells) <- vapply(COLUMNS, function(cc) cc$label, character(1))
    cells
}

# Table 2 highlights the best recall and the best precision in each row, excluding the
# oracle column -- the oracle is a ceiling, not a competitor, so including it would make
# every row bold in the same place and say nothing.
fmt.row <- function(vals, digits = 3) {
    competitors <- !vapply(COLUMNS, function(cc) cc$oracle, logical(1))
    out <- character(0)
    for (metric in c("recall", "prec")) {
        v <- vapply(vals, function(x) unname(x[[metric]]), numeric(1))
        best <- suppressWarnings(max(v[competitors], na.rm = TRUE))
        f <- ifelse(is.na(v), "--", formatC(v, format = "f", digits = digits))
        mark <- competitors & !is.na(v) & abs(v - best) < 1e-9
        f[mark] <- paste0("**", f[mark], "**")
        out <- rbind(out, f)
    }
    # interleave recall, precision per column
    as.vector(rbind(out[1, ], out[2, ]))
}

md.table <- function(cells) {
    hdr <- unlist(lapply(names(cells), function(l)
        c(paste0(l, "<br>Recall"), paste0(l, "<br>Prec."))))
    lines <- c(paste("|  | ", paste(hdr, collapse = " | "), "|"),
               paste("| ---: |", paste(rep("---:", length(hdr)), collapse = " | "), "|"))
    for (t in TRUTH.LEVELS) {
        vals <- lapply(cells, function(cc) cc$model[[t]])
        lines <- c(lines, paste("|", MODEL.ROWS[[t]], "|",
                                paste(fmt.row(vals), collapse = " | "), "|"))
    }
    vals <- lapply(cells, function(cc) cc$edge)
    c(lines, paste("| **T1 − T2 edge** |", paste(fmt.row(vals), collapse = " | "), "|"))
}

csv.rows <- function(cells, scope) {
    do.call(rbind, lapply(names(cells), function(l) {
        cc <- cells[[l]]
        rbind(
            do.call(rbind, lapply(TRUTH.LEVELS, function(t) data.frame(
                scope = scope, column = l, row = t,
                recall = unname(cc$model[[t]]["recall"]),
                precision = unname(cc$model[[t]]["prec"]),
                trios = cc$n, stringsAsFactors = FALSE))),
            data.frame(scope = scope, column = l, row = "T1-T2 edge",
                       recall = unname(cc$edge["recall"]),
                       precision = unname(cc$edge["prec"]),
                       trios = cc$n, stringsAsFactors = FALSE))
    }))
}

ALL     <- SAMPLE.SIZES
MATCHED <- SAMPLE.SIZES[SAMPLE.SIZES <= 300]

cells.all     <- build(ALL)
cells.matched <- build(MATCHED)

n.of <- function(cells) vapply(cells, function(cc) cc$n, numeric(1))

# the arm-invariant MR-GGI figures, for the note
ggi.raw.all     <- edge.row("mrggi", "mrggi", "edge", ALL)
ggi.raw.matched <- edge.row("mrggi", "mrggi", "edge", MATCHED)

lines <- c(
"# Overall inference performance: four methods side by side",
"",
paste("Generated by `results_scripts/make_table2.R` from `tables/confusion_counts_long.csv`.",
      "Do not edit by hand; rerun the script."),
"",
paste("Rebuilds the layout of **Table 2 of `MRGN_v8.pdf`** on the revised simulation, with",
      "MR-GGI added. Rows are the five generating models plus the T1−T2 edge; each method",
      "column carries Recall and Precision."),
"",
"**Definitions.** For a generating model `M`:",
"",
"```",
"recall(M)    = trios generated under M that were labelled M  /  trios generated under M",
"precision(M) = trios labelled M that were generated under M  /  trios labelled M",
"```",
"",
paste("The `T1 − T2 edge` row scores the **edge-present** class: recall is the share of",
      "genuine T1−T2 edges recovered, precision the share of called edges that are real.",
      "Model labels are coarse -- `M0.1` and `M0.2` both count as `M0`."),
"",
paste("**Bold marks the best recall and the best precision in each row**, excluding the",
      "MRGN oracle column. The oracle is a ceiling rather than a competitor: it is MRGN",
      "handed the true confounder set, shown to separate the cost of confounder *selection*",
      "from the cost of the method."),
"",
"---",
"",
"## Table 2a. All available trios",
"",
md.table(cells.all),
"",
paste(sprintf("**Trios per column:** %s.",
              paste(sprintf("%s %d", names(n.of(cells.all)), n.of(cells.all)),
                    collapse = "; "))),
"",
paste("**MRPC is scored on 900 trios against the others' 1,500.** Its n = 670 and n = 1000",
      "groups were not run -- MRPC times out on 61% and 75% of trios there even in its cheap",
      "arm, at ~10 h per group (`inference_config.R`). Its column is therefore not strictly",
      "comparable to its neighbours here; Table 2b removes that difference."),
"",
"---",
"",
"## Table 2b. Matched — all four methods on the same 900 trios (n ≤ 300)",
"",
md.table(cells.matched),
"",
paste(sprintf("**Trios per column:** %s.",
              paste(sprintf("%s %d", names(n.of(cells.matched)), n.of(cells.matched)),
                    collapse = "; "))),
"",
paste("Every column is now scored on identical trios, so differences between columns are",
      "differences between methods. Note this is the harder end of the sample-size range for",
      "all four -- confounder selection has no power at n = 50 -- so the absolute levels are",
      "below Table 2a for every method that lost groups."),
"",
"---",
"",
"## Notes",
"",
paste("**GMAC has no model rows, and this is structural.** It returns a mediation call --",
      "Cis Mediated / Trans Mediated / No Mediation / Undirected -- which is not one of the",
      "five models, so it cannot be scored on them and is shown as `--`. The original",
      "Table 2 treats it identically: named in the caption and marked on the T1−T2 row, but",
      "with no column."),
"",
paste("**MR-GGI's model rows are reported because they fail.** `mrggi.model` is built by",
      "`MRGN::class.vec()` from MR-GGI's own pairwise estimates together with the",
      "instrument-gene and marginal tests. With a single instrument the Wald-ratio p-value",
      "reduces to the instrument→outcome t-test, so the two causal tests are not independent",
      "of the two marginals: the six-vector `class.vec()` receives carries four distinct",
      "tests, and **M0 and M2 are unreachable**. Read the model row for MR-GGI as a",
      "measurement of that limit, not as a competitive score. `MRGGI_METHODS.md` §5.2."),
"",
paste("**MR-GGI's arms are not confounder adjustments.** `MRggi()` has no covariate",
      "argument; its estimator is strictly pairwise, so the raw-p T1−T2 call is *identical*",
      "in all four arms and only the FDR-adjusted call varies. The `MR-GGI (CS-q)` column",
      "above is the FDR-adjusted call. The arm-invariant raw-p call -- the one directly",
      sprintf("comparable with GMAC and MRGN -- gives recall %.3f / precision %.3f over all",
              ggi.raw.all["recall"], ggi.raw.all["prec"]),
      sprintf("trios and %.3f / %.3f matched.",
              ggi.raw.matched["recall"], ggi.raw.matched["prec"])),
"",
paste("**MRGN's oracle arm is a genuine ceiling; MRPC's is not.** For MRGN, extra covariates",
      "cost residual degrees of freedom, so the true-confounder arm bounds what selection",
      "could achieve. For MRPC the same arm is a handicap -- CS-q beats it at every sample",
      "size, because MRPC runs a PC algorithm whose cost grows with the node count whether",
      "or not the covariates are the right ones. That is why no MRPC oracle column is shown",
      "here; see `INFERENCE_PERFORMANCE.md` §3.1 and `METHODS.md` Table 12c."),
"",
paste("**A no-call is not a wrong call.** Trios landing in `Other` (no topology matched),",
      "`Failed` (MRPC timed out) or `Weak instrument` (MR-GGI's first-stage F gate) depress",
      "recall without inflating any other model's precision. Per-method no-call rates are in",
      "`INFERENCE_PERFORMANCE.md` §2."),
"",
paste("Machine-readable form: `table2_inference_summary.csv`. Full matrices behind every",
      "cell: `INFERENCE_PERFORMANCE.md` §4 and `tables/confusion_counts_long.csv`."),
"")

writeLines(lines, OUT.MD, useBytes = TRUE)
utils::write.csv(rbind(csv.rows(cells.all, "all"), csv.rows(cells.matched, "matched")),
                 OUT.CSV, row.names = FALSE)

cat(sprintf("wrote %s | %d lines\n", OUT.MD, length(lines)))
cat(sprintf("wrote %s\n", OUT.CSV))
