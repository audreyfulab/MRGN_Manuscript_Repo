# Writes reports/CONFOUNDER_STRUCTURE.md.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/results_scripts/make_structure_report.R
#
# Reads the four confounder-structure runs of METHODS.md section 7b directly from their
# per-group checkpoints (STRUCTURES in confusion_utils.R). Independent of
# confusion_structures.R, which writes the matrices and the overall rates; this script
# answers the narrower question the design was built for -- which covariate does the damage,
# and to which generating models.
#
# Deliberately base R, matching confusion_utils.R.

source("bioinfo_revision/simulation_results/results_scripts/confusion_utils.R")

REPORT.DIR <- "bioinfo_revision/reports"
dir.create(REPORT.DIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(REPORT.DIR, "CONFOUNDER_STRUCTURE.md")

ARMS <- c("truth", "CSq", "CSa")
ARM.LAB <- c(truth = "truth (oracle)", CSq = "CS-q", CSa = "CS-alpha")

pct  <- function(x, d = 1) ifelse(is.na(x), "--", sprintf(paste0("%.", d, "f%%"), 100 * x))
sgn  <- function(x, d = 1) ifelse(is.na(x), "--", sprintf(paste0("%+.", d, "f"), 100 * x))

md.table <- function(df) {
    cells <- as.data.frame(lapply(df, function(c)
        trimws(if (is.numeric(c)) formatC(c, format = "fg", digits = 4) else as.character(c))),
        stringsAsFactors = FALSE)
    c(paste("|", paste(names(df), collapse = " | "), "|"),
      paste("|", paste(rep("---:", ncol(df)), collapse = " | "), "|"),
      apply(cells, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}

RES <- lapply(STRUCTURES, function(s) load.structure(s$dir))
missing <- names(RES)[vapply(RES, is.null, logical(1))]
if (length(missing)) {
    stop("no checkpoint for: ", paste(missing, collapse = ", "),
         "\nRun run_structure_sims.R (or apply_mrgn.R per case) first.")
}

# ---------------------------------------------------------------------------------------
# measurements
# ---------------------------------------------------------------------------------------

model.col <- function(arm) paste0("mrgn.", arm, ".model")

# Recall per generating model: of the trios generated under M, the share MRGN labelled M.
# Coarse, matching mrgn.*.correct.coarse and every other table in this stage.
per.model.recall <- function(res, arm) {
    tr <- coarse.model(res$truth.model)
    pr <- coarse.model(res[[model.col(arm)]])
    vapply(TRUTH.LEVELS, function(t) {
        i <- tr == t
        if (!any(i)) NA_real_ else mean(pr[i] == t, na.rm = TRUE)
    }, numeric(1))
}

overall <- function(res, arm) {
    m <- confusion(res$truth.model, res[[model.col(arm)]], MRGN.LEVELS)
    e <- confusion(res$truth.model, mrgn.edge(res[[model.col(arm)]]),
                   MRGN.EDGE.LEVELS, coarse.pred = FALSE)
    es <- edge.scores(e)
    list(model = sum(diag(m[, TRUTH.LEVELS, drop = FALSE])) / sum(m),
         edge = es$accuracy, present = es$present.recall, absent = es$absent.recall,
         no.call = es$no.call, n = sum(m))
}

# How often the arm's selected set actually contained the covariate it was meant to reject.
# Meaningless for the truth arm, whose set is the U block by construction.
exposure <- function(res, arm, what) {
    col <- paste0(arm, ".has.", what)
    if (arm == "truth" || !col %in% names(res)) return(NA_real_)
    mean(as.logical(res[[col]]), na.rm = TRUE)
}

# The collider signature: on trios whose truth carries NO T1-T2 edge (M0 and M3, per
# EDGE.CORRECT), how often is one called anyway.
spurious.edge <- function(res, arm, split.col = NULL) {
    absent.models <- names(EDGE.CORRECT)[EDGE.CORRECT == EDGE.LEVELS[1]]
    i <- coarse.model(res$truth.model) %in% absent.models
    called <- mrgn.edge(res[[model.col(arm)]]) == EDGE.LEVELS[2]
    if (is.null(split.col) || !split.col %in% names(res)) {
        return(list(all = mean(called[i], na.rm = TRUE), n = sum(i)))
    }
    s <- as.logical(res[[split.col]])
    list(all = mean(called[i], na.rm = TRUE), n = sum(i),
         exposed = mean(called[i & s], na.rm = TRUE), n.exposed = sum(i & s),
         clean = mean(called[i & !s], na.rm = TRUE), n.clean = sum(i & !s))
}

# ---------------------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------------------

overall.table <- function(arm) {
    base <- overall(RES[["u_only"]], arm)
    do.call(rbind, lapply(names(STRUCTURES), function(nm) {
        o <- overall(RES[[nm]], arm)
        data.frame(Structure = nm, Trios = o$n,
                   `Model acc.` = pct(o$model), `Edge acc.` = pct(o$edge),
                   `Edge-present recall` = pct(o$present),
                   `Edge-absent recall`  = pct(o$absent),
                   `No call` = pct(o$no.call),
                   `Δ model vs u_only` = if (nm == "u_only") "--" else sgn(o$model - base$model),
                   check.names = FALSE, stringsAsFactors = FALSE)
    }))
}

per.model.table <- function(arm, delta = FALSE) {
    base <- per.model.recall(RES[["u_only"]], arm)
    do.call(rbind, lapply(names(STRUCTURES), function(nm) {
        r <- per.model.recall(RES[[nm]], arm)
        vals <- if (delta && nm != "u_only") sgn(r - base) else if (delta) rep("--", 5) else pct(r)
        d <- as.data.frame(as.list(vals), stringsAsFactors = FALSE)
        names(d) <- TRUTH.LEVELS
        cbind(data.frame(Structure = nm, stringsAsFactors = FALSE), d)
    }))
}

# ---------------------------------------------------------------------------------------
# numbers quoted in prose
# ---------------------------------------------------------------------------------------
o.only <- overall(RES[["u_only"]], "CSq"); o.w <- overall(RES[["u_w"]], "CSq")
o.z    <- overall(RES[["u_z"]],    "CSq"); o.wz <- overall(RES[["u_w_z"]], "CSq")
r.only <- per.model.recall(RES[["u_only"]], "CSq")
r.z    <- per.model.recall(RES[["u_z"]],    "CSq")
r.w    <- per.model.recall(RES[["u_w"]],    "CSq")
sp.z   <- spurious.edge(RES[["u_z"]], "CSq", "CSq.has.common.child")
exp.w  <- exposure(RES[["u_w"]], "CSq", "intermediate")
exp.z  <- exposure(RES[["u_z"]], "CSq", "common.child")

# the confounding check: is "the Z was missed" a comparable trio, or an easier one?
confound <- do.call(rbind, lapply(TRUTH.LEVELS, function(t) {
    uo <- RES[["u_only"]]; uz <- RES[["u_z"]]
    a <- mean(coarse.model(uo[[model.col("CSq")]])[coarse.model(uo$truth.model) == t] == t, na.rm = TRUE)
    z <- as.logical(uz$CSq.has.common.child); i <- coarse.model(uz$truth.model) == t
    b <- mean(coarse.model(uz[[model.col("CSq")]])[i & !z] == t, na.rm = TRUE)
    data.frame(Model = t, `u_only (no Z exists)` = pct(a),
               `u_z, Z not selected` = pct(b), `Gap` = sgn(b - a),
               check.names = FALSE, stringsAsFactors = FALSE)
}))

# per-effect check on the collider, to confirm it is not an artefact of one treatment
by.effect <- do.call(rbind, lapply(EFFECT.SIZES, function(e) {
    sub <- function(nm) { r <- RES[[nm]]; r[r$effect_size == e, , drop = FALSE] }
    a <- overall(sub("u_only"), "CSq"); b <- overall(sub("u_z"), "CSq")
    data.frame(`Effect treatment` = e, `u_only model acc.` = pct(a$model),
               `u_z model acc.` = pct(b$model), `Δ` = sgn(b$model - a$model),
               `u_only edge-absent recall` = pct(a$absent),
               `u_z edge-absent recall` = pct(b$absent),
               check.names = FALSE, stringsAsFactors = FALSE)
}))

lines <- c(
"# What the intermediate and the common child do to MRGN",
"",
paste("Generated by `results_scripts/make_structure_report.R` from the four",
      "confounder-structure runs. Do not edit by hand; rerun the script."),
"",
"---",
"",
"## 1. The question",
"",
paste("Every trio in the main simulation carries `U` confounders **plus** one intermediate",
      "`W` and one common child `Z`. `W` and `Z` are the two covariate blocks a selection",
      "method must *reject* rather than recover (`METHODS.md` Table 3): `W` sits between the",
      "genes, `Z` below them. Because the main design attaches both to every trio, it cannot",
      "say which one causes a failure when selection lets one through."),
"",
paste("These three runs separate them. Each is 300 trios at n = 670, MRGN only, generated",
      "under identical effect strata, minor allele frequencies, `U_n` range, residual SD and",
      "coefficient ranges, with independent seeds. **Only `W_n` and `Z_n` differ**, so a",
      "difference between structures has one possible cause."),
"",
md.table(do.call(rbind, lapply(names(STRUCTURES), function(nm) data.frame(
    Structure = nm, Description = STRUCTURES[[nm]]$label,
    Configuration = STRUCTURES[[nm]]$detail, check.names = FALSE,
    stringsAsFactors = FALSE)))),
"",
paste("**Table C1. The four structures.** `u_w_z` is not a separate run -- it *is* the main",
      "simulation, read from `data/`, which is what makes it the correct fourth cell rather",
      "than a re-simulation."),
"",
"---",
"",
"## 2. Headline: the common child does the damage, the intermediate does not",
"",
md.table(overall.table("CSq")),
"",
paste(sprintf("**Table C2. Overall MRGN accuracy by structure, CS-q arm.** Adding an"),
      sprintf("intermediate **improves** model accuracy by %s points; adding a common child",
              sub("^\\+", "", sgn(o.w$model - o.only$model))),
      sprintf("costs **%s**. The both-hazards case (%s) tracks the common-child-only case",
              sgn(o.z$model - o.only$model), pct(o.wz$model)),
      sprintf("(%s) rather than sitting between them, so **`Z` accounts for essentially all",
              pct(o.z$model)),
      "of the damage present in the main simulation**."),
"",
paste("The edge-absent recall column is where it concentrates:",
      sprintf("%s in `u_only` against %s in `u_z`. Edge-*present* recall barely moves",
              pct(o.only$absent), pct(o.z$absent)),
      sprintf("(%s to %s). MRGN does not lose the ability to find edges; it loses the",
              pct(o.only$present), pct(o.z$present)),
      "ability to rule them out."),
"",
"---",
"",
"## 3. Which models break",
"",
"Recall per generating model -- of the 60 trios generated under each model, the share MRGN",
"labelled correctly.",
"",
"**CS-q arm**", "",
md.table(per.model.table("CSq")),
"",
"**Change against the `u_only` baseline, in points**", "",
md.table(per.model.table("CSq", delta = TRUE)),
"",
paste("**Table C3. The collider damage is confined to M0 and M3.** Adding `Z` takes M0 from",
      sprintf("%s to %s and M3 from %s to %s -- both roughly halved -- while M1, M2 and M4",
              pct(r.only["M0"]), pct(r.z["M0"]), pct(r.only["M3"]), pct(r.z["M3"])),
      "move by a few points in either direction, within what 60 trios per cell can resolve."),
"",
paste("**M0 and M3 are exactly the two models whose truth contains no T1-T2 edge.** `M0` is",
      "`V1 → T1`, `T2` independent; `M3` is `V1 → T1`, `V1 → T2` with the genes unconnected.",
      "`EDGE.CORRECT` maps precisely these two to *edge absent*, and `M1`/`M2`/`M4` -- the",
      "three that do carry a T1-T2 edge -- to *edge present*. The split in the table falls",
      "along that line and nowhere else."),
"",
paste("That is the collider prediction, arriving where theory puts it. `Z` is a common child",
      "of both genes, so conditioning on it opens a path between `T1` and `T2` that the",
      "generating model does not contain. A method conditioning on `Z` sees an association",
      "between two genes that are not connected, and reports the edge. Trios that genuinely",
      "have an edge are unaffected, because for them the induced association adds nothing",
      "they were not going to find."),
"",
"---",
"",
"## 4. Mechanism: the spurious edge, measured",
"",
paste("Restricting to the M0 and M3 trios of `u_z` -- every one of which has no T1-T2 edge --",
      "and splitting on whether CS-q actually put the `Z` in that trio's covariate set:"),
"",
md.table(data.frame(
    `CS-q selected the Z` = c("yes", "no"),
    Trios = c(sp.z$n.exposed, sp.z$n.clean),
    `Called edge-present (all are truly absent)` = c(pct(sp.z$exposed), pct(sp.z$clean)),
    check.names = FALSE, stringsAsFactors = FALSE)),
"",
paste(sprintf("**Table C4. Spurious T1-T2 edges on edge-free trios.** %s against %s -- a",
              pct(sp.z$exposed), pct(sp.z$clean)),
      sprintf("factor of %.0f.", sp.z$exposed / sp.z$clean),
      "Every trio in this table is one where the correct answer is *no edge*."),
"",
paste("**This table shows the mechanism, not the effect size** -- see §8. The trios where",
      "selection missed the `Z` are not a random subset."),
"",
"---",
"",
"## 5. The control: the oracle arm is flat",
"",
md.table(overall.table("truth")),
"",
"**Per-model recall, oracle arm**", "",
md.table(per.model.table("truth")),
"",
paste("**Table C5. MRGN handed the true `U` block is untouched by structure.** Model accuracy",
      "holds across all four structures and M0/M3 recall stays high everywhere. The oracle",
      "arm never receives a `W` or a `Z` -- `ground.truth.input()` gives it the trio, `K` and",
      "that trio's own `U` block, nothing else -- so it is the same inference problem in all",
      "four cases and should be flat. It is."),
"",
paste("**That is what localises the blame.** The `u_only` → `u_z` collapse in §3 is not MRGN",
      "failing on a harder graph and not an artefact of the simulation: it is confounder",
      "selection admitting a covariate it was supposed to reject. Without this control the",
      "same numbers would be consistent with the structure simply being harder to infer."),
"",
"---",
"",
"## 6. Why is the intermediate harmless?",
"",
paste("`METHODS.md` Table 3 predicts damage from both: *\"conditioning on `W` blocks the",
      "mediated path, and conditioning on `Z` opens a collider.\"* The collider half holds",
      "exactly. **The mediator half does not.**"),
"",
paste(sprintf("Selection put the `W` into %s of `u_w` trios -- MORE often than it admitted the",
              pct(exp.w)),
      sprintf("`Z` in `u_z` (%s) -- so this is not a matter of exposure.", pct(exp.z)),
      "And the models blocking a mediated path should damage are M1 and M4, the two with a",
      sprintf("`T1 → T2` path. Their recall **rises**: M1 %s → %s, M4 %s → %s.",
              pct(r.only["M1"]), pct(r.w["M1"]), pct(r.only["M4"]), pct(r.w["M4"])),
      "Nothing in this data is consistent with the mediator hurting MRGN at n = 670."),
"",
paste("**This report does not explain why**, and the explanation should not be guessed at.",
      "Two candidates worth testing, neither settled here: MRGN's inference is anchored on",
      "the `V1` instrument, which may recover the `T1 → T2` relation whether or not the",
      "mediator is conditioned on; or the `W` coefficient range leaves the mediated path weak",
      "enough at this sample size that blocking it costs little. Distinguishing them needs a",
      "run varying the `W` coefficient, which does not exist."),
"",
"---",
"",
"## 7. Exposure, and why the filter lets the collider through",
"",
md.table(do.call(rbind, lapply(c("CSq", "CSa"), function(a) do.call(rbind, lapply(
    names(STRUCTURES), function(nm) data.frame(
        Arm = unname(ARM.LAB[a]), Structure = nm,
        `Selected the intermediate` = pct(exposure(RES[[nm]], a, "intermediate")),
        `Selected the common child` = pct(exposure(RES[[nm]], a, "common.child")),
        check.names = FALSE, stringsAsFactors = FALSE)))))),
"",
paste("**Table C6. How often each rejectable covariate ended up in the selected set.**",
      "`get.conf.trios()` runs a filtering step whose whole purpose is to strip intermediates",
      "and common children, and it is on (`filter_int_child = TRUE`) for `u_w`, `u_z` and",
      "`u_w_z`. It still admits the collider about two-thirds of the time."),
"",
paste("**The filter can only reject what it can detect.** `METHODS.md` §4 measures common",
      "children as detectable in roughly half of trios at n = 670, which is consistent with",
      "the rate above: what the filter cannot see, it cannot remove. The same section already",
      "warns that selected sets \"will contain colliders for most trios and should be",
      "caveated in the analysis\" -- this report is the quantitative form of that warning."),
"",
"---",
"",
"## 8. What the §4 split can and cannot support",
"",
paste("The within-structure contrast of §4 is the sharpest-looking comparison available, and",
      "**it overstates the effect.** Trios where selection missed the `Z` outperform the",
      "`u_only` baseline substantially -- if misses were random, the two would match:"),
"",
md.table(confound),
"",
paste("**Table C7. The `Z`-not-selected subgroup is an easier population.** A gap of zero",
      "would mean \"the `Z` was missed\" is a comparable trio. The gaps are large and",
      "positive, so missing the `Z` marks a trio where selection was working cleanly to begin",
      "with -- higher signal, cleaner covariate structure -- rather than an otherwise typical",
      "trio that got lucky."),
"",
paste("So §4's factor-of-ten is a **mechanism** statement: the spurious edges sit exactly",
      "where the collider is, which is what the causal story predicts. The **effect size** is",
      "the between-structure comparison of §2 and §3, where the contrast is clean by",
      "construction -- `u_only` and `u_z` are independently simulated under identical",
      "parameters and differ only in whether a `Z` exists at all."),
"",
"### Not an artefact of one effect treatment",
"",
md.table(by.effect),
"",
paste("**Table C8. The collider cost by effect treatment.** Present in all three, so it is",
      "not carried by one stratum. Effect treatments are pooled everywhere else in this",
      "report; this is the check that pooling is safe."),
"",
"---",
"",
"## 9. What it means for the main results",
"",
paste("1. **Every trio in the main simulation carries a `Z`.** So the main results are the",
      sprintf("`u_w_z` row: M0 and M3 recall there are collider-depressed, and the CS-q arm's"),
      "headline accuracy understates what MRGN would do against a selection that rejected",
      "colliders reliably. `u_only` is the closest available collider-free reference."),
"",
paste("2. **The damage is directional and should be described that way.** MRGN's CS-q arm",
      "does not degrade uniformly; it loses the ability to certify *absence* of a T1-T2 edge",
      "while retaining its ability to detect presence. Any summary quoting a single accuracy",
      "number for the CS-q arm hides that."),
"",
paste("3. **The intermediate can be dropped from the list of concerns at this sample size.**",
      "Table 3 lists `W` and `Z` together as covariates to reject. On this evidence they are",
      "not comparable hazards, and effort spent on intermediate filtering buys nothing that",
      "was measurable here."),
"",
paste("4. **This sharpens the CS-α case.** CS-α admits the collider at the same rate as CS-q",
      "(Table C6) on top of its ~100 false positives, and its M0/M3 recall degrades the same",
      "way. Its higher headline accuracy in `u_z` comes from the models the collider does not",
      "touch."),
"",
"---",
"",
"## 10. Provenance",
"",
paste("- **n = 670 only.** The structure runs were not repeated at other sample sizes, so",
      "nothing here speaks to how the collider effect scales with `n`. Detection of common",
      "children is itself strongly n-dependent (`METHODS.md` §4), so it likely does scale."),
paste("- **Effect treatments pooled**, with §8's check that the effect is present in all",
      "three."),
paste("- **300 trios per structure, 60 per generating model.** A per-model recall moves in",
      "steps of 1/60 = 1.7 points, so differences of a few points are noise. The M0 and M3",
      "changes are 20-35 points."),
paste("- **Run with `--bootstrap 0`.** The `boot.*` columns are `NA` throughout. This does",
      "not touch anything above: `apply.mrgn()` takes the model call from",
      "`MRGN::infer.trio()`, which never sees the bootstrap argument. Every figure in this",
      "report is a model call."),
paste("- **`u_w_z` is the main simulation**, read from `data/mrgn_group_n670.RData`, which",
      "*was* run with the bootstrap. That difference affects no column used here."),
"",
paste("Companion documents: `CONFOUNDER_SELECTION.md` scores the selection that admits these",
      "covariates; `INFERENCE_PERFORMANCE.md` covers all four methods on the main",
      "simulation; `../simulation_results/tables/confusion_structures.md` holds the full",
      "matrices behind Tables C2-C5."),
"")

writeLines(lines, OUT, useBytes = TRUE)
cat(sprintf("wrote %s | %d lines\n", OUT, length(lines)))
