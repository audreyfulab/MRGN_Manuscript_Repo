# Changelog: legacy MRGN GTEx results vs the rerun without the known covariates.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/gtex_rerun/make_mrgn_changelog.R
#
# Run bioinfo_revision/gtex_rerun/rerun_gtex_no_known_covs.R first.
#
# Legacy  = MRGN with pcr, platform, sex as known covariates, no permutation test:
#             CS-q     Permutation_test_analysis/results/no.perm.all.trios.WB.RData
#             CS-alpha Permutation_test_analysis/results/no.perm.all.trios.WB.liberal.confs.alpha05.RData
#           (14 x 3248; row 14 is identical to MRGN.Inferred.Model.CSq / .Csalpha, asserted below)
# Rerun   = the same trios and confounders with the three covariates dropped:
#             updated_results/mrgn_noKC_CSq.RData, mrgn_noKC_CSalpha.RData (18 x 3248)
# The rerun reproduced the legacy calls on 100/100 trios per arm with the covariates kept,
# so every difference below is due to dropping the covariates.
#
# Writes
#   bioinfo_revision/gtex_rerun/MRGN_CHANGELOG.md                     narrative changelog
#   updated_results/mrgn_changelog_trios.csv                          every trio whose model or
#                                                                     any edge indicator changed
#   updated_results/mrgn_changelog_all_trios.csv                      old vs new for all trios
#   updated_results/mrgn_transition_<arm>.csv                         old x new model counts
#   updated_figures/fig_mrgn_legacy_vs_rerun.{pdf,png}

suppressMessages({
  library(MRGN)
  library(ggplot2)
  library(ggpubr)
})
source('Manuscript/scripts/MRGN_write_up_helper_functions.R')

ROOT    <- "bioinfo_revision/gtex_rerun"
RES.DIR <- file.path(ROOT, "updated_results")
FIG.DIR <- file.path(ROOT, "updated_figures")
dir.create(FIG.DIR, recursive = TRUE, showWarnings = FALSE)

master <- read.csv(file.path(ROOT, "TableS3_GTEx_all_trios_master_updatedOct1.csv"), check.names = FALSE)
stopifnot(nrow(master) == 3248)

ARMS <- list(
  CSq     = list(label = "CS-q",
                 old.file = "Permutation_test_analysis/results/no.perm.all.trios.WB.RData",
                 new.file = file.path(RES.DIR, "mrgn_noKC_CSq.RData"),
                 old.col = "MRGN.Inferred.Model.CSq", new.col = "MRGN.Inferred.Model.CSq.noKC",
                 pcs.col = "MRGN.number.of.PCs"),
  CSalpha = list(label = "CS-α",
                 old.file = "Permutation_test_analysis/results/no.perm.all.trios.WB.liberal.confs.alpha05.RData",
                 new.file = file.path(RES.DIR, "mrgn_noKC_CSalpha.RData"),
                 old.col = "MRGN.Inferred.Model.Csalpha", new.col = "MRGN.Inferred.Model.Csalpha.noKC",
                 pcs.col = "MRGN.libconf.alpha05.number.of.PCs"))

IND    <- c("b11", "b12", "b21", "b22", "V1:T1", "V1:T2")
PV     <- c("pb11", "pb12", "pb21", "pb22", "pV1:T1", "pV1:T2")
MODELS <- c("M0.1", "M0.2", "M1.1", "M1.2", "M2.1", "M2.2", "M3", "M4", "Other")

num.rows <- function(res, rows) {
  m <- sapply(rows, function(r) as.numeric(unlist(res[r, ])))
  colnames(m) <- rows
  m
}
# some genes have no BioMart name in the master table; fall back to the Ensembl ID
gene.label <- function(name, id) ifelse(is.na(name) | name == "", id, name)
t1t2.edge <- function(model, ind) {
  vapply(seq_along(model), function(i)
    ind.med.edge(get.adj.from.class(model[i], reg.vec = ind[i, ])), numeric(1))
}

# ---------------------------------------------------------------------------------------
# per-arm diff
# ---------------------------------------------------------------------------------------

diffs <- lapply(names(ARMS), function(a) {
  arm <- ARMS[[a]]
  old <- loadRData(arm$old.file)
  new <- loadRData(arm$new.file)
  stopifnot(ncol(old) == 3248, ncol(new) == 3248)
  old.mod <- as.character(unlist(old["Inferred.Model", ]))
  new.mod <- as.character(unlist(new["Inferred.Model", ]))
  stopifnot(all(old.mod == master[[arm$old.col]]), all(new.mod == master[[arm$new.col]]))

  old.ind <- num.rows(old, IND); new.ind <- num.rows(new, IND)
  old.p   <- num.rows(old, PV);  new.p   <- num.rows(new, PV)
  flipped <- old.ind != new.ind

  d <- data.frame(
    arm = a,
    trio.index = seq_len(3248),
    SNP = master$SNP,
    cis.gene = gene.label(master$Cis.Gene.Name, master$Cis.Gene.ID),
    trans.gene = gene.label(master$Trans.Gene.Name, master$Trans.Gene.ID),
    n.PCs = master[[arm$pcs.col]],
    old.model = old.mod,
    new.model = new.mod,
    old.class = convert.cats(old.mod),
    new.class = convert.cats(new.mod),
    model.changed = old.mod != new.mod,
    class.changed = convert.cats(old.mod) != convert.cats(new.mod),
    indicators.changed = apply(flipped, 1, function(f) paste(IND[f], collapse = ";")),
    old.T1T2.edge = t1t2.edge(old.mod, old.ind),
    new.T1T2.edge = t1t2.edge(new.mod, new.ind),
    check.names = FALSE, stringsAsFactors = FALSE)
  for (k in seq_along(IND)) {
    d[[paste0("old.", IND[k])]] <- old.ind[, k]
    d[[paste0("new.", IND[k])]] <- new.ind[, k]
  }
  for (k in seq_along(PV)) {
    d[[paste0("old.", PV[k])]] <- old.p[, k]
    d[[paste0("new.", PV[k])]] <- new.p[, k]
  }
  d$any.change <- d$model.changed | d$indicators.changed != ""
  d
})
names(diffs) <- names(ARMS)
all.d <- do.call(rbind, diffs)

write.csv(all.d, file.path(RES.DIR, "mrgn_changelog_all_trios.csv"), row.names = FALSE)
write.csv(subset(all.d, any.change), file.path(RES.DIR, "mrgn_changelog_trios.csv"), row.names = FALSE)

trans <- lapply(diffs, function(d) table(legacy = factor(d$old.model, MODELS),
                                         rerun  = factor(d$new.model, MODELS)))
for (a in names(trans)) {
  m <- as.data.frame.matrix(trans[[a]])
  write.csv(cbind(legacy = rownames(m), m), file.path(RES.DIR, paste0("mrgn_transition_", a, ".csv")),
            row.names = FALSE)
}

# ---------------------------------------------------------------------------------------
# console summary
# ---------------------------------------------------------------------------------------

summ <- do.call(rbind, lapply(names(diffs), function(a) {
  d <- diffs[[a]]
  data.frame(arm = a,
             model.changed = sum(d$model.changed),
             class.changed = sum(d$class.changed),
             subtype.only  = sum(d$model.changed & !d$class.changed),
             indicator.only = sum(!d$model.changed & d$indicators.changed != ""),
             any.change    = sum(d$any.change),
             edge.gained   = sum(d$old.T1T2.edge == 0 & d$new.T1T2.edge == 1),
             edge.lost     = sum(d$old.T1T2.edge == 1 & d$new.T1T2.edge == 0))
}))
print(summ)

# ---------------------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------------------

INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3df"; DIAG <- "#f0efec"
SEQ <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
GAIN <- "#2a78d6"; LOSS <- "#e34948"

base.theme <- theme_minimal(base_size = 13) +
  theme(text = element_text(colour = INK),
        axis.text = element_text(colour = INK2),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = INK2, size = 11),
        plot.background = element_rect(fill = "white", colour = NA))

heat <- function(a) {
  m <- as.data.frame(trans[[a]], responseName = "n")
  m$diag <- as.character(m$legacy) == as.character(m$rerun)
  max.off <- max(m$n[!m$diag], 1)
  m$fill <- ifelse(m$diag, NA, ifelse(m$n == 0, NA, m$n))
  m$lab.col <- ifelse(!m$diag & m$n > max.off * 0.55, "white", ifelse(m$diag, INK2, INK))
  m$lab <- ifelse(m$n == 0, "", format(m$n, big.mark = ","))
  d <- diffs[[a]]
  ggplot(m, aes(rerun, legacy)) +
    geom_tile(data = subset(m, diag), fill = DIAG, colour = "white", linewidth = 1) +
    geom_tile(data = subset(m, !diag & n > 0), aes(fill = fill), colour = "white", linewidth = 1) +
    geom_text(aes(label = lab, colour = I(lab.col)), size = 3.6) +
    scale_fill_gradientn(colours = SEQ, limits = c(1, max.off), name = "trios changed",
                         breaks = function(l) unique(round(pretty(l)))[unique(round(pretty(l))) >= 1]) +
    scale_y_discrete(limits = rev(MODELS)) +
    scale_x_discrete(position = "top") +
    coord_fixed() +
    labs(title = paste0("MRGN, ", ARMS[[a]]$label),
         subtitle = sprintf("%d of 3,248 trios changed model (%.1f%%); gray diagonal = unchanged",
                            sum(d$model.changed), 100 * mean(d$model.changed)),
         x = "Rerun without pcr, platform, sex", y = "Legacy (with covariates)") +
    base.theme +
    theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"),
          legend.title = element_text(size = 10, colour = INK2),
          axis.title.x.top = element_text(margin = margin(b = 6)))
}

flips <- function(a) {
  d <- diffs[[a]]
  f <- do.call(rbind, lapply(IND, function(k) {
    o <- d[[paste0("old.", k)]]; n <- d[[paste0("new.", k)]]
    data.frame(indicator = k,
               direction = c("gained (0 → 1)", "lost (1 → 0)"),
               n = c(sum(o == 0 & n == 1), -sum(o == 1 & n == 0)))
  }))
  f$indicator <- factor(f$indicator, levels = rev(IND))
  f$direction <- factor(f$direction, levels = c("lost (1 → 0)", "gained (0 → 1)"))
  lim <- max(abs(f$n), 1) * 1.25
  ggplot(f, aes(n, indicator, fill = direction)) +
    geom_vline(xintercept = 0, colour = INK2, linewidth = 0.4) +
    geom_col(width = 0.6) +
    geom_text(data = subset(f, n != 0), aes(label = abs(n), hjust = ifelse(n > 0, -0.3, 1.3)),
              colour = INK, size = 3.6) +
    scale_fill_manual(values = c("lost (1 → 0)" = LOSS, "gained (0 → 1)" = GAIN), name = NULL,
                      drop = FALSE) +
    scale_x_continuous(limits = c(-lim, lim), labels = function(x) abs(x)) +
    labs(title = paste0("Edge indicators flipped, ", ARMS[[a]]$label),
         subtitle = sprintf("%d trios with at least one flipped indicator", sum(d$indicators.changed != "")),
         x = "trios", y = NULL) +
    base.theme +
    theme(legend.position = "bottom",
          panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3))
}

grDevices::pdf(NULL)
fig <- ggarrange(heat("CSq"), heat("CSalpha"), flips("CSq"), flips("CSalpha"),
                 labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, heights = c(1.35, 1),
                 font.label = list(size = 16, face = "bold", color = INK))
invisible(grDevices::dev.off())

cairo_pdf(file.path(FIG.DIR, "fig_mrgn_legacy_vs_rerun.pdf"), width = 14, height = 13)
plot(fig); invisible(dev.off())
png(file.path(FIG.DIR, "fig_mrgn_legacy_vs_rerun.png"), width = 14, height = 13, units = "in", res = 300)
plot(fig); invisible(dev.off())

# ---------------------------------------------------------------------------------------
# markdown changelog
# ---------------------------------------------------------------------------------------

md.table <- function(df) {
  df[] <- lapply(df, as.character)
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep(" --- ", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
count.tab <- function(a) {
  d <- diffs[[a]]
  lv <- c("M0", "M1", "M2", "M3", "M4", "Other")
  o <- table(factor(d$old.class, lv)); n <- table(factor(d$new.class, lv))
  data.frame(class = lv, legacy = as.vector(o), rerun = as.vector(n),
             net = sprintf("%+d", as.vector(n) - as.vector(o)))
}
trans.md <- function(a) {
  m <- as.data.frame.matrix(trans[[a]])
  keep.r <- rowSums(m) > 0; keep.c <- colSums(m) > 0
  m <- m[keep.r, keep.c, drop = FALSE]
  for (i in rownames(m)) for (j in colnames(m)) if (i == j) m[i, j] <- paste0("*", m[i, j], "*")
  md.table(cbind(`legacy \\ rerun` = rownames(m), m))
}
trio.md <- function(a) {
  d <- subset(diffs[[a]], model.changed)
  d <- d[order(d$old.model, d$new.model, d$trio.index), ]
  md.table(data.frame(trio = d$trio.index, SNP = d$SNP,
                      `cis gene` = d$cis.gene, `trans gene` = d$trans.gene, PCs = d$n.PCs,
                      legacy = d$old.model, rerun = d$new.model,
                      `indicators flipped` = d$indicators.changed,
                      `T1-T2 edge` = paste0(d$old.T1T2.edge, " → ", d$new.T1T2.edge),
                      check.names = FALSE))
}

s <- function(a, col) summ[summ$arm == a, col]
md <- c(
  "# MRGN GTEx changelog: legacy vs rerun without known covariates",
  "",
  "Generated by `make_mrgn_changelog.R`. Do not edit by hand; rerun the script instead.",
  "",
  "## What changed in the analysis",
  "",
  "| | legacy | rerun |",
  "| --- | --- | --- |",
  "| Known covariates in every trio | pcr, platform, sex | none (already removed by PEER normalization) |",
  "| Confounder PCs | CS-q / CS-α sets | identical sets, not reselected |",
  "| MRGN call | `infer.trio(use.perm = FALSE)` | identical |",
  "| Trios | 3,248 whole-blood trios | identical |",
  "| Source | `Permutation_test_analysis/results/no.perm.all.trios.WB*.RData` | `updated_results/mrgn_noKC_*.RData` |",
  "",
  "With the covariates kept, the rerun reproduced the legacy calls on 100/100 sampled trios",
  "in each arm, so every difference listed here comes from dropping the three covariates.",
  "",
  "## Summary",
  "",
  md.table(data.frame(
    ` ` = c("trios whose model changed", "... of which the model class changed",
            "... of which only the subtype changed (e.g. M0.1 → M0.2)",
            "trios with a flipped edge indicator but the same model",
            "trios with any change", "T1-T2 edge gained", "T1-T2 edge lost"),
    `CS-q` = unlist(summ[summ$arm == "CSq", -1]),
    `CS-α` = unlist(summ[summ$arm == "CSalpha", -1]), check.names = FALSE)),
  "",
  "Figure: `updated_figures/fig_mrgn_legacy_vs_rerun.png` (A, B transition matrices; C, D",
  "edge indicators that flipped).",
  "")
for (a in names(ARMS)) {
  md <- c(md,
    paste0("## ", ARMS[[a]]$label),
    "",
    "### Model class counts",
    "",
    md.table(count.tab(a)),
    "",
    "### Transitions (rows legacy, columns rerun; diagonal in italics = unchanged)",
    "",
    trans.md(a),
    "",
    sprintf("### Every trio whose model changed (%d)", s(a, "model.changed")),
    "",
    trio.md(a),
    "")
}
md <- c(md,
  "## Files",
  "",
  "- `updated_results/mrgn_changelog_trios.csv`: every trio (per arm) whose model or any",
  "  edge indicator changed, with old and new indicators and p-values.",
  "- `updated_results/mrgn_changelog_all_trios.csv`: the same columns for all 3,248 trios per arm.",
  "- `updated_results/mrgn_transition_CSq.csv`, `mrgn_transition_CSalpha.csv`: old × new model counts.",
  "")
writeLines(md, file.path(ROOT, "MRGN_CHANGELOG.md"), useBytes = TRUE)

cat("wrote", file.path(ROOT, "MRGN_CHANGELOG.md"), "\n")
cat("wrote", file.path(FIG.DIR, "fig_mrgn_legacy_vs_rerun.{pdf,png}"), "\n")
