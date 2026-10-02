# GTEx Whole Blood: rerun MRGN and GMAC inference WITHOUT the known covariates.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/gtex_rerun/rerun_gtex_no_known_covs.R [--cores=4]
#
# The GTEx expression was PEER-normalized with the GTEx-reported covariates already regressed
# out, so adding pcr, platform and sex again as known confounders adjusts for them twice.
# This script repeats the INFERENCE step only -- confounder selection is not rerun:
#
#   MRGN  CS-q   GTEx/data/data.with.PCs.WholeBlood.RData, columns 4:6 dropped
#   MRGN  CS-a   GTEx/trios_data_prlnc_liberal_conf_sel/data_with_confs/trios.with.confs.WholeBlood.RData,
#                columns 4:6 dropped
#                Both with infer.trio(use.perm = FALSE), the settings that produced the
#                MRGN.Inferred.Model.CSq / .Csalpha columns (checked by a gate per arm).
#   GMAC         GTEx/results/GMAC/GMAC_input_list_for_WholeBlood.RData with an EMPTY known
#                confounder matrix; settings as GTEx/results/GMAC/gtex.analysis.gmac.R.
#                GMAC's own PC selection never reads the known confounders, so the selected
#                PCs are unchanged (asserted below) and only the mediation tests differ.
#
# MRPC is not rerun here: at n = 670 with 20-50 confounders per trio most MRPC fits did not
# finish within 10 minutes in a test, so the MRPC-ADDIS column is left as it was.
#
# Memory: the two trio lists are ~0.5 and ~0.9 GB on disk and several GB in memory, so the
# arms are loaded, run and released one at a time. Every stage checkpoints to
# bioinfo_revision/gtex_rerun/updated_results/ and is skipped if its checkpoint exists, so a
# crashed run resumes where it stopped.
#
# Writes the new columns (suffix .noKC) to the end of
# bioinfo_revision/gtex_rerun/TableS3_GTEx_all_trios_master_updatedOct1.csv, plus a copy and
# the raw outputs in updated_results/.

suppressMessages({
  library(MRGN)
  library(parallel)
})

# the adapted GMAC code (gmac())
source("adapted_GMAC_func/GMAC_moded/R/GMAC.R")

OUT.DIR   <- "bioinfo_revision/gtex_rerun/updated_results"
TABLE     <- "bioinfo_revision/gtex_rerun/TableS3_GTEx_all_trios_master_updatedOct1.csv"
KC.NAMES  <- c("pcr", "platform", "sex")
args      <- commandArgs(trailingOnly = TRUE)
arg.val   <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[1]) else default
}
N.CORES   <- as.integer(arg.val("--cores", 4))
dir.create(OUT.DIR, recursive = TRUE, showWarnings = FALSE)

stamp <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
ckpt  <- function(name) file.path(OUT.DIR, name)

master <- read.csv(TABLE, check.names = FALSE)
stopifnot(nrow(master) == 3248)

cl <- makeCluster(N.CORES)
root <- normalizePath(getwd(), winslash = "/")
clusterExport(cl, "root", envir = environment())
invisible(clusterEvalQ(cl, {
  setwd(root)
  suppressMessages(library(MRGN))
  source("adapted_GMAC_func/GMAC_moded/R/GMAC.R")
}))
clusterSetRNGStream(cl, 2026)

# ---------------------------------------------------------------------------------------
# stages 1-2 -- MRGN, one arm at a time. For each arm:
#   gate:  with the covariates KEPT, infer.trio(use.perm = FALSE) must reproduce the stored
#          calls on 100 trios, or the old-vs-new comparison means nothing
#   drop:  columns 4:6 must be exactly pcr, platform, sex; remove them
#   run:   18 x 3248 (MRGN 1.3.9 adds coef11..coef22 to the old 14-row reg.res layout:
#          rows 1:6 indicators, 7:12 p-values, 13 minor allele freq, 18 Inferred.Model)
# ---------------------------------------------------------------------------------------

run.mrgn.arm <- function(data.file, old.calls, n.pcs, file, label) {
  if (file.exists(ckpt(file))) {
    stamp("stage 2: MRGN ", label, " loaded from checkpoint")
    return(loadRData(ckpt(file)))
  }
  trios <- loadRData(data.file)
  stopifnot(length(trios) == 3248)
  stamp("stage 1: ", label, " loaded")

  set.seed(2026)
  gate.idx <- sort(sample(length(trios), 100))
  gate <- unlist(parLapply(cl, trios[gate.idx],
                           function(x) infer.trio(x, use.perm = FALSE)$Inferred.Model))
  stamp("stage 1 gate: ", label, " reproduced ", sum(gate == old.calls[gate.idx]), "/100")
  stopifnot(all(gate == old.calls[gate.idx]))

  trios <- lapply(trios, function(x) {
    stopifnot(identical(colnames(x)[4:6], KC.NAMES))
    x[, -(4:6), drop = FALSE]
  })
  invisible(gc())
  # the trios line up with the master table, and only the three covariates went
  stopifnot(all(vapply(trios, function(x) colnames(x)[1], "") == master$SNP))
  stopifnot(all(vapply(trios, ncol, 1L) - 3 == n.pcs))
  stopifnot(!any(KC.NAMES %in% unlist(lapply(trios, colnames))))

  res <- simplify2array(parLapply(cl, trios, infer.trio, use.perm = FALSE), higher = FALSE)
  stopifnot(ncol(res) == 3248)
  save(res, file = ckpt(file))
  rm(trios); invisible(gc())
  stamp("stage 2: MRGN ", label, " done")
  res
}

mrgn.q <- run.mrgn.arm("GTEx/data/data.with.PCs.WholeBlood.RData",
                       master$MRGN.Inferred.Model.CSq, master$MRGN.number.of.PCs,
                       "mrgn_noKC_CSq.RData", "CS-q")
mrgn.a <- run.mrgn.arm("GTEx/trios_data_prlnc_liberal_conf_sel/data_with_confs/trios.with.confs.WholeBlood.RData",
                       master$MRGN.Inferred.Model.Csalpha, master$MRGN.libconf.alpha05.number.of.PCs,
                       "mrgn_noKC_CSalpha.RData", "CS-alpha")

# ---------------------------------------------------------------------------------------
# stage 4 -- GMAC, cis gene as mediator then trans gene as mediator, no known confounders
# ---------------------------------------------------------------------------------------

input.list <- loadRData("GTEx/results/GMAC/GMAC_input_list_for_WholeBlood.RData")
n.samp <- ncol(input.list$exp.dat)
no.kc  <- matrix(numeric(0), nrow = 0, ncol = n.samp)   # gmac() transposes it to n x 0

run.gmac <- function(trios.idx, file) {
  if (file.exists(ckpt(file))) return(loadRData(ckpt(file)))
  out <- gmac(cl = cl, known.conf = no.kc, cov.pool = input.list$cov.pool,
              exp.dat = input.list$exp.dat, snp.dat.cis = input.list$snp.dat.cis,
              trios.idx = trios.idx, nperm = 1000, nominal.p = TRUE,
              fdr = 0.05, fdr_filter = 0.1)
  save(out, file = ckpt(file))
  out
}
gmac.cis   <- run.gmac(input.list$trios.idx, "gmac_noKC_cis.RData")
stamp("stage 4: GMAC cis done")
gmac.trans <- run.gmac(input.list$trios.idx[, c(1, 3, 2)], "gmac_noKC_trans.RData")
stamp("stage 4: GMAC trans done")

# GMAC's PC selection never reads the known confounders: the selected PCs must be the ones
# the original run used
old.cis <- loadRData("GTEx/results/GMAC/GMAC_cis_inference_for_tissue_WholeBlood.RData")
stopifnot(identical(unname(rowSums(gmac.cis$sel.conf.ind)), unname(rowSums(old.cis$cov.indicator.list))))

# mediation call per trio -- verbatim from GTEx/scripts/make_mrgn_triotables.R
replace.tf <- function(tf.cols = NULL) {
  type.of.med <- NULL
  for (i in 1:dim(tf.cols)[1]) {
    if (isTRUE(all.equal(tf.cols[i, ], c(TRUE, TRUE), check.attributes = FALSE))) {
      type.of.med[i] <- "Undirected"
    } else if (isTRUE(all.equal(tf.cols[i, ], c(TRUE, FALSE), check.attributes = FALSE))) {
      type.of.med[i] <- "Cis Mediated"
    } else if (isTRUE(all.equal(tf.cols[i, ], c(FALSE, TRUE), check.attributes = FALSE))) {
      type.of.med[i] <- "Trans Mediated"
    } else if (isTRUE(all.equal(tf.cols[i, ], c(FALSE, FALSE), check.attributes = FALSE))) {
      type.of.med[i] <- "No Mediation"
    }
  }
  return(type.of.med)
}
qsig <- function(p) qvalue::qvalue(p, fdr.level = .1,
                                   lambda = seq(0.05, max(p, na.rm = T), 0.05))$significant

p.cis   <- gmac.cis$pvals[, "Known_sel_pool"]
p.trans <- gmac.trans$pvals[, "Known_sel_pool"]
gmac.tab <- data.frame(input.list$trio.ref,
                       cis_pval_Known_sel_pool = p.cis,
                       trans_pval_Known_sel_pool = p.trans,
                       cis_effect_change_Known_sel_pool = gmac.cis$beta.change[, "Known_sel_pool"],
                       trans_effect_change_Known_sel_pool = gmac.trans$beta.change[, "Known_sel_pool"],
                       type.of.med.qval = replace.tf(cbind(qsig(p.cis), qsig(p.trans))),
                       type.of.med.at.05.cutoff = replace.tf(cbind(p.cis < 0.05, p.trans < 0.05)),
                       type.of.med.at.01.cutoff = replace.tf(cbind(p.cis < 0.01, p.trans < 0.01)),
                       GMAC.number.of.PCs = rowSums(gmac.cis$sel.conf.ind))

# GMAC row order is the input list's trio.ref; map it onto the master table by trio IDs.
# make.unique() put ".1"-style suffixes on repeated gene columns, so strip them to compare.
strip <- function(x) sub("[.][0-9]+$", "", x)
key.gmac   <- paste(strip(gmac.tab$snp), strip(gmac.tab$cis), strip(gmac.tab$trans))
key.master <- paste(master$SNP, master$Cis.Gene.ID, master$Trans.Gene.ID)
stopifnot(!anyDuplicated(key.gmac), !anyDuplicated(key.master))
ord <- match(key.master, key.gmac)
stopifnot(!anyNA(ord))
gmac.tab <- gmac.tab[ord, ]

# the same mapping applied to the ORIGINAL run must reproduce the stored GMAC calls
old.trans <- loadRData("GTEx/results/GMAC/GMAC_trans_inference_for_tissue_WholeBlood.RData")
old.call <- replace.tf(cbind(old.cis$output.table[ord, 2] < 0.01, old.trans$output.table[ord, 2] < 0.01))
stopifnot(all(old.call == master$GMAC.Inferred.Mediation.Type.alpha01))
stopifnot(all(gmac.tab$GMAC.number.of.PCs == master$GMAC.number.of.PCs))
save(gmac.tab, file = ckpt("gmac_noKC_table.RData"))
stamp("stage 4: GMAC post-processed and aligned to the master table")

# ---------------------------------------------------------------------------------------
# stage 5 -- write the new columns and the old-vs-new summary
# ---------------------------------------------------------------------------------------

master$MRGN.Inferred.Model.CSq.noKC              <- unlist(mrgn.q["Inferred.Model", ])
master$MRGN.Inferred.Model.Csalpha.noKC          <- unlist(mrgn.a["Inferred.Model", ])
master$GMAC.Inferred.Mediation.Type.Qval.noKC    <- gmac.tab$type.of.med.qval
master$GMAC.Inferred.Mediation.Type.alpha01.noKC <- gmac.tab$type.of.med.at.01.cutoff
master$GMAC.Inferred.Mediation.Type.alpha05.noKC <- gmac.tab$type.of.med.at.05.cutoff
master$GMAC.cis.pval.noKC                        <- gmac.tab$cis_pval_Known_sel_pool
master$GMAC.trans.pval.noKC                      <- gmac.tab$trans_pval_Known_sel_pool

write.csv(master, TABLE, row.names = FALSE)
write.csv(master, file.path(OUT.DIR, basename(TABLE)), row.names = FALSE)
stamp("stage 5: wrote ", TABLE, " (", ncol(master), " columns)")

pairs <- list(
  "MRGN CS-q"            = c("MRGN.Inferred.Model.CSq", "MRGN.Inferred.Model.CSq.noKC"),
  "MRGN CS-alpha"        = c("MRGN.Inferred.Model.Csalpha", "MRGN.Inferred.Model.Csalpha.noKC"),
  "GMAC qval"            = c("GMAC.Inferred.Mediation.Type.Qval", "GMAC.Inferred.Mediation.Type.Qval.noKC"),
  "GMAC alpha01"         = c("GMAC.Inferred.Mediation.Type.alpha01", "GMAC.Inferred.Mediation.Type.alpha01.noKC"),
  "GMAC alpha05"         = c("GMAC.Inferred.Mediation.Type.alpha05", "GMAC.Inferred.Mediation.Type.alpha05.noKC"))
summ <- do.call(rbind, lapply(names(pairs), function(m) {
  old <- master[[pairs[[m]][1]]]; new <- master[[pairs[[m]][2]]]
  lv <- sort(unique(c(old, new)), na.last = TRUE)
  data.frame(method = m, call = ifelse(is.na(lv), "NA", lv),
             with.known.covs = vapply(lv, function(v) sum(old %in% v), 1L),
             without.known.covs = vapply(lv, function(v) sum(new %in% v), 1L),
             row.names = NULL)
}))
agree <- data.frame(method = names(pairs),
                    trio.agreement = vapply(pairs, function(p) mean(master[[p[1]]] == master[[p[2]]], na.rm = TRUE), 1))
write.csv(summ, file.path(OUT.DIR, "summary_old_vs_new.csv"), row.names = FALSE)
write.csv(agree, file.path(OUT.DIR, "agreement_old_vs_new.csv"), row.names = FALSE)
print(summ); print(agree)
stopCluster(cl)
stamp("done.")
