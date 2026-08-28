# Apply MR-GGI to the CS-i confounder set ONLY.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/apply_mrggi_csi.R --cores 2
#
# Writes mrggi_csi_group_n<size>.RData per group, then inference_mrggi_csi.RData/.csv.
# merge_csi.R then joins the mrggi.CSi.* columns into inference_mrggi.RData.
#
# WHY A SEPARATE PASS: re-running apply_mrggi.R with CSi added to mrggi.arms would recompute
# the CS-alpha arm too, and that arm alone is 9.4 h of the 10.4 h the full MR-GGI run costs
# -- it carries a median of 82-106 covariates, so MRggi computes ~5,800 gene pairs per trio.
# CS-i carries 2-24, i.e. between the `none` arm (1 pair) and CS-q, so this pass should cost
# on the order of the CS-q arm's ~12 min rather than the CS-alpha arm's 9.4 h.
#
# The T1-T2 raw p-value is ARM-INVARIANT by construction -- MRggi's estimator is pairwise,
# so .TSLS() sees only T1, T2 and the instrument regardless of what else rides along in y.
# The CS-i arm therefore cannot change mrggi.*.edge; it can only change edge.fdr, via the
# multiplicity correction across each gene's pairs. confusion_mrggi.R asserts that
# invariance, and adding a fifth arm gives it a fifth column to check it against.
#
# Needs a CS.i block in the selection cache: run backfill_csi.R first.

library(MRGN)
library(MRggi)
library(parallel)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

cat("=== MR-GGI / CS-i ===  alpha:", mrggi.alpha, "| cor.thr:", mrggi.cor.thr,
    "| min first-stage F:", mrggi.min.F, "| cores:", n.cores,
    "| sizes:", if (is.null(sample.sizes)) "all" else paste(sample.sizes, collapse = ","), "\n")

cl <- parallel::makeCluster(n.cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
    library(MRGN)
    library(MRggi)
}))

# Same export list and same reasoning as apply_mrggi.R: the workers evaluate
# mrggi.run.arms() on a self-contained per-trio payload and reach these as globals. They are
# exported by name rather than re-sourced, because inference_config.R reads commandArgs()
# and a PSOCK worker's command line is its own, not this process's.
parallel::clusterExport(cl, c(
    "mrggi.run.arms", "mrggi.one.trio", "mrggi.fields", "prefixed", "safely",
    "mrggi.arms", "mrggi.cor.thr", "mrggi.alpha", "mrggi.min.F", "mrggi.p.adjust"))

run.method.groups(
    method = "mrggi_csi",
    needs.selection = TRUE,
    cl = cl,
    runner = function(datasets, sel, cov.pool, known.conf, cl) {
        run.mrggi.group(datasets, sel, colnames(cov.pool), cl = cl,
                        arms = "CSi", verbose = TRUE)
    })

cat("\ncombining MR-GGI / CS-i:\n")
invisible(combine.method("mrggi_csi", sample.sizes = sample.sizes))
