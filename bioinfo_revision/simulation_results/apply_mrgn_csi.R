# Apply MRGN to the CS-i confounder set ONLY.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/apply_mrgn_csi.R --cores 2
#
# Writes mrgn_csi_group_n<size>.RData per group, then inference_mrgn_csi.RData/.csv.
# merge_csi.R then joins the mrgn.CSi.* columns into inference_mrgn.RData.
#
# WHY A SEPARATE PASS rather than adding CSi to apply_mrgn.R and re-running with
# --rerun-inference 1: that would recompute truth, CS-q and CS-alpha as well, and for MRGN
# specifically it would redo the bootstrap, which is n.bootstrap resamples x 3 arms x 1,500
# trios. The three existing arms are already on disk and unchanged by this work. Running
# CS-i alone and merging the columns costs one arm instead of four.
#
# The "_csi" in method = makes method.checkpoint() and combine.method() write to their own
# filenames, so nothing here can overwrite the existing mrgn_group_n*.RData.
#
# Needs a CS.i block in the selection cache: run backfill_csi.R first.

library(MRGN)
library(parallel)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

cat("=== MRGN / CS-i ===  cores:", n.cores,
    "| bootstrap:", if (mrgn.bootstrap) n.bootstrap else "OFF (boot.* columns will be NA)",
    "| sizes:", if (is.null(sample.sizes)) "all" else paste(sample.sizes, collapse = ","), "\n")

cl <- parallel::makeCluster(n.cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
invisible(parallel::clusterEvalQ(cl, library(MRGN)))
parallel::clusterSetRNGStream(cl, 234)

run.method.groups(
    method = "mrgn_csi",
    needs.selection = TRUE,
    cl = cl,
    runner = function(datasets, sel, cov.pool, known.conf, cl) {
        run.mrgn.group(datasets, sel, colnames(cov.pool), cl = cl,
                       bootstrap = mrgn.bootstrap, number_of_samples = n.bootstrap,
                       arms = "CSi", verbose = TRUE)
    })

cat("\ncombining MRGN / CS-i:\n")
invisible(combine.method("mrgn_csi", sample.sizes = sample.sizes))
