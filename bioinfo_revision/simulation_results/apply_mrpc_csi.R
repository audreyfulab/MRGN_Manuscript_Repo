# Apply MRPC to the CS-i confounder set ONLY.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/apply_mrpc_csi.R --cores 2
#
# Writes mrpc_csi_group_n<size>.RData per group, then inference_mrpc_csi.RData/.csv.
# merge_csi.R then joins the mrpc.CSi.* columns into inference_mrpc.RData.
#
# MRPC() does not parallelise -- apply.mrpc() wraps a single call in
# R.utils::withTimeout() -- so --cores is accepted and ignored here, exactly as in
# apply_mrpc.R. The 2 cores this is launched with buy nothing for MRPC; they are the
# per-method budget the other three runs are held to.
#
# COVERAGE: this follows the existing MRPC run and covers n = 50/150/300 only. The n = 670
# and n = 1000 groups have never been run under the 180 s cap (see inference_config.R), and
# adding a CS-i column for them here would produce a file the existing mrpc rows have no
# counterpart for. Pass --sizes to change that deliberately.
#
# CS-i sits between CS-q and the oracle in set size (median 2-24 covariates against CS-q's
# 0-23 and truth's 25-29), so its timeout risk is closer to CS-q's -- which never timed out
# at n <= 300 -- than to the truth arm's, which lost 103 of 900 fits at n = 300.
#
# Needs a CS.i block in the selection cache: run backfill_csi.R first.

library(MRGN)
library(MRPC)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

if (is.null(sample.sizes)) sample.sizes <- c(50, 150, 300)

cat("=== MRPC / CS-i ===  timeout:", mrpc.timeout, "s per fit |",
    "sizes:", paste(sample.sizes, collapse = ","), "\n")

run.method.groups(
    method = "mrpc_csi",
    needs.selection = TRUE,
    cl = NULL,
    runner = function(datasets, sel, cov.pool, known.conf, cl) {
        run.mrpc.group(datasets, sel, colnames(cov.pool), timeout = mrpc.timeout,
                       arms = "CSi", verbose = TRUE)
    })

cat("\ncombining MRPC / CS-i:\n")
invisible(combine.method("mrpc_csi", sample.sizes = sample.sizes))
