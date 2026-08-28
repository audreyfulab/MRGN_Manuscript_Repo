# Apply MRPC to the simulated trios.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/apply_mrpc.R
#
# Which confounder sets are attempted is set by `mrpc.arms` in inference_config.R, and is
# currently the true confounders and CS-q: the CS-alpha set is too large for MRPC to fit,
# timing out on 79% of trios at n = 150 and heading to ~100% at larger n. See the note on
# mrpc.arms for the measurements. A disabled arm still gets its columns, with the reason in
# mrpc.<arm>.error, so every group's checkpoint keeps the same schema.
#
# THE TRUTH ARM is the same oracle MRGN gets -- trio + K + that trio's own U block. It is
# what lets the MRPC and MRGN columns be read against each other: the gap between the oracle
# and CS-q is the cost of confounder selection rather than of the method. It is also the arm
# most likely to time out, carrying a median of 25-29 confounders at n = 670 against CS-q's
# 16, where CS-q already times out 61% of the time. Smoke-test it with
# --sizes 670 --max-per-group 20 before committing to a full run.
#
# Writes mrpc_group_n<size>.RData per sample-size group, then inference_mrpc.RData/.csv.
#
# MRPC is single threaded -- MRPC() itself does not parallelise and apply.mrpc() wraps one
# call in R.utils::withTimeout() -- so this script never builds a cluster. That is also
# why run_all_inference.R gives it a single core and hands the rest to MRGN and GMAC.
# The timeout is the thing to watch: mrpc.timeout seconds per fit, two fits per trio.

library(MRGN)
library(MRPC)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

cat("=== MRPC ===  arms:", paste(mrpc.arms, collapse = ", "),
    "| timeout:", mrpc.timeout, "s per fit |",
    "sizes:", if (is.null(sample.sizes)) "all" else paste(sample.sizes, collapse = ","), "\n")

run.method.groups(
    method = "mrpc",
    needs.selection = TRUE,
    cl = NULL,
    runner = function(datasets, sel, cov.pool, known.conf, cl) {
        run.mrpc.group(datasets, sel, colnames(cov.pool),
                       timeout = mrpc.timeout, verbose = TRUE)
    })

cat("\ncombining MRPC:\n")
invisible(combine.method("mrpc", sample.sizes = sample.sizes))
