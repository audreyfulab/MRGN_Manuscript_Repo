# Apply MRPC to the simulated trios.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/apply_mrpc.R
#
# Two confounder settings per trio, the CS-q and CS-alpha selected sets. MRPC gets no
# ground-truth arm: the comparison of interest is against MRGN under the same selected
# confounders.
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

cat("=== MRPC ===  timeout:", mrpc.timeout, "s per fit |",
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
