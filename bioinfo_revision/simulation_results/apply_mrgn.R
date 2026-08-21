# Apply MRGN to the simulated trios.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/apply_mrgn.R
#
# Three confounder settings per trio:
#   the TRUE confounders (trio + K + that trio's own U block) -- the ceiling any
#     selection method is working toward
#   the CS-q selected confounders
#   the CS-alpha selected confounders
#
# Writes mrgn_group_n<size>.RData per sample-size group, then inference_mrgn.RData/.csv.
# Nothing here touches the MRPC or GMAC outputs, so all three can run concurrently --
# see run_all_inference.R.
#
# This is the expensive method: the bootstrap is n.bootstrap replicates x 3 settings x
# every trio, parallelised over the cluster.

library(MRGN)
library(parallel)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

cat("=== MRGN ===  cores:", n.cores, "| bootstrap:", n.bootstrap,
    "| sizes:", if (is.null(sample.sizes)) "all" else paste(sample.sizes, collapse = ","), "\n")

cl <- parallel::makeCluster(n.cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
invisible(parallel::clusterEvalQ(cl, library(MRGN)))
parallel::clusterSetRNGStream(cl, 234)

run.method.groups(
    method = "mrgn",
    needs.selection = TRUE,
    cl = cl,
    runner = function(datasets, sel, cov.pool, known.conf, cl) {
        run.mrgn.group(datasets, sel, colnames(cov.pool), cl = cl,
                       bootstrap = TRUE, number_of_samples = n.bootstrap, verbose = TRUE)
    })

cat("\ncombining MRGN:\n")
invisible(combine.method("mrgn", sample.sizes = sample.sizes))
