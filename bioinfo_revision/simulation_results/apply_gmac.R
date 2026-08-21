# Apply GMAC to the simulated trios.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/apply_gmac.R
#
# GMAC selects its own confounders, so this script never reads the CS-q/CS-alpha cache --
# run.gmac.all() does the selection and the mediation test in one batch call per group,
# and gmac.* columns record which covariates it chose.
#
# Writes gmac_group_n<size>.RData per sample-size group, then inference_gmac.RData/.csv.
#
# Cost is gmac.nperm permutations x 2 batch calls (cis mediator, then trans) per group,
# parallelised over the cluster. gmac() dispatches child.p/conf.fdr/getp.func through
# parLapply and those helpers are sourced rather than packaged, so the workers have to
# source them for themselves.

library(MRGN)
library(parallel)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

cat("=== GMAC ===  cores:", n.cores, "| nperm:", gmac.nperm,
    "| sizes:", if (is.null(sample.sizes)) "all" else paste(sample.sizes, collapse = ","), "\n")

cl <- parallel::makeCluster(n.cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
root <- normalizePath(getwd(), winslash = "/")
parallel::clusterExport(cl, "root", envir = environment())
invisible(parallel::clusterEvalQ(cl, {
    setwd(root)
    library(MRGN)
    source("adapted_GMAC_func/gmac_one_trio.R")
    source("adapted_GMAC_func/GMAC_moded/R/GMAC.R")
    NULL
}))
parallel::clusterSetRNGStream(cl, 234)

run.method.groups(
    method = "gmac",
    needs.selection = FALSE,      # GMAC picks its own confounders
    cl = cl,
    runner = function(datasets, sel, cov.pool, known.conf, cl) {
        run.gmac.group(datasets, cov.pool, known.conf, cl = cl,
                       nperm = gmac.nperm, verbose = TRUE)
    })

cat("\ncombining GMAC:\n")
invisible(combine.method("gmac", sample.sizes = sample.sizes))
