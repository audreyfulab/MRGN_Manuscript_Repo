# Apply GMAC's mediation test to the CS-i confounder set.
#
# NOTE: Set your working directory to the repository root before running this script.
#
#   Rscript bioinfo_revision/simulation_results/apply_gmac_csi.R --cores 2
#
# Writes gmac_csi_group_n<size>.RData per group, then inference_gmac_csi.RData/.csv.
# merge_csi.R then joins the gmac.CSi.* columns into inference_gmac.RData.
#
# ---------------------------------------------------------------------------------------
# THIS ARM IS A CHECK, NOT A NEW RESULT
# ---------------------------------------------------------------------------------------
#
# GMAC selects its own confounders, so unlike MRGN/MRPC/MR-GGI it does not need CS-i to run.
# The reason to run it anyway is that CS-i is claimed to BE GMAC's selection rule:
# get.conf.trios(adjust_by = "individual") corrects q-values within each covariate across
# trios, which is the family GMAC's conf.fdr() applies (adapted_GMAC_func/gmac_get_conf.R:98).
# The counts already agree closely -- 2.6 selected per trio at n = 50 against GMAC's 2.6,
# 9.6 at n = 300 against GMAC's 9.6 -- and this arm tests the identification where it
# actually matters: if the two selections agree, GMAC's mediation call on CS-i should track
# its call on its own set. A large divergence would mean the count agreement is coincidence
# and the sets differ, which is what run.gmac.csi.group()'s comment records as the open
# residual (GMAC restricts a covariate's family to trios passing its child filter; MRGN
# drops filtered covariates from the design instead).
#
# NO BATCH gmac() CALL. run.gmac.group() spends its time on the batch call that performs
# GMAC's own selection across the group, and a fixed confounder set makes it irrelevant.
# This is one apply.gmac() per trio -- gmac.nperm permutations x 2 (cis, then trans) --
# spread over the cluster. At ~0.79 s per trio that is roughly 20 min of core time for
# 1,500 trios.
#
# Needs a CS.i block in the selection cache: run backfill_csi.R first.

library(MRGN)
library(parallel)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

cat("=== GMAC / CS-i ===  cores:", n.cores, "| nperm:", gmac.nperm,
    "| sizes:", if (is.null(sample.sizes)) "all" else paste(sample.sizes, collapse = ","), "\n")

cl <- parallel::makeCluster(n.cores)
on.exit(parallel::stopCluster(cl), add = TRUE)

# The workers call apply.gmac(), which is defined in inference_utils.R and reaches
# gmacOneTrio() from adapted_GMAC_func. Both are sourced on each worker rather than
# exported: apply.gmac()'s closure environment is the master's global env, which is not
# serialised with it, so exporting the function alone would leave gmacOneTrio() unfound.
#
# setwd(root) first, because the sources are relative paths and a PSOCK worker starts in
# its own working directory.
root <- normalizePath(getwd(), winslash = "/")
parallel::clusterExport(cl, "root", envir = environment())
invisible(parallel::clusterEvalQ(cl, {
    setwd(root)
    library(MRGN)
    source("bioinfo_revision/simulation/simulation_utils.R")
    source("bioinfo_revision/simulation_results/inference_config.R")
    source("bioinfo_revision/simulation_results/inference_utils.R")
    NULL
}))
parallel::clusterSetRNGStream(cl, 234)

run.method.groups(
    method = "gmac_csi",
    needs.selection = TRUE,       # unlike apply_gmac.R: the CS-i set comes from the cache
    cl = cl,
    runner = function(datasets, sel, cov.pool, known.conf, cl) {
        run.gmac.csi.group(datasets, sel, colnames(cov.pool), cl = cl,
                           nperm = gmac.nperm, known.conf = known.conf, verbose = TRUE)
    })

cat("\ncombining GMAC / CS-i:\n")
invisible(combine.method("gmac_csi", sample.sizes = sample.sizes))
