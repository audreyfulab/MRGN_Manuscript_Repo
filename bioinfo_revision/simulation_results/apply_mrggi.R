# Apply MR-GGI to the simulated trios.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/apply_mrggi.R
#
# MR-GGI is a Mendelian randomisation method: it uses V1 as a genetic instrument to estimate
# the causal effect between the two genes, rather than adjusting for confounders. See
# MRGGI_METHODS.md.
#
# ONE MRggi() CALL PER TRIO PER ARM. MRggi() takes a whole gene matrix, so the trio goes in
# as y = (T1, T2, covariates) instead of being called once per edge. The arms in mrggi.arms
# are `none` (bare trio, what this script did before), `truth`, `CSq` and `CSa`, differing
# only in which covariates ride along in y.
#
# THE ARMS ARE NOT A CONFOUNDER ADJUSTMENT and must not be read as one. MRggi's estimator is
# pairwise -- .TSLS() sees only the two genes and their instruments -- so the T1-T2 Wald
# ratio is identical in all four arms. What the covariates change is the multiplicity
# correction: MRggi adjusts each gene's p-values across that gene's pairs, so the T1-T2
# p-value is corrected for T1's pairs with every covariate as well. That is the whole
# difference, and inference_config.R's mrggi.arms comment records the measurement.
#
# Three adaptations are forced by a trio having one variant -- only T1 is instrumented, X
# must be positionally aligned with y, and colnames(y) must be set. All three are derived
# and measured in mrggi_feasibility.R and documented in MRGGI_METHODS.md and in the MR-GGI
# section of inference_utils.R.
#
# Writes mrggi_group_n<size>.RData per sample-size group, then inference_mrggi.RData/.csv.
#
# THIS SCRIPT NOW BUILDS A CLUSTER. It used to be single threaded on the assumption that
# MR-GGI was a handful of lm() calls per trio, which was true of the bare trio and is not
# true of the CS-alpha arm: that arm carries a median of 82-106 covariates, so MRggi
# computes ~5,800 gene pairs per trio and the arm costs 9.4 h of the 10.4 h total across
# every group. Trios are independent, so the group is spread over the cluster.

library(MRGN)
library(MRggi)
library(parallel)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

cat("=== MR-GGI ===  arms:", paste(mrggi.arms, collapse = ", "),
    "| alpha:", mrggi.alpha, "| cor.thr:", mrggi.cor.thr,
    "| min first-stage F:", mrggi.min.F, "| cores:", n.cores,
    "| sizes:", if (is.null(sample.sizes)) "all" else paste(sample.sizes, collapse = ","), "\n")

cl <- parallel::makeCluster(n.cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
    library(MRGN)
    library(MRggi)
}))

# The workers evaluate mrggi.one.dataset(), which reaches the helpers and the mrggi.*
# settings as globals. They are exported by name rather than by re-sourcing the config on
# each worker: inference_config.R reads commandArgs(), which on a PSOCK worker is the
# worker's own command line and not this process's, so a --out-dir / --sim-file run would
# silently get the defaults back AND dir.create() the default out.dir as a side effect.
#
# The list is short because the workers only ever run mrggi.run.arms() on a self-contained
# per-trio payload -- the id and selection-score columns are built on the master, so
# nothing here needs `sel`, `datasets` or the row builders.
parallel::clusterExport(cl, c(
    "mrggi.run.arms", "mrggi.one.trio", "mrggi.fields", "prefixed", "safely",
    # settings referenced as defaults inside those functions
    "mrggi.arms", "mrggi.cor.thr", "mrggi.alpha", "mrggi.min.F", "mrggi.p.adjust"))

run.method.groups(
    method = "mrggi",
    needs.selection = TRUE,      # truth/CSq/CSa arms all need the selection cache
    cl = cl,
    runner = function(datasets, sel, cov.pool, known.conf, cl) {
        run.mrggi.group(datasets, sel, colnames(cov.pool), cl = cl,
                        arms = mrggi.arms, verbose = TRUE)
    })

cat("\ncombining MR-GGI:\n")
invisible(combine.method("mrggi", sample.sizes = sample.sizes))
