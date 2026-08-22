# Apply MR-GGI to the simulated trios.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation_results/apply_mrggi.R
#
# MR-GGI is a Mendelian randomisation method: it uses V1 as a genetic instrument to
# estimate the causal effect between the two genes, rather than selecting confounders and
# adjusting for them. It therefore takes NO confounder set -- the instrument is what is
# supposed to handle confounding, which is the method's premise. See MRGGI_METHODS.md.
#
# Two adaptations are required for trios, both because there is only one variant: the
# outcome gene is given no instrument (a column of zeros), and only the cis -> trans
# direction is used for edge calls. Both are derived and measured in
# mrggi_feasibility.R and documented in MRGGI_METHODS.md.
#
# Writes mrggi_group_n<size>.RData per sample-size group, then inference_mrggi.RData/.csv.
#
# Single threaded -- the estimator is a handful of lm() calls per trio -- so no cluster is
# built and run_all_inference.R gives this one core, as it does MRPC.

library(MRGN)
library(MRggi)
source("bioinfo_revision/simulation/simulation_utils.R")
source("bioinfo_revision/simulation_results/inference_config.R")
source("bioinfo_revision/simulation_results/inference_utils.R")

cat("=== MR-GGI ===  alpha:", mrggi.alpha, "| cor.thr:", mrggi.cor.thr,
    "| min first-stage F:", mrggi.min.F,
    "| sizes:", if (is.null(sample.sizes)) "all" else paste(sample.sizes, collapse = ","), "\n")

run.method.groups(
    method = "mrggi",
    needs.selection = FALSE,      # MR-GGI uses the instrument, not a confounder set
    cl = NULL,
    runner = function(datasets, sel, cov.pool, known.conf, cl) {
        run.mrggi.group(datasets, verbose = TRUE)
    })

cat("\ncombining MR-GGI:\n")
invisible(combine.method("mrggi", sample.sizes = sample.sizes))
