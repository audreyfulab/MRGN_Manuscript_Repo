# Confounder-structure simulations: what MRGN does when the trio's covariate structure
# changes.
#
# NOTE: Set your working directory to the repository root before running this script.
# e.g., setwd("path/to/MRGN_Manuscript_Repo")
#
#   Rscript bioinfo_revision/simulation/confounder_structure_simulation.R
#
# updated_data_simulation.R only ever generates ONE structure: confounders (U) plus one
# intermediate (W) plus one common child (Z). Reviewers asked what happens under the other
# three, so this script generates them:
#
#   u_only   confounders only                       W_n = 0, Z_n = 0
#   u_w      confounders + 1 intermediate           W_n = 1, Z_n = 0
#   u_z      confounders + 1 common child           W_n = 0, Z_n = 1
#   (u_w_z)  confounders + both -- ALREADY COVERED by simulated_trios.RData, not repeated
#
# n = 670 ONLY, and MRGN only. 670 is where the K block of real clinical covariates exists,
# and it is the sample size the manuscript's Table 2 expansion is about. Each case is
# 5 models x 3 effect sizes x 20 replicates = 300 trios, matching one sample-size group of
# the main simulation exactly, so the four structures are directly comparable.
#
# Everything else is held identical to updated_data_simulation.R -- the same effect-size
# strata, the same minor allele frequency draw, the same U_n range, the same residual SD,
# the same conf.coef.ranges. The ONLY thing that varies between the three cases is W_n and
# Z_n. That is the point: a difference in the results has one possible cause.
#
# simData.from.graph() accepts empty W and Z blocks -- verified for all four combinations of
# (W_n, Z_n) in {0,1}^2 -- so simulation_utils.R needs no change and is reused unmodified.
#
# Writes three files into simulated_data/:
#   simulated_trios_u_only.RData
#   simulated_trios_u_w.RData
#   simulated_trios_u_z.RData
#
# Inference is run by simulation_results/run_structure_sims.R, which points apply_mrgn.R at
# each of these with --sim-file / --out-dir / --filter-int-child.

library(MRGN)
library(mvtnorm)

source("./bioinfo_revision/simulation/simulation_utils.R")

# One seed per case, all distinct from updated_data_simulation.R's 234, so the three cases
# are independent draws rather than the same trios with columns removed. Sharing a seed
# would make the U blocks identical across cases and turn a comparison between structures
# into a paired test that the analysis does not treat as one.
CASES <- list(
    u_only = list(W_n = 0, Z_n = 0, seed = 2341,
                  file = "simulated_trios_u_only.RData",
                  label = "confounders only"),
    u_w    = list(W_n = 1, Z_n = 0, seed = 2342,
                  file = "simulated_trios_u_w.RData",
                  label = "confounders + 1 intermediate"),
    u_z    = list(W_n = 0, Z_n = 1, seed = 2343,
                  file = "simulated_trios_u_z.RData",
                  label = "confounders + 1 common child"))

number_of_replicates <- 20
SAMPLE.SIZE <- 670

# identical to updated_data_simulation.R; see the long note there for why the strata are
# contiguous, continuous, and why b.snp gets the wider range
effect_sizes <- list(
    b.snp = list(small  = c(0.1, 0.5),
                 medium = c(0.5, 1.0),
                 large  = c(1.0, 1.5)),
    b.med = list(small  = c(0.1, 0.3),
                 medium = c(0.3, 0.5),
                 large  = c(0.5, 1.0)))

clinical.covs <- loadRData("./GTEx/data/kclist_top5_tiss.RData")
kc.sample.size <- nrow(clinical.covs$WholeBlood)
kc.names <- colnames(clinical.covs$WholeBlood)

if (SAMPLE.SIZE != kc.sample.size) {
    stop("SAMPLE.SIZE is ", SAMPLE.SIZE, " but the clinical covariates have ",
         kc.sample.size, " rows. The K block only exists at the Whole Blood donor count; ",
         "these cases are defined at that sample size.")
}


build.scenarios <- function(case) {
    scenarios <- expand.grid(model = c("model0", "model1", "model2", "model3", "model4"),
                             sample.size = SAMPLE.SIZE,
                             effect_size = c("small", "medium", "large"),
                             replicate = 1:number_of_replicates,
                             stringsAsFactors = FALSE)
    n <- nrow(scenarios)

    scenarios$minor.freq <- sample(seq(0.01, 0.5, 0.01), n, replace = TRUE)
    scenarios <- draw.effect.sizes(scenarios, effect_sizes)
    scenarios$SD <- 1

    # the whole point of this script
    scenarios$W_n <- case$W_n
    scenarios$Z_n <- case$Z_n
    scenarios$K_n <- length(kc.names)

    # same range and same guard as the main simulation: the confounder count must leave
    # room for T1, T2, G and the intercept
    scenarios$U_n <- pmin(sample(1:50, n, replace = TRUE),
                          scenarios$sample.size - scenarios$Z_n - scenarios$W_n - 4)
    scenarios
}


root <- getwd()
for (nm in names(CASES)) {
    case <- CASES[[nm]]
    cat("\n=== ", nm, ": ", case$label,
        "  (W_n = ", case$W_n, ", Z_n = ", case$Z_n, ") ===\n", sep = "")

    set.seed(case$seed)
    scenarios <- build.scenarios(case)
    cat("  ", nrow(scenarios), " trios at n = ", SAMPLE.SIZE, "\n", sep = "")

    # simulate.all.datasets() writes relative to the working directory, as in
    # updated_data_simulation.R
    setwd("./bioinfo_revision/simulation/simulated_data/")
    result <- simulate.all.datasets(
        scenarios,
        clinical.covs,
        kc.names = kc.names,
        verbose = FALSE,
        save = TRUE,
        filename = case$file)
    setwd(root)

    cols <- colnames(result[[1]]$data)
    cat("  wrote ", case$file, " | first trio columns: ",
        paste(head(cols, 6), collapse = ", "), " ... ",
        paste(tail(cols, 2), collapse = ", "), "\n", sep = "")
    cat("  W columns present: ", any(grepl("^W[0-9]", cols)),
        " | Z columns present: ", any(grepl("^Z[0-9]", cols)), "\n", sep = "")
}
cat("\ndone. Run inference with:",
    "\n  Rscript bioinfo_revision/simulation_results/run_structure_sims.R\n")
