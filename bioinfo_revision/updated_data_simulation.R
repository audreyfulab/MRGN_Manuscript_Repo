library(MRGN)
library(mvtnorm)
library(gridExtra)

# simulation utilities (draw.effect.sizes, simulate.dataset, name.trio.columns,
# conf.r.squared, write.rdata, simulate.all.datasets). Paths here are relative to the
# repository root, which is where this script is meant to be run from.
source("./bioinfo_revision/simulation_utils.R")

# set simulation conditions
# -----------------------------
# seed
set.seed(234)
# number of replicates of each scenario. 5 models x 5 sample sizes x 3 effect sizes x 50
# replicates = 3750 datasets (the pre-revision simulation ran 1500 at a single sample size)
number_of_replicates <- 50

# a list, not c(): c(small = c(0.1, 0.3), ...) flattens and names the elements
# "small1", "small2", ..., so effect_sizes[["small"]] would be a subscript error.
effect_sizes <- list(small = c(0.1, 0.3), medium = c(0.3, 0.5), large = c(0.5, 1))

# stringsAsFactors = FALSE matters here: gen.graph.skel() dispatches on the model with
# switch(model, model0 = ..., ...), and a factor is silently treated as its integer
# code. That happens to give the right topology only while the alphabetical level order
# lines up with the branch order; keeping model a character string removes the trap.
scenarios <- expand.grid(model =c("model0", "model1", "model2", "model3", "model4"),
                         sample.size = c(50, 150, 300, 670, 1000),
                         effect_size = c("small", "medium", "large"),
                         replicate = 1:number_of_replicates,
                         stringsAsFactors = FALSE)

# load clinical covariates for use in simulating datasets with sample size 670
clinical.covs = loadRData("./GTEx/data/kclist_top5_tiss.RData")
# the clinical covariates (pcr, platform, sex) are only observed for the 670 Whole
# Blood donors, so the K confounders can only be included at that sample size. Every
# other sample size is simulated with the K block excluded (K_n = 0).
kc.sample.size <- nrow(clinical.covs$WholeBlood)
kc.names <- colnames(clinical.covs$WholeBlood)

# count number of scenarios
n <- nrow(scenarios)

# add additional parameters to scenarios sampled from uniform distributions
scenarios$minor.freq <- sample(seq(0.01, 0.5, 0.01), n, replace = TRUE)

# draw the SNP and mediation effect sizes for every scenario. Without this, b.snp and
# b.med never make it into scenarios and simData.from.graph() fails on B[1, 2] <- 1 * b.snp
scenarios <- draw.effect.sizes(scenarios, effect_sizes)

# residual SD, fixed at 1 following Yang et al. 2017.
#
# Do not scale this by b.med the way Simulation/sim_data.R did: b.med is drawn from the
# same small/medium/large stratum as b.snp, so scaling the noise by it holds b.snp/SD
# roughly constant across strata and flattens the effect-size factor (measured: partial
# |cor(V1,T1)| of 0.55 / 0.63 / 0.48 for small / medium / large). The original script had
# no effect-size strata, so its formula does not transfer.
scenarios$SD <- 1

# number of each type of confounding variable, drawn per scenario so the whole design
# is recorded in one table.
# set number of intermediate and common child variables based on condition
scenarios$W_n <- 1
scenarios$Z_n <- 1
# K confounders are the observed clinical covariates: available at n = 670 only
scenarios$K_n <- ifelse(scenarios$sample.size == kc.sample.size,
                        length(kc.names), 0)

# U confounders are unobserved, and the number is drawn from a uniform distribution. range is adjusted
# so that the total number of confounders does not exceed the sample size minus 4 (for T1, T2, G, and the intercept).
scenarios$U_n <- pmin(sample(1:50, n, replace = TRUE),
                      scenarios$sample.size - scenarios$Z_n - scenarios$W_n - 4)


root = getwd()
setwd("./bioinfo_revision/simulated_data/")
result <- simulate.all.datasets(
    scenarios,
    clinical.covs,
    kc.names = kc.names,
    verbose = TRUE,
    save = TRUE,
    filename = "simulated_trios.RData"
    )
setwd(root)
