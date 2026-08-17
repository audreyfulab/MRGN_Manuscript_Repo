library(MRGN)
library(mvtnorm)
library(gridExtra)

# set simulation conditions
# -----------------------------
# seed
set.seed(234)
# number of replicates of each scenario
# (previously 300, dropping to 100 now with more scenarios)
number_of_replicates <- 100

effect_sizes <- c(small = c(0.1, 0.3), medium = c(0.3, 0.5), large = c(0.5, 1))

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

same.effects <- function(scenarios) {
    
    for(effect_scene in c("small", "medium", "large")){
        # get the indices of the scenarios with this effect size
        idx <- which(scenarios$effect_size == effect_scene)
        # sample b.snp from a uniform distribution over the range of effect sizes for this scenario
        scenarios[idx, ]$b.snp <- sample(
            seq(
                from = effect_sizes[[effect_scene]][1],
                to = effect_sizes[[effect_scene]][2],
                by = 0.05),
            size = length(idx), 
            replace = TRUE
        )
        # sample b.med from the same range as b.snp
        scenarios[idx, ]$b.med <- sample(
            seq(
                from = effect_sizes[[effect_scene]][1],
                to = effect_sizes[[effect_scene]][2],
                by = 0.05),
            size = length(idx), 
            replace = TRUE
        )
    }
    return(scenarios)
}

scenarios$SD <- 1 #just use 1 for now (same as Yang et al. 2017)

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


simulate.dataset <- function(settings, clinical.covs, verbose = TRUE) {

    # simulate dataset based on settings

    # initalize
    init=0
    attempts=0
    # number of each type of confounder for this scenario
    K_n <- settings$K_n
    U_n <- settings$U_n
    W_n <- settings$W_n
    Z_n <- settings$Z_n

    # simulate confounders from multivariate normal distribution
    Uvars <- rmvnorm(n = settings$sample.size,
                     mean = rep(0, U_n),
                     sigma = diag(U_n))

    # simData.from.graph() fills the K and U nodes by matching the node names generated
    # by gen.graph.skel() ("K1".."K3", "U1".."Ud") against colnames(conf.mat). An
    # unnamed conf.mat matches nothing, and it fails two different ways: as a matrix it
    # silently returns all-NA confounder columns that propagate into T1/T2 as NaN, and
    # as a data.frame it errors with "undefined columns selected". These names are
    # load bearing.
    conf.mat <- as.data.frame(Uvars)
    colnames(conf.mat) <- paste0("U", 1:U_n)
    # if sample size is 670, include clinical covariates in the confounders (original sim)
    if (K_n > 0) {
        conf.mat <- cbind.data.frame(clinical.covs$WholeBlood, conf.mat)
        colnames(conf.mat) <- c(paste0("K", 1:K_n), paste0("U", 1:U_n))
    }
    row.names(conf.mat) <- 1:nrow(conf.mat)

    if(verbose){
        print(paste0("printing params for sim trio:"))
        print(settings)
    }
    # simulate dataset until all genotypes are represented
    while(length(init)<=2){
        attempts = attempts + 1
        # simulate dataset using simData.from.graph function
        sim_data = simData.from.graph(model = settings$model,
                            theta = settings$minor.freq,
                            b0.1 = 0,
                            b.snp = settings$b.snp,
                            b.med = settings$b.med,
                            sd.1 = settings$SD,
                            conf.num.vec = c(K = K_n, U = U_n,
                                                W = W_n,
                                                Z = Z_n),
                            simulate.confs = FALSE,
                            conf.mat = conf.mat,
                            plot.graph = FALSE,
                            # all four entries must stay in this list even when a block
                            # is empty: gen.graph.skel indexes it positionally (K,U,W,Z)
                            conf.coef.ranges = list(K = c(0, 0),
                                                    U = c(0.05, 0.5), # follows from actual pc effect range
                                                    W = c(0.05, 0.5), # follows from actual pc effect range
                                                    Z = c(1, 1.5))) #GMAC setting
        init=unique(sim_data$data$V1)
    }
    cat("resampled", attempts, "times before 3 genotypes were represented\n")
    return(list(data = sim_data$data, attempts = attempts))
}


name.trio.columns <- function(data, index, K_n, kc.names) {
    # match the column naming of the original simulation: the K confounders carry the
    # clinical covariate names, and every simulated confounder (U, W, Z) is suffixed
    # with the index of the dataset it came from, so the columns stay unique once the
    # confounders of all datasets are pooled into a single matrix
    nms <- colnames(data)
    if (K_n > 0) {
        nms[3 + seq_len(K_n)] <- kc.names
    }
    sim.conf.idx <- seq(4 + K_n, length(nms))
    nms[sim.conf.idx] <- paste0(nms[sim.conf.idx], ".", index)
    colnames(data) <- nms
    return(data)
}


conf.r.squared <- function(data, K_n, U_n) {
    # proportion of the variance of T1 and T2 explained by the U confounders: the
    # simulated analogue of regressing a GTEx gene on its selected PCs. Recorded so the
    # realized confounding strength can be checked against the real data, where the
    # selected PCs explain a median R2 of 0.41 (cis) and 0.31 (trans) in Whole Blood.
    # Only the U block is used -- W is an intermediate and Z a common child of the
    # trio, so regressing on those would not be the same quantity.
    if (U_n > (nrow(data) - 2)) {
        # saturated fit (at least as many confounders as observations): no usable R2
        return(c(R2.T1.U = NA_real_, R2.T2.U = NA_real_))
    }
    U.block <- data[, 3 + K_n + seq_len(U_n), drop = FALSE]
    c(R2.T1.U = summary(lm(data[, 2] ~ ., data = U.block))$r.squared,
      R2.T2.U = summary(lm(data[, 3] ~ ., data = U.block))$r.squared)
}


write.rdata <- function(data, filename) {
    if (is.null(filename)) {
        path <- file.path(getwd(), "simulated_trios.RData")
    } else {
        path <- filename
    }
    save(data, file = path)
}

simulate.all.datasets <- function(scenarios, clinical.covs, verbose = TRUE, save = FALSE, filename = NULL) {
    # simulate datasets for all scenarios
    datasets <- vector("list", nrow(scenarios))
    cis.conf.effects <- NULL
    trans.conf.effects <- NULL

    for (i in 1:nrow(scenarios)) {
        settings <- scenarios[i, ]
        sim <- simulate.dataset(settings, clinical.covs, verbose = verbose)
        # every parameter that generated this dataset, plus what the draw realized
        params <- cbind.data.frame(dataset = i,
                                   settings,
                                   n.resamples = sim$attempts,
                                   as.list(conf.r.squared(sim$data,
                                                          K_n = settings$K_n,
                                                          U_n = settings$U_n)),
                                   row.names = NULL)

        cis.conf.effects = summary(lm(T1 ~ ., data = sim$data[, -c(1,3)]))$coefficients
        trans.conf.effects = summary(lm(T2 ~ ., data = sim$data[, -c(1,2)]))$coefficients

        datasets[[i]] <- list(data = name.trio.columns(sim$data, index = i,
                                                       K_n = settings$K_n,
                                                       kc.names = kc.names),
                              params = params,
                              conf.effects = list(cis = cis.conf.effects,
                                                  trans = trans.conf.effects))
    }
    if (save) {
        write.rdata(datasets, filename)
    }
    return(datasets)
}

root = getwd()
setwd("./bioinfo_revision/simulated_data/")
result <- simulate.all.datasets(
    scenarios,
    clinical.covs,
    verbose = TRUE,
    save = TRUE,
    filename = "simulated_trios.RData"
    )
setwd(root)

