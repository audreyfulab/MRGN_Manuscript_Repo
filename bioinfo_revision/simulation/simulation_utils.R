# Utility functions for the revised trio simulation.
# Sourced by updated_data_simulation.R; kept free of top-level side effects so the file
# can be sourced from anywhere without running a simulation.
#
# Requires the MRGN (simData.from.graph) and mvtnorm (rmvnorm) packages, attached by
# the calling script.


# b.snp and b.med are drawn from SEPARATE ranges, both indexed by the scenario's
# effect_size stratum. The SNP effect spans (0, 1.5] and the mediation effect (0, 1],
# each split into three strata, so `effect_sizes` is a two-level list:
#
#   list(b.snp = list(small = c(lo, hi), medium = ..., large = ...),
#        b.med = list(small = c(lo, hi), medium = ..., large = ...))
#
# Two properties the strata must have, both checked below.
#
# CONTIGUOUS. Each stratum's upper bound is the next one's lower bound, so the three
# together tile the whole range with no holes. An earlier version used disjoint
# intervals -- small [0.05, 0.50], medium [0.55, 1.00] -- which made every value in
# (0.50, 0.55) unreachable by construction, an artificial hole in the middle of the
# effect-size axis.
#
# CONTINUOUS. Values are drawn with runif() over the interval rather than sampled from a
# seq(..., by = 0.05) grid. On a grid only multiples of the step exist, so a trio could
# never have b.snp = 0.51 or 0.52 no matter how the intervals line up.
#
# Neither parameter may be exactly 0: b.snp = 0 makes V1 -> T1 a null edge, so a "model0"
# dataset would carry no edges at all and its truth label M0.1 would be wrong, and
# b.med = 0 breaks model1, model2 and model4 the same way. A lower bound of 0 is allowed
# in the spec -- runif() does not return its lower bound -- but the draws are floored to
# keep that guarantee explicit rather than relying on the RNG.
draw.effect.sizes <- function(scenarios, effect_sizes) {

    if (!all(c("b.snp", "b.med") %in% names(effect_sizes))) {
        stop("effect_sizes must be a list with 'b.snp' and 'b.med' entries, each holding ",
             "small/medium/large ranges")
    }

    strata <- c("small", "medium", "large")

    # the columns have to exist before they can be filled in per effect size:
    # scenarios[idx, ]$b.snp <- ... cannot create a column that isn't there yet
    scenarios$b.snp <- NA_real_
    scenarios$b.med <- NA_real_

    for (param in c("b.snp", "b.med")) {

        rngs <- effect_sizes[[param]][strata]
        if (any(sapply(rngs, is.null))) {
            stop("effect_sizes$", param, " is missing one of small/medium/large")
        }
        if (rngs[["small"]][1] < 0) {
            stop("effect_sizes$", param, "$small starts below 0")
        }
        for (k in 1:2) {
            if (!isTRUE(all.equal(rngs[[k]][2], rngs[[k + 1]][1]))) {
                stop("effect_sizes$", param, ": ", strata[k], " ends at ", rngs[[k]][2],
                     " but ", strata[k + 1], " starts at ", rngs[[k + 1]][1],
                     " -- strata must be contiguous or values in the gap are unreachable")
            }
        }

        for (effect_scene in strata) {
            rng <- rngs[[effect_scene]]
            # indices of the scenarios in this stratum
            idx <- which(scenarios$effect_size == effect_scene)
            # b.snp and b.med are drawn independently, so within a stratum a trio can
            # have a SNP effect larger or smaller than its mediation effect
            scenarios[[param]][idx] <- pmax(runif(length(idx), rng[1], rng[2]),
                                            .Machine$double.eps)
        }
    }
    return(scenarios)
}


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
                                                    # real PC effects on cis/trans genes
                                                    # in Whole Blood run to +/- 0.20 at
                                                    # the central 95% and +/- 0.27 at the
                                                    # central 99% (see pc_distribution_
                                                    # invest/README.md); gen.conf.coefs
                                                    # draws |a| from this interval and
                                                    # flips the sign at neg.freq = 0.5.
                                                    # The old upper bound of 0.5 sat well
                                                    # past the real distribution and gave
                                                    # R2(T1 | U) = 0.65 against 0.41 in
                                                    # real data; 0.3 lands the realized
                                                    # R2 near the real value.
                                                    U = c(0, 0.3),
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


simulate.all.datasets <- function(scenarios, clinical.covs, kc.names, verbose = TRUE, save = FALSE, filename = NULL) {
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
