#######################################
# GMAC method for confounder selection
#######################################
gmac_get_conf <- function(cl=NULL, cov.pool=NULL, exp.dat, snp.dat.cis, trios.idx, fdr=0.05, fdr_filter=0.1) {
    
  triomatrix <- array(NA, c(dim(exp.dat)[2], dim(trios.idx)[1], 3))
  
  print ("organizing data")
  for (i in 1:dim(trios.idx)[1]) {
    triomatrix[, i, ] <- cbind(round(snp.dat.cis[trios.idx[i, 1], ], digits = 0),
                               exp.dat[trios.idx[i, 2], ], exp.dat[trios.idx[i, 3], ])
  }
  
  ## obtain pc confounder candidates
  num_trio <- dim(triomatrix)[2]
  
  ## use cov.pool covariates
  pool_cov <- t(cov.pool)
  ## obtain pool confounder candidates and estimate pool confounders
  num_pool <- dim(pool_cov)[2]
  print ("obtain pool confounder candidates")
  if (is.null(cl)) {
    p_value_child_pool <- matrix(unlist(lapply(1:num_trio, child.p, tripletmatrix = triomatrix,
                                               covariates = pool_cov), use.names = F),
                                 byrow = TRUE, ncol = num_pool)
  } else {
    p_value_child_pool <- matrix(unlist(parLapply(cl, 1:num_trio, child.p,
                                                  tripletmatrix = triomatrix, covariates = pool_cov),
                                        use.names = F), byrow = TRUE, ncol = num_pool)
  }
  q_child_pool <- qvalue(p = as.vector(p_value_child_pool), fdr.level = fdr_filter)$qvalues
  conf_candidates_pool <- matrix(0, nrow = nrow(p_value_child_pool), ncol = num_pool)
  conf_candidates_pool[which(q_child_pool >= fdr_filter)] <- 1
  
  
  ## Estimate pool confounders
  print ("estimate pool confounders")
  if (is.null(cl)) {
    est_conf_pool_idx <- matrix(unlist(lapply(1:num_pool, conf.fdr, tripletmatrix = triomatrix,
                                              covariates = pool_cov, conf_candidates = conf_candidates_pool,
                                              fdr = fdr), use.names = F), ncol = num_pool)
  } else {
    est_conf_pool_idx <- matrix(unlist(parLapply(cl, 1:num_pool, conf.fdr,
                                                 tripletmatrix = triomatrix, covariates = pool_cov,
                                                 conf_candidates = conf_candidates_pool, fdr = fdr), use.names = F), ncol = num_pool)
  }
  
  print ("complete")
  output <- est_conf_pool_idx
  return(output)
}


## This function uses stratified fdr to figure out, for each locus, the
## list of covariates that do not play roles as child/intermediate mediator.
child.p <- function(i, tripletmatrix, covariates) {
  treatment <- tripletmatrix[, i, 1]
  n_obs <- dim(tripletmatrix)[1]
  n_cov <- dim(covariates)[2]
  p_value_child <- rep(NA, n_cov)
  cov_length <- dim(covariates)[1]
  if (cov_length == (n_obs + 1)) {
    covariates <- covariates[-cov_length, ]
  }
  for (j in 1:n_cov) {
    p_value_child[j] <- summary(lm(covariates[, j] ~ treatment))$coef[2, 4]
  }
  return(p_value_child)
}

## This function uses stratified fdr to figure out, for each covariate,
## the list of trios where the covariate plays a role as a confounder.
conf.fdr <- function(i, tripletmatrix, covariates, conf_candidates, fdr) {
  n_obs <- dim(tripletmatrix)[1]
  n_tri <- dim(tripletmatrix)[2]
  p_value <- rep(NA, n_tri)
  cov <- covariates[, i]
  cov_length <- length(cov)
  str_fdr <- rep(0, n_tri)
  candidate_trio_id <- sort(which(conf_candidates[, i] == 1))
  if (length(cov) == (n_obs + 1)) {
    pi0 = 1 - cov[cov_length]
    cov <- cov[-cov_length]
    for (j in 1:n_tri) {
      f <- summary(lm(cov ~ tripletmatrix[, j, 2] + tripletmatrix[, j, 3]))$fstatistic
      p_value[j] <- pf(f[1], f[2], f[3], lower.tail = F)
    }
    
    q <- pi0 * p_value[candidate_trio_id] * length(candidate_trio_id)/rank(p_value[candidate_trio_id])
    str_fdr[candidate_trio_id[which(q <= fdr)]] <- 1
  } else if (cov_length == n_obs) {
    for (j in 1:n_tri) {
      #print (j)
      f <- summary(lm(cov ~ tripletmatrix[, j, 2] + tripletmatrix[, j, 3]))$fstatistic
      p_value[j] <- pf(f[1], f[2], f[3], lower.tail = F)
    }
    
    q <- qvalue(p = p_value[candidate_trio_id], fdr.level = fdr)$qvalues
    str_fdr[candidate_trio_id[which(q <= fdr)]] <- 1
  }
  return(str_fdr)
}

##qvalue is the function to estimate q values given a vector of p values. The code for this function is extracted from the "qvalue" package written by Alan Dabney and John Storey.
# We include the functions of qvalue here since it is reported than some operating system cannot isntall it, and we are using version 2.8.0.

qvalue <- function(p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, pi0 = NULL, ...) {
  # Argument checks
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  } else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    stop("'fdr.level' must be in (0, 1].")
  }

  # Calculate pi0 estimate
  if (is.null(pi0)) {
    pi0s <- pi0est(p, ...)
  } else {
    if (pi0 > 0 && pi0 <= 1)  {
    pi0s = list()
    pi0s$pi0 = pi0
    } else {
      stop("pi0 is not (0,1]")
    }
  }

  # Calculate q-value estimates
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m / (i * (1 - (1 - p[o]) ^ m))))[ro]
  } else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m /i ))[ro]
  }
  qvals_out[rm_na] <- qvals
  # Calculate local FDR estimates
  if (lfdr.out) {
    lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  } else {
    lfdr_out <- NULL
  }

  # Return results
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level,
                   significant = (qvals <= fdr.level),
                   pi0.lambda = pi0s$pi0.lambda, lambda = pi0s$lambda,
                   pi0.smooth = pi0s$pi0.smooth)
  } else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
                   pvalues = p_in, lfdr = lfdr_out, pi0.lambda = pi0s$pi0.lambda,
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}

##pi0est is a function to estimates the proportion of true null p-values. The code for this function is extracted from the "qvalue" package written by Alan Dabney and John Storey.
pi0est <- function(p, lambda = seq(0.05,0.95,0.05), pi0.method = c("smoother", "bootstrap"),
                   smooth.df = 3, smooth.log.pi0 = FALSE, ...) {
  # Check input arguments
  rm_na <- !is.na(p)
  p <- p[rm_na]
  pi0.method = match.arg(pi0.method)
  m <- length(p)
  lambda <- sort(lambda) # guard against user input

  ll <- length(lambda)
  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  } else if (ll > 1 && ll < 4) {
    stop(cat("ERROR:", paste("length(lambda)=", ll, ".", sep=""),
             "If length of lambda greater than 1,",
             "you need at least 4 values."))
  } else if (min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }
  # Determines pi0
  if (ll == 1) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0.lambda <- pi0
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL
  } else {
    pi0 <- sapply(lambda, function(l) mean(p >= l) / (1 - l))
    pi0.lambda <- pi0
    # Smoother method approximation
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      } else {
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    } else if (pi0.method == "bootstrap") {
      # Bootstrap method closed form solution by David Robinson
      minpi0 <- quantile(pi0, prob = 0.1)
      W <- sapply(lambda, function(l) sum(p >= l))
      mse <- (W / (m ^ 2 * (1 - lambda) ^ 2)) * (1 - W / m) + (pi0 - minpi0) ^ 2
      pi0 <- min(pi0[mse == min(mse)], 1)
      pi0Smooth <- NULL
    } else {
      stop('ERROR: pi0.method must be one of "smoother" or "bootstrap".')
    }
  }
  if (pi0 <= 0) {
    stop("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.")
  }
  return(list(pi0 = pi0, pi0.lambda = pi0.lambda,
              lambda = lambda, pi0.smooth = pi0Smooth))
}

##lfdr is a function to estimate the local FDR values from p-values. The code for this function is extracted from the "qvalue" package written by Alan Dabney and John Storey.
lfdr <- function(p, pi0 = NULL, trunc = TRUE, monotone = TRUE,
                 transf = c("probit", "logit"), adj = 1.5, eps = 10 ^ -8, ...) {
  # Check inputs
  lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("P-values not in valid range [0,1].")
  } else if (is.null(pi0)) {
    pi0 <- pi0est(p, ...)$pi0
  }
  n <- length(p)
  transf <- match.arg(transf)
  # Local FDR method for both probit and logit transformations
  if (transf == "probit") {
    p <- pmax(p, eps)
    p <- pmin(p, 1 - eps)
    x <- qnorm(p)
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    lfdr <- pi0 * dnorm(x) / y
  } else {
    x <- log((p + eps) / (1 - p + eps))
    myd <- density(x, adjust = adj)
    mys <- smooth.spline(x = myd$x, y = myd$y)
    y <- predict(mys, x)$y
    dx <- exp(x) / (1 + exp(x)) ^ 2
    lfdr <- (pi0 * dx) / y
  }
  if (trunc) {
    lfdr[lfdr > 1] <- 1
  }
  if (monotone) {
    o <- order(p, decreasing = FALSE)
    ro <- order(o)
    lfdr <- cummax(lfdr[o])[ro]
  }
  lfdr_out[rm_na] <- lfdr
  return(lfdr_out)
}
