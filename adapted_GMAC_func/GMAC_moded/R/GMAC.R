#' Genomic Mediation analysis with Adaptive Confounding adjustment
#'
#' The gmac function performs genomic mediation analysis with adaptive confounding adjustment. It tests for mediation effects for a set of user specified mediation trios (e.g., eQTL, cis- and trans-genes) in the genome with the assumption of the presence of cis-association. The gmac function considers either a user provided pool of potential confounding variables, real or constructed by other methods, or all the PCs based on the expression data as the potential confounder pool. It returns the mediation p-values and the proportions mediated (e.g., the percentage of reduction in trans-effects after accounting for cis-mediation), based on the mediation tests i) adjusting for known confounders only, and ii) adjusting for known confounders and adaptively selected potential confounders for each mediation trio. It also provides plots of mediation p-values (in the negative of log base of 10) versus the proportions mediated based on the above two adjustments.
#'
#' In genomic studies, a large number of mediation tests are often performed, and it is challenging to adjust for unmeasured confounding effects for the cis- and trans-genes (i.e., mediator-outcome) relationship. The current function adaptively selects the variables to adjust for each mediation trio given a large pool of constructed or real potential confounding variables. The function allows the input of variables known to be potential cis- and trans-genes (mediator-outcome) confounders in all mediation tests (\code{known.conf}), and the input of the pool of candidate confounders from which potential confounders for each mediation test will be adaptively selected (\code{cov.pool}). When no pool is provided (\code{cov.pool = NULL}), all the PCs based on expression data (\code{exp.dat}) will be constructed as the potential confounder pool.
#'
#' The algorithm assumes the presence of cis-association (treatment-mediator association), random eQTL (treatment) and the standard identification assumption in causal mediation literature that no effect of eQTL (treatment) that confounds the cis- and trans-genes (mediator-outcome) relationship. The algorithm will first filter out common child (Figure 1.B) and intermediate variables (Figure 1.C) from \code{cov.pool} for each mediation trio at a pre-specified significance threshold of FDR (\code{fdr_filter}) by utilizing their associations with the eQTL (treatment). Then, confounder (Figure 1.A) set for each mediation trio will be selected from the retained pool of candidate variables using a stratified FDR approach. Specifically, for each trio, the p-values of association for each candidate variable to the cis-gene (mediator) and trans-gene (outcome) pairs are obtained based on the F-test for testing the joint association to either the cis-gene (mediator) or the trans-gene (outcome). For each candidate variable, a pre-specified FDR (\code{fdr}) threshold is applied to the p-values corresponding to the joint associations of this variable to all the potential mediation trios. Lastly, mediation is tested for each mediation trio. Adjusting for the adaptively selected confounder set, we calculate the mediation statistic as the Wald statistic for testing the indirect mediation effect \eqn{H_0: \beta_1 = 0} based on the regression \eqn{T_j = \beta_0+\beta_1 C_i+\beta_2 L_i + \tau X_{ij}+\epsilon} where \eqn{L_i}, \eqn{C_i}, \eqn{T_j} and \eqn{X_{ij}} are the eQTL genotype (treatment), the cis-gene expression level (mediator), the trans-gene expression level (outcome) and the selected set of potential confounding variables. P-values are calculated based on within-genotype group permutation on the cis-gene expression level which maintains the cis- and trans-associations while breaks the potential mediation effect from the cis- to the trans-gene transcript.
#'
#'\if{html}{\figure{Figure1.png}{options: width=5in}} \if{latex}{\figure{Figure1.png}{options: width=5in}}
#'
#'Figure 1. Graphical illustrations of (A) a potential mediation relationship among an eQTL \eqn{L_i}, its cis-gene transcript \eqn{C_i}, and a trans-gene transcript \eqn{T_j}, with confounders \eqn{X_{ij}}(i.e., variables affecting both \eqn{C_i} and \eqn{T_j}), allowing \eqn{L_i} to affect \eqn{T_j}  via a pathway independent of \eqn{C_i}. For the mediation effect tests to have a causal interpretation, adjustment must be made for the confounders. (B) A potential mediation trio with common child variables, \eqn{Z_{ij}} (i.e., variables affected by both \eqn{C_i} and \eqn{T_j}). Adjusting for common child variables in mediation analysis would ``marry" \eqn{C_i} and \eqn{T_j} and make \eqn{C_i} appearing to be regulating \eqn{T_j} even if there is no such effect. (C) A potential mediation trio with intermediate variables \eqn{W_{ij}} (i.e., variables affected by \eqn{C_i} and affecting \eqn{T_j}). Adjusting for intermediate variables in mediation analysis would prevent the detection of the true mediation effect from \eqn{C_i} to \eqn{T_j}.
#'
#'The algorithm returns the mediation p-values (\code{pvals}) and the proportions mediated (\code{beta.change}, i.e., the percentage of reduction in trans-effects after accounting for cis-mediation), based on the mediation tests i) adjusting for known confounders only, and ii) adjusting for known confounders and adaptively selected potential confounders for each mediation trio. It also returns indicator matrix for the selected potential confounders (\code{sel.conf.ind}). Plots of mediation p-values (in the negative of log base of 10) versus the proportions mediated based on the adjustments i) and ii) are provided. The plot could further be used as a diagnostic check for sufficiency in confounding adjustment in scenarios such as cis-gene mediating trans-gene regulation pattern, where we expect the trios with very significant mediation p-values to have positive proportions mediated. Therefore, a J shape pattern is expected when most if not all confounding effects have been well adjusted, whereas a U shape pattern may indicate the presence of unadjusted confounders.
#'
#' @param cl Parallel backend if it is set up. It is used for parallel computing.
#' @param known.conf A known confounders matrix which is adjusted in all mediation tests. Each row is a confounder, each column is a sample.
#' @param cov.pool The pool of candidate confounding variables from which potential confounders are adaptively selected to adjust for each mediation test. Each row is a covariate, each column is a sample.
#' @param exp.dat A gene expression matrix. Each row is for one gene, each column is a sample.
#' @param snp.dat.cis The cis-eQTL genotype matrix. Each row is an eQTL, each column is a sample.
#' @param trios.idx A matrix of selected trios indexes (row numbers) for mediation tests. Each row consists of the index (i.e., row number) of the eQTL in \code{snp.dat.cis}, the index of cis-gene transcript in \code{exp.dat}, and the index of trans-gene in \code{exp.dat}. The dimension is the number of trios by three.
#' @param nperm The number of permutations for testing mediation.
#' @param fdr The false discovery rate to select confounders. We set \code{fdr}=0.05 as default.
#' @param fdr_filter The false discovery rate to filter common child and intermediate variables. We set \code{fdr_filter}=0.1 as default.
#' @param nominal.p An option to obtain the nominal p-value or permutation-based p-value, which is the default.
#' @return The algorithm will return a list of p-values, beta changes, and indicator matrix for confounders selected.
#' \item{pvals}{The mediation p-values. A matrix with dimension of the number of trios by two ("Adjust Known Covariates Only", "Adjust Known + Selected Covariates").}
#' \item{beta.change}{The proportions mediated. A matrix with dimension of the number of trios by two ("Adjust Known Covariates Only", "Adjust Known + Selected Covariates").}
#' \item{sel.conf.ind}{An indicator matrix with dimension of the number of trios by the number of covariates in \code{cov.pool} or \code{pc.matrix} if the principal components (PCs) based on expression data are used as the pool of potential confounders.}
#' \item{pc.matrix}{PCs will be returned if the PCs based on expression data are used as the pool of potential confounders. Each column is a PC.}
#' @references Fan Yang, Jiebiao Wang, the GTEx consortium, Brandon L. Pierce, and Lin S. Chen. Identifying cis-mediators for trans-eQTLs across many human tissues using genomic mediation analysis. Genome Research. (Pending acceptance). \doi{10.1101/078683}
#' @export
#' @importFrom qvalue qvalue
#' @importFrom parallel parLapply
#' @importFrom stats lm pf pnorm sd prcomp
#' @importFrom graphics par plot
#' @examples data(example)
#'
#' # a fast example with only 100 permutations
#' output <- gmac(known.conf=dat$known.conf, cov.pool=dat$cov.pool, exp.dat=dat$exp.dat, snp.dat.cis=dat$snp.dat.cis,
#'                    trios.idx=dat$trios.idx[1:40,], nperm=100, nominal.p=TRUE)
#'
#' plot(output)
#'
#'
#' \dontrun{
#' ## the construction of PCs as cov.pool
#' pc <- prcomp(t(dat$exp.dat), scale = T)
#' cov.pool <- t(pc$x)
#'
#'
#' ## generate a cluster with 2 nodes for parallel computing
#' cl <- makeCluster(2)
#' output <- gmac(cl=cl, known.conf=dat$known.conf, cov.pool=cov.pool, exp.dat=dat$exp.dat, snp.dat.cis=dat$snp.dat.cis,
#' trios.idx=dat$trios.idx, nominal.p=TRUE)
#' stopCluster(cl)
#' }

gmac <- function(cl = NULL, known.conf, cov.pool=NULL, exp.dat, snp.dat.cis, trios.idx,
                     nperm=10000, fdr=0.05, fdr_filter=0.1, nominal.p = FALSE) {

  known_confounders <- t(known.conf)

  triomatrix <- array(NA, c(dim(exp.dat)[2], dim(trios.idx)[1], 3))
  print ("organizing data")
  for (i in 1:dim(trios.idx)[1]) {
    triomatrix[, i, ] <- cbind(round(snp.dat.cis[trios.idx[i, 1], ], digits = 0),
                               exp.dat[trios.idx[i, 2], ], exp.dat[trios.idx[i, 3], ])
  }

  ## obtain pc confounder candidates
  num_trio <- dim(triomatrix)[2]

  use.PC = F
  ## if cov.pool not provided, construct pc as cov.pool
  if(is.null(cov.pool)) {
    use.PC = T
    ## normalize expression data so that every gene contributes equally to
    ## the construction of PCs
    exp.dat2 <- t(apply(exp.dat, 1, function(x) (x - mean(x))/sd(x)))

    pc <- prcomp(t(exp.dat2))
    all_pc <- pc$x
    eigen <- pc$sdev^2
    expvar <- eigen/sum(eigen)
    cov.pool <- t(rbind(all_pc, expvar))
  }

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
  #print(p_value_child_pool)
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
  ## Final mediation analysis for each trio by within group permutation of mediators
  print("Running final trio analysis")
  if (is.null(cl)) {
    start = Sys.time()
    output <- lapply(1:num_trio, getp.func, triomatrix = triomatrix,
                     known_confounders = known_confounders, pool_cov = pool_cov,
                     est_conf_pool_idx = est_conf_pool_idx, nperm = nperm,
                     nominal.p = nominal.p, use.PC = use.PC)
    end = Sys.time()
    time.to.compute.trios = c(start = start, end = end)

  } else {
    # the timing has to be recorded on this branch too: comp.time is read unconditionally
    # when the output list is assembled below, so a parallel run errored with
    # "object 'time.to.compute.trios' not found"
    start = Sys.time()
    output <- parLapply(cl, 1:num_trio, getp.func, triomatrix = triomatrix,
                        known_confounders = known_confounders, pool_cov = pool_cov,
                        est_conf_pool_idx = est_conf_pool_idx, nperm = nperm,
                        nominal.p = nominal.p, use.PC = use.PC)
    end = Sys.time()
    time.to.compute.trios = c(start = start, end = end)
  }
  pvals <- matrix(unlist(lapply(output, function(x) x$pvals), use.names = FALSE),
                  byrow = T, ncol = 2)
  beta.change <- matrix(unlist(lapply(output, function(x) x$beta.change),
                               use.names = FALSE), byrow = T, ncol = 2)

  colnames(pvals) <- colnames(beta.change) <- c("Known", "Known_sel_pool")
  print("done!")
  if(use.PC){
	output <- list(pvals=pvals, beta.change = beta.change, pc.matrix = all_pc, sel.conf.ind = est_conf_pool_idx)
  } else{
		output <- list(pvals=pvals, beta.change = beta.change, sel.conf.ind = est_conf_pool_idx,
		               comp.time = time.to.compute.trios)
  }


  class(output) <- 'GMAC'
  return(output)
}


#' @export
plot.GMAC <- function(x, ...) {
  Pvalues_mediation <- x$pvals
  beta.change <- x$beta.change

  fig.title=c("Adjust Known Covariates Only", "Adjust Known + Selected Covariates")
  par(mfrow=c(1,2))
  for(k in 1:2) plot(y=-log10(Pvalues_mediation[,k]+1e-16),x=beta.change[,k],cex.main=0.8, cex.lab=0.8,
                     pch=20, main=fig.title[k],xlim=c(-3,3),xlab="Proportion of beta change",
                     ylab="-log10(P-value of mediation)")
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
      f <- summary(lm(cov ~ tripletmatrix[, j, 2] + tripletmatrix[,
                                                                  j, 3]))$fstatistic
      p_value[j] <- pf(f[1], f[2], f[3], lower.tail = F)
    }

    q <- pi0 * p_value[candidate_trio_id] * length(candidate_trio_id)/rank(p_value[candidate_trio_id])
    str_fdr[candidate_trio_id[which(q <= fdr)]] <- 1
  } else if (cov_length == n_obs) {
    for (j in 1:n_tri) {
      f <- summary(lm(cov ~ tripletmatrix[, j, 2] + tripletmatrix[,
                                                                  j, 3]))$fstatistic
      p_value[j] <- pf(f[1], f[2], f[3], lower.tail = F)
    }
    #print(p_value)
    q <- qvalue(p = p_value[candidate_trio_id], fdr.level = fdr)$qvalues
    str_fdr[candidate_trio_id[which(q <= fdr)]] <- 1
  }
  return(str_fdr)
}

my.solve <- function(X) {
  if (!is.matrix(X))
    X <- matrix(X, nrow = sqrt(length(X)))
  ss <- svd(X)
  Xinv <- ss$u %*% diag(1/ss$d, nrow = nrow(X), ncol = nrow(X)) %*% t(ss$v)
  return(Xinv)
}


## Function: indirect effect
Indirect <- function(mediator, treatment, outcome, confounderset, return.t.only = TRUE) {
  n = length(treatment)
  x <- cbind(mediator, treatment, 1, confounderset)
  p_x <- dim(x)[2]
  inverse_xx <- my.solve(t(x) %*% x)
  beta <- inverse_xx %*% t(x) %*% outcome
  var_beta <- as.numeric(1/(n - p_x) * (sum(outcome^2) - t(outcome) %*%
                                          x %*% inverse_xx %*% t(x) %*% outcome)) * inverse_xx
  t_stat <- beta[1]/sqrt(var_beta[1, 1])
  if (return.t.only) {
    return(t_stat)
  } else {
    x2 = cbind(treatment, 1, confounderset)
    inverse_xx2 <- my.solve(t(x2) %*% x2)
    beta.total <- inverse_xx2 %*% t(x2) %*% outcome
    return(list(t_stat = t_stat, beta = beta[2], beta.total = beta.total[1]))
  }
}


nominal.pfun <- function(stat, stat0) {
  2 * (1 - pnorm(abs((stat - mean(stat0))/sd(stat0))))
}


get.beta.change <- function(beta_direct, beta_total) {
  beta_change <- (beta_total - beta_direct)/beta_total
  return(beta_change)
}


getp.func <- function(i, triomatrix, known_confounders, pool_cov, est_conf_pool_idx, nperm, nominal.p, use.PC) {
  treatment <- triomatrix[, i, 1]
  mediator <- triomatrix[, i, 2]
  outcome <- triomatrix[, i, 3]
  if(use.PC) pool_cov = pool_cov[-nrow(pool_cov),]
  print(which(est_conf_pool_idx[i, ] == 1))
  if(length(which(est_conf_pool_idx[i, ] == 1)) == 0){
    cov_known_sel_pool <- known_confounders
  }else{
    cov_known_sel_pool <- cbind(known_confounders, pool_cov[, which(est_conf_pool_idx[i, ] == 1)])
  }
  
  mediator_perm <- matrix(rep(mediator, each = nperm), byrow = TRUE, ncol = nperm)
  for (j in 0:2) {
    ind <- which(treatment == j)
    if (length(ind) > 1) {
      mediator_perm[ind, ] <- apply(mediator_perm[ind, ], 2, sample)
    }
  }

  indirect.known <- Indirect(mediator = mediator, treatment = treatment,
                             outcome = outcome, confounderset = known_confounders, return.t.only = FALSE)
  t_known <- indirect.known$t_stat
  beta_change_known <- get.beta.change(beta_direct = indirect.known$beta,
                                       beta_total = indirect.known$beta.total)

  indirect.sel.pool <- Indirect(mediator = mediator, treatment = treatment,
                                outcome = outcome, confounderset = cov_known_sel_pool, return.t.only = FALSE)
  t_known_sel_pool <- indirect.sel.pool$t_stat
  beta_change_sel_pool <- get.beta.change(beta_direct = indirect.sel.pool$beta,
                                          beta_total = indirect.sel.pool$beta.total)

  t_perm_known <- apply(mediator_perm, 2, Indirect, treatment = treatment,
                        outcome = outcome, confounderset = known_confounders)
  t_perm_known_sel_pool <- apply(mediator_perm, 2, Indirect, treatment = treatment,
                                 outcome = outcome, confounderset = cov_known_sel_pool)

  if (nominal.p) {
    pvalue_known <- nominal.pfun(stat = t_known, stat0 = t_perm_known)
    pvalue_known_sel_pool <- nominal.pfun(stat = t_known_sel_pool,
                                          stat0 = t_perm_known_sel_pool)
  } else {
    pvalue_known <- mean(abs(t_known) <= abs(t_perm_known))
    pvalue_known_sel_pool <- mean(abs(t_known_sel_pool) <=
                                    abs(t_perm_known_sel_pool))
  }

  pvals <- c(pvalue_known, pvalue_known_sel_pool)
  beta.change <- c(beta_change_known, beta_change_sel_pool)

  return(list(pvals = pvals, beta.change = beta.change))
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

