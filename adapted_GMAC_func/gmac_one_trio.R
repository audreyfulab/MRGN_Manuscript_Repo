###########################################################
# function to apply the GMAC method to a single trio with
# its associated confounders
#
# trio: nx3 data matrix of a trio with treatment, mediator
#       and outcome in three columns
# confounders: nxk data matrix of k confounders
###########################################################
gmacOneTrio <- function(trio, confounders, nperm=NULL, nominal.p) {
  treatment <- trio[, 1]
  mediator <- trio[, 2]
  outcome <- trio[, 3]
  #if(use.PC) pool_cov = pool_cov[-nrow(pool_cov),]
  #cov_known_sel_pool <- cbind(known_confounders, unknown_cov)
  
  # perform permutation
  mediator_perm <- matrix(rep(mediator, each = nperm), byrow = TRUE, ncol = nperm)
  for (j in 0:2) {
    ind <- which(treatment == j)
    if (length(ind) > 1) {
      mediator_perm[ind, ] <- apply(mediator_perm[ind, ], 2, sample)
    }
  }
  
  
  indirect.known <- Indirect(mediator = mediator, treatment = treatment,
                             outcome = outcome, confounderset = confounders, return.t.only = FALSE)
  t_known <- indirect.known$t_stat
  beta_change_known <- get.beta.change(beta_direct = indirect.known$beta,
                                       beta_total = indirect.known$beta.total)
  
  #indirect.sel.pool <- Indirect(mediator = mediator, treatment = treatment,
  #                              outcome = outcome, confounderset = cov_known_sel_pool, return.t.only = FALSE)
  #t_known_sel_pool <- indirect.sel.pool$t_stat
  #beta_change_sel_pool <- get.beta.change(beta_direct = indirect.sel.pool$beta,
  #                                        beta_total = indirect.sel.pool$beta.total)
  
  t_perm_known <- apply(mediator_perm, 2, Indirect, treatment = treatment,
                        outcome = outcome, confounderset = confounders)
  #t_perm_known_sel_pool <- apply(mediator_perm, 2, Indirect, treatment = treatment,
  #                               outcome = outcome, confounderset = cov_known_sel_pool)
  
  if (nominal.p) {
    pvalue_known <- nominal.pfun(stat = t_known, stat0 = t_perm_known)
  #  pvalue_known_sel_pool <- nominal.pfun(stat = t_known_sel_pool,
  #                                        stat0 = t_perm_known_sel_pool)
  } else {
    pvalue_known <- mean(abs(t_known) <= abs(t_perm_known))
  #  pvalue_known_sel_pool <- mean(abs(t_known_sel_pool) <=
  #                                  abs(t_perm_known_sel_pool))
  }
  
  #pvals <- c(pvalue_known, pvalue_known_sel_pool)
  #beta.change <- c(beta_change_known, beta_change_sel_pool)
  
  return(list(pval = pvalue_known, beta_change = beta_change_known))
}

get.beta.change <- function(beta_direct, beta_total) {
  beta_change <- (beta_total - beta_direct)/beta_total
  return(beta_change)
}

nominal.pfun <- function(stat, stat0) {
  2 * (1 - pnorm(abs((stat - mean(stat0))/sd(stat0))))
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

my.solve <- function(X) {
  if (!is.matrix(X))
    X <- matrix(X, nrow = sqrt(length(X)))
  ss <- svd(X)
  Xinv <- ss$u %*% diag(1/ss$d, nrow = nrow(X), ncol = nrow(X)) %*% t(ss$v)
  return(Xinv)
}
