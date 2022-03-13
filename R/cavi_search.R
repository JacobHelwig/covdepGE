## -----------------------------------------------------------------------------
## -----------------------------cavi_search-------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## Performs CAVI and grid search for a fixed data matrix and response for n
## linear regressions, where the l-th regression is weighted with respect to the
## l-th individual. Returns matrix of posterior inclusion probabilities, details
## on the variational updates, a vector of warnings, and the number of final
## CAVI that DNC
## -----------------------------ARGUMENTS---------------------------------------
## X_mat: n x p matrix; predictors
##
## Z: n x p' matrix; extraneous covariates
##
## D: n x n matrix; weights (k,l entry is the weight of the k-th individual
## with respect to the l-th individual using the l-th individual's bandwidth)
##
## y: vector of length n; response
##
## alpha: numeric in [0, 1]; global initialization value for the variational
## parameters alpha_matrices (approximates probabilities of inclusion). 0.2 by
## default
##
## mu: numeric; global initialization value for the variational parameters
## mu_matrices (approximates posterior mean of regression coefficients). 0 by
## default
##
## sigmasq_vec: vector of length n_param with positive entries; candidate
## values of sigmasq, the error term variance. NULL by default
##
## sigmabetasq_vec: vector of length n_param with positive entries; candidate
## values of sigmabeta_sq, the slab variance. NULL by default
##
## pi_vec: vector of length n_param with entries in (0, 1); candidate values of
## pi. 0.1 by default
##
## tolerance: positive numeric; end CAVI when the Frobenius norm of the
## iteration-to-iteration change in the alpha matrix are within tolerance.
## `1e-12` by default
##
## max_iter_grid: positive integer; during the grid search, if the
## tolerance criteria has not been met by `max_iter_grid` iterations, end the
## CAVI. `1e4` by default
##
## max_iter_final: positive integer; for the final CAVI, if the tolerance
## criteria has not been met by `max_iter_final` iterations, end the CAVI. `1e4`
## by default
##
## warnings: logical; if T, convergence and grid warnings will be
## displayed. Convergence warnings occur when the tolerance exit condition has
## not been met by max_iter_grid or max_iter_final iterations. Grid warnings
## occur when, for either sigmabetasq_vec or pi_vec, the grid is longer than 2
## candidates, and the final CAVI uses a candidate value on the grid boundary.
## T by default
##
## resp_index: integer in {1,...,p + 1}; the index of the column that is y in
## data_mat
##
## CS: logical; if T, pi_vec and sigma_sq will be selected
## according to Carbonetto-Stephens. F by default
## -----------------------------RETURNS-----------------------------------------
## Returns `list` with the following values:
##
## 1. alpha_matrix: n x p matrix; the l, j entry is the variational
## approximation to the posterior inclusion probability of the j-th variable in
## a regression with the y fixed as the response with weightings taken with
#  respect to the l-th individual
##
## 2. CAVI_details: list with the following values:
##  - sigmasq, sigmabetasq, pi: numerics; the values of the hyperparameters
##  that maximized the ELBO for the j-th variable
##  - ELBO: numeric; the maximum value of ELBO for the final CAVI
##  - converged_iter: integer; the number of iterations to attain convergence
##  for the final CAVI
##  - hyperparameters: n_param x 4 matrix; each of the hyperparameter grid
##  points with the resulting ELBO
##
## 3. warnings_vec: character vector; Vector of convergence warnings to be
## displayed in covdepGE_main
##
## 4. final_DNC: integer; number of final CAVIs that did not converge
##
## 5. sigmasq: n x (p + 1) matrix; fitted error term variances for each
## individual in the final model. Column j corresponds to the regression with
## the j-th variable fixed as the response.
##
## 6. sigmabeta_sq: n x (p + 1) matrix; fitted slab variances for each
## individual in the final model. Column j corresponds to the regression with
## the j-th variable fixed as the response
## -----------------------------------------------------------------------------
cavi_search <- function(X, Z, D, y, alpha, mu, ssq, sbsq, pip, nssq, nsbsq, npip,
                        ssq_upper_mult, var_lower, elbo_tol, alpha_tol,
                          max_iter, warnings, resp_index, CS, R, grid_search){

  # get the dimensions of the data and hyperparameter candidates
  n <- nrow(X)
  p <- ncol(X)

  # instantiate initial values of variational parameters; the l, j entry is
  # the variational approximation to the j-th parameter in a regression with
  # the resp_index predictor fixed as the response with weightings taken with
  # respect to the l-th individual
  alpha <- matrix(alpha, n, p)
  mu <- matrix(mu, n, p)

  # get hyperparameter values from varbvs
  if (CS){
    set.seed(resp_index)
    idmod <- varbvs::varbvs(X, y, Z = Z[ , 1], verbose = FALSE)
    probs <- 1 / (1 + exp(-idmod$logodds)) # need to convert to log base 10
    ssq <- rep(mean(idmod$sigma), length(sbsq))
    pip <- rep(mean(probs), length(sbsq))
  }

  # if the hyperparameter grid has not been supplied, create it
  if (any(is.null(c(ssq, sbsq, pip)))){

    # find the upper bound for the grid of ssq
    ssq_upper <- ssq_upper_mult * var(y)

    # find the lasso estimate to the proportion of non-zero coefficients
    lasso <- glmnet::cv.glmnet(X, y)
    non0 <- sum(coef(lasso, s = "lambda.min")[-1] != 0)
    non0 <- max(non0, 1)
    pi_hat <- mean(non0) / p

    # find the sum of the variances for each of the columns of X_mat
    s2_sum <- sum(apply(X, 2, var))

    # find the upper bound for the grid of sbsq
    sbsq_upper <- 25 / (pi_hat * s2_sum)

    # create the grid candidates for ssq and sbsq
    ssq <- exp(seq(log(var_lower), log(ssq_upper), length.out = nssq))
    sbsq <- exp(seq(log(var_lower), log(sbsq_upper), length.out = nsbsq))

    # create posterior inclusion probability grid
    pip <- seq(0.05, 0.45, length.out = npip)

    # create the grid
    hp <- expand.grid(pip = pip, ssq = ssq, sbsq = sbsq)
  }else{
    hp <- data.frame(pip = pip, ssq = ssq, sbsq = sbsq)
  }

  # loop to optimize hyperparameters; run CAVI for each grid points; store the
  # resulting ELBO
  if (R){
    out <- grid_search_R(y, D, X, mu, alpha, hp$ssq, hp$sbsq, hp$pip, elbo_tol,
                         alpha_tol, max_iter, grid_search)
    best_pip <- out$pip
    # use the best hyperparameters from the grid search to perform a proper CAVI
    # until max_iter is reached or alpha converges
    out <- cavi_R(y, D, X, out$mu, out$alpha, out$ssq, out$sbsq, best_pip,
                  elbo_tol, alpha_tol, max_iter, F)
  }else{
    out <- grid_search_c(y, D, X, mu, alpha, hp$ssq, hp$sbsq, hp$pip, elbo_tol,
                         alpha_tol, max_iter, grid_search)
    best_pip <- out$pip
    # use the best hyperparameters from the grid search to perform a proper CAVI
    # until max_iter is reached or alpha converges
    out <- cavi_c(y, D, X, out$mu, out$alpha, out$ssq, out$sbsq, best_pip,
                  elbo_tol, alpha_tol, max_iter, F)
  }

  # save CAVI details
  cavi_details <- list(ELBO = out$elbo, iterations = out$converged_iter,
                       converged = out$converged_iter < max_iter)

  # save hyperparameter details
  hyp <- list(ssq = out$ssq, sbsq = out$sbsq, pip = best_pip, pip_cands = pip,
              ssq_cands = ssq, sbsq_cands = sbsq, grid_sz = nrow(hp))

  return(list(alpha_matrix = out$alpha, mu_matrix = out$mu,
              cavi_details = cavi_details, hyperparameters = hyp))
}
