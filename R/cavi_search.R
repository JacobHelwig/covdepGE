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
cavi_search <- function(X, Z, D, y, alpha, mu, ssq, ssq_p, ssq_q, sbsq, sbsq_p,
                        sbsq_q, pip, pip_p, pip_q, nssq, nsbsq, npip,
                        ssq_upper_mult, var_lower, grid_search, elbo_tol,
                        alpha_tol, max_iter, warnings, resp_index, CS, R){

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
    ssq <- matrix(mean(idmod$sigma), nrow(sbsq), p + 1)
    pip <- matrix(mean(probs), nrow(sbsq), p + 1)
  }

  # if the hyperparameter grid has not been fully supplied, create it
  if (any(is.null(c(ssq, sbsq, pip)))){

    # find the upper bound for the grid of ssq
    ssq_upper <- ssq_upper_mult * var(y)

    # find the lasso estimate to the proportion of non-zero coefficients
    lasso <- glmnet::cv.glmnet(X, y)
    non0 <- sum(coef(lasso, s = "lambda.min")[-1] != 0)
    non0 <- max(non0, 1)
    pi_hat <- non0 / p

    # find the sum of the variances for each of the columns of X_mat
    s2_sum <- sum(apply(X, 2, var))

    # find the upper bound for the grid of sbsq
    sbsq_upper <- 25 / (pi_hat * s2_sum)

    # create the grid candidates for ssq and sbsq
    ssq <- exp(seq(log(var_lower), log(ssq_upper), length.out = nssq))
    sbsq <- exp(seq(log(var_lower), log(sbsq_upper), length.out = nsbsq))

    # create posterior inclusion probability grid
    pip <- exp(seq(log(0.01), log(min(0.9, 2 * pi_hat)), length.out = npip))

    # create the grid
    hp <- expand.grid(pip = pip, ssq = ssq, sbsq = sbsq)

    # take the importance density and prior density to be uniform; then, the ratio
    # of the prior to the importance density for each hyperparameter setting is
    # one
    hp$ratio <- 1

  }else{

    # otherwise, the user has supplied the hyperparameters; take the resp_index
    # column of each as the current hyperparameters
    hp <- data.frame(pip = pip[ , resp_index], ssq = ssq[ , resp_index],
                     sbsq = sbsq[ , resp_index])

    # if any of the prior or importance densities have not been specified, take
    # them to be uniform
    if (is.null(ssq_p)) ssq_p <- matrix(1, nrow(ssq), p + 1)
    if (is.null(sbsq_p)) sbsq_p <- matrix(1, nrow(sbsq), p + 1)
    if (is.null(pip_p)) pip_p <- matrix(1, nrow(pip), p + 1)
    if (is.null(ssq_q)) ssq_q <- matrix(1, nrow(ssq), p + 1)
    if (is.null(sbsq_q)) sbsq_q <- matrix(1, nrow(sbsq), p + 1)
    if (is.null(pip_q)) pip_q <- matrix(1, nrow(pip), p + 1)

    # calculate the joint prior and joint importance densities for each of the
    # hyperparameter settings
    hp_prior <- ssq_p[ , resp_index] * sbsq_p[ , resp_index] * pip_p[ , resp_index]
    hp_impt <- ssq_q[ , resp_index] * sbsq_q[ , resp_index] * pip_q[ , resp_index]

    # calculate the ratio of the joint prior density to the joint importance
    # density for each of the hyperparamter settings
    hp$ratio <- hp_prior / hp_impt
  }

  # perform CAVI for each of the hyperparameter settings
  if (grid_search){

    # use grid search to select hyperparameters
    if (R){
      out <- grid_search_R(y, D, X, mu, alpha, hp$ssq, hp$sbsq, hp$pip, elbo_tol,
                           alpha_tol, max_iter, grid_search)
    }else{
      out <- grid_search_c(y, D, X, mu, alpha, hp$ssq, hp$sbsq, hp$pip, elbo_tol,
                           alpha_tol, max_iter, grid_search)
    }

    # use the best hyperparameters from the grid search to perform a proper
    # CAVI until max_iter is reached or alpha converges
    if (R){
      out <- cavi_R(y, D, X, out$mu, out$alpha, out$ssq, out$sbsq, out$pip,
                    elbo_tol, alpha_tol, max_iter, F)
    }else{
      out <- cavi_c(y, D, X, out$mu, out$alpha, out$ssq, out$sbsq, out$pip,
                    elbo_tol, alpha_tol, max_iter, F)
    }

    alpha <- out$alpha
    mu <- out$mu

    # save CAVI details
    cavi_details <- list(ELBO = out$elbo, iterations = out$converged_iter,
                         converged = out$converged_iter < max_iter)

    # save hyperparameter details
    hyp <- list(ssq = out$ssq, sbsq = out$sbsq, pip = out$pip, pip_cands = pip,
                ssq_cands = ssq, sbsq_cands = sbsq, grid_sz = nrow(hp))

    # save progress of alpha and elbo
    prog <- list(alpha = out$alpha_prog, elbo = out$elbo_prog)

    }else{

      # apply importance sampling to average over hyperparameters

      # create lists/ vector for storing the alpha matrices and mu matrices/
      # elbo
      alpha_theta <- mu_theta <- vector("list", nrow(hp))
      elbo_theta <- rep(NA, nrow(hp))

      # iterate over each of the hyperparameter settings
      for (j in 1:nrow(hp)){

        # fix hyperparameter setting
        hp_j <- hp[j , c("ssq", "sbsq", "pip")]

        # perform CAVI for the hyperparameter setting
        if (R){
          out <- cavi_R(y, D, X, mu, alpha, hp_j$ssq, hp_j$sbsq, hp_j$pip,
                        elbo_tol, alpha_tol, max_iter, F)
        } else{
          out <- cavi_c(y, D, X, mu, alpha, hp_j$ssq, hp_j$sbsq, hp_j$pip,
                        elbo_tol, alpha_tol, max_iter, F)
        }

        # save the alpha and mu matrices and elbo
        alpha_theta[[j]] <- out$alpha
        mu_theta[[j]] <- out$mu
        elbo_theta[j] <- out$elbo
      }

      # calculate the weights for each individual
      weights <- exp(elbo_theta - max(elbo_theta)) * hp$ratio
      weights <- weights / sum(weights)

      # calculate the weighted average of the alpha and mu matrices and elbo
      alpha_theta <- lapply(1:nrow(hp), function(theta_ind) weights[theta_ind] *
                              alpha_theta[[theta_ind]])
      mu_theta <- lapply(1:nrow(hp), function(theta_ind) weights[theta_ind] *
                           mu_theta[[theta_ind]])
      elbo_theta_w <- weights * elbo_theta

      # sum across the alpha and mu matrices to obtain the final alpha and
      # mu matrices
      alpha <- Reduce("+", alpha_theta)
      mu <- Reduce("+", mu_theta)

      # sum across the ELBO to obtain the final ELBO
      elbo <- sum(elbo_theta_w)

      # save CAVI details
      cavi_details <- list(ELBO = elbo, converged = T)

      # save hyperparameter details
      hyp <- list(grid_sz = nrow(hp), weights = weights, elbo = elbo_theta,
                  hyperparameter_grid = hp)

      # save progress of alpha and elbo
      prog <- NA
    }

  return(list(alpha_matrix = alpha, mu_matrix = mu,
              cavi_details = cavi_details, hyperparameters = hyp,
              progress = prog))
}
