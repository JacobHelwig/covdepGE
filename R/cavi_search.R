## -----------------------------------------------------------------------------
## -----------------------------cavi_search---------------------------------------
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
cavi_search <- function(X_mat, Z, D, y, alpha, mu, sigmasq_vec, update_sigmasq,
                        sigmabetasq_vec, update_sigmabetasq, pi_vec, tolerance,
                        max_iter_grid, max_iter_final, warnings, resp_index, CS,
                        R = F){

  # get the dimensions of the data and hyperparameter candidates
  n <- nrow(X_mat)
  p <- ncol(X_mat)
  n_param <- length(pi_vec)

  if (length(unique(as.numeric(sigmasq_vec))) == 1){
    sigmasq_vec <- matrix(var(y), n, n_param)
  }

  # instantiate initial values of variational parameters; the l, j entry is
  # the variational approximation to the j-th parameter in a regression with
  # the resp_index predictor fixed as the response with weightings taken with
  # respect to the l-th individual
  alpha_mat <- matrix(alpha, n, p)
  mu_mat <- matrix(mu, n, p)

  # get hyperparameter values from varbvs
  set.seed(resp_index)
  if (CS){
    idmod <- varbvs::varbvs(X_mat, y, Z = Z[ , 1], verbose = FALSE)
    probs <- 1 / (1 + exp(-idmod$logodds)) # need to convert to log base 10
    sigmasq_vec <- matrix(mean(idmod$sigma), n, ncol(sigmabetasq_vec))
    pi_vec <- rep(mean(probs), ncol(sigmabetasq_vec))
  }

  # loop to optimize sigmabeta_sq; run CAVI for each grid points; store the
  # resulting ELBO
  if (R){
    out <- grid_search_R(y, D, X_mat, mu_mat, alpha_mat, sigmasq_vec,
                         update_sigmasq,  sigmabetasq_vec, update_sigmabetasq,
                         pi_vec, tolerance, max_iter_grid)
  } else{
    grid_search_out <- grid_search_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq_vec,
                                     sigmabetasq_vec, pi_vec, tolerance,
                                     max_iter_grid)
  }

  # if (update_sigmasq | update_sigmabetasq){
  #
  #   # if at least one of the variance hyperparameters is being updated via MAPE,
  #   # repeat the grid search with the initial values of the hyperparameter(s)
  #   # being updated as the final values from the first grid search, as well as
  #   # the initial values for the variational parameters as the final values from '
  #   # the first grid search
  #   sigmasq_vec <- matrix(sigmasq, n, n_param)
  #   sigmabetasq_vec <- matrix(sigmabetasq_vec, n, n_param)
  #   grid_search_out <- grid_search_R(y, D, X_mat, mu_mat, alpha_mat, sigmasq_vec,
  #                                    update_sigmasq, sigmabetasq_vec,
  #                                    update_sigmabetasq, pi_vec, tolerance,
  #                                    max_iter_grid)
  # }

  # save CAVI details
  cavi_details <- list(ELBO = out$elbo, iterations = out$iterations,
                       converged = out$iterations < max_iter_grid)

  # save hyperparameter details
  hyperparameters <- list(sigmasq = out$sigmasq, sigmabeta_sq = out$sigmabeta_sq,
                          pi = out$pi, pi_grid = pi_vec)

  # alpha is an n by p matrix; the i,j-th entry is the probability of inclusion
  # for the i-th individual for the j-th variable according to the regression on
  # y
  alpha_matrix <- out$alpha

  return(list(alpha_matrix = alpha_matrix, cavi_details = cavi_details,
              hyperparameters = hyperparameters))
}
