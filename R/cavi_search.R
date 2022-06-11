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
    ssq <- matrix(mean(idmod$sigma), nrow(sbsq), p + 1)
    pip <- matrix(mean(probs), nrow(sbsq), p + 1)
  }

  # if the hyperparameter grid has not been fully supplied, create it
  if (any(is.null(c(ssq, sbsq, pip)))){

    # if either sbsq or pip has not been supplied, then use LASSO to estimate
    # the proportion of non-zero coefficients
    if (any(is.null(c(sbsq, pip)))){
      lasso <- glmnet::cv.glmnet(X, y)
      non0 <- sum(coef(lasso, s = "lambda.min")[-1] != 0)
      non0 <- max(non0, 1)
      pi_hat <- non0 / p
    }

    # if ssq has not been supplied, create the grid; otherwise, find the unique
    # values of the supplied ssq
    if (is.null(ssq)){

      # find the upper bound for the grid of ssq
      ssq_upper <- ssq_upper_mult * var(y)

      # create the grid candidates for ssq
      ssq <- exp(seq(log(var_lower), log(ssq_upper), length.out = nssq))

    }else{
      ssq <- unique(ssq[ , resp_index])
    }

    # if sbsq has not been supplied, create the grid; otherwise, find the unique
    # values of the supplied sbsq
    if (is.null(sbsq)){

      # find the sum of the variances for each of the columns of X_mat
      s2_sum <- sum(apply(X, 2, var))

      # find the upper bound for the grid of sbsq
      sbsq_upper <- 25 / (pi_hat * s2_sum)

      # create the grid candidates for ssq and sbsq
      sbsq <- exp(seq(log(var_lower), log(sbsq_upper), length.out = nsbsq))

    }else{
      sbsq <- unique(sbsq[ , resp_index])
    }

    # if pip has not been supplied, create the grid; otherwise, find the unique
    # values of the supplied pip
    if (is.null(pip)){

      # create posterior inclusion probability grid
      pip <- exp(seq(log(0.01), log(min(0.9, 2 * max(pi_hat, 0.05))), length.out = npip))


    }else{
      pip <- unique(pip[ , resp_index])
    }

    # create the grid
    hp <- expand.grid(pip = pip, ssq = ssq, sbsq = sbsq)

  }else{

    # otherwise, the user has supplied the hyperparameters; take the resp_index
    # column of each as the current hyperparameters
    hp <- data.frame(pip = pip[ , resp_index], ssq = ssq[ , resp_index],
                     sbsq = sbsq[ , resp_index])

  }

  # perform grid search to select the best hyperparameters
  if (!(grid_search %in% c("model_average", "hybrid"))){

    # perform CAVI for each of the hyperparameter settings
    # use grid search to select hyperparameters
    out_grid <- grid_search_c(y, D, X, mu, alpha, hp$ssq, hp$sbsq, hp$pip,
                              elbo_tol, alpha_tol, max_iter, T)

    # add the ELBO for each of the grid points to the hyperparameter grid
    hp$elbo <- out_grid$elbo_vec

    # use the best hyperparameters from the grid search to perform a proper
    # CAVI until max_iter is reached or alpha converges
    out <- cavi_c(y, D, X, out_grid$mu, out_grid$alpha, out_grid$ssq,
                  out_grid$sbsq, out_grid$pip, elbo_tol, alpha_tol, max_iter, F)


    # save CAVI details
    cavi_details <- list(ELBO = out$elbo,
                         iterations = out$converged_iter + out_grid$iterations,
                         converged = out$converged_iter < max_iter)

    # save hyperparameter details
    hyp <- list(ssq = out_grid$ssq, sbsq = out_grid$sbsq, pip = out_grid$pip,
                grid = hp, grid_sz = nrow(hp))

    # save progress of alpha and elbo
    prog <- list(alpha = out$alpha_prog, elbo = out$elbo_prog)

    return(list(alpha_matrix = out$alpha, mu_matrix = out$mu,
                cavi_details = cavi_details, hyperparameters = hyp,
                progress = prog))

  }else if(grid_search == "model_average"){

    # apply importance sampling to average over hyperparameters

    # create lists for storing the alpha matrices, mu matrices and elbos
    elbo_theta <- alpha_theta <- mu_theta <- vector("list", nrow(hp))

    # iterate over each of the hyperparameter settings
    for (j in 1:nrow(hp)){

      # fix hyperparameter setting
      hp_j <- hp[j, ]

      out <- cavi_c(y, D, X, mu, alpha, hp_j$ssq, hp_j$sbsq, hp_j$pip, elbo_tol,
                    alpha_tol, max_iter, F)

      # save the alpha and mu matrices
      alpha_theta[[j]] <- out$alpha
      mu_theta[[j]] <- out$mu

      # calculate the elbo for each individual under the current hyperparameter
      # setting and save
      elbo_l <- rep(NA, n)
      for (l in 1:n){
        elbo_l[l] <- ELBO_calculator_c(y, D[ , l], X, t(out$ssq_var[l, ]),
                                       t(out$mu[l, ]), t(out$alpha[l, ]),
                                       hp_j$ssq, hp_j$sbsq, hp_j$pip)
      }

      # save the ELBO
      elbo_theta[[j]] <- elbo_l
    }

    # vector for saving the average elbo for each individual
    elbo_avg <- rep(NA, n)

    # matrix for saving importance weights
    weights <- NA
    if (p < 10) weights <- matrix(NA, n, nrow(hp))

    # calculate the weights for each individual
    for (l in 1:n){

      # retrieve the elbo for the l-th individual for each hyperparameter
      # setting
      elbo_l <- sapply(elbo_theta, `[`, l)

      # subtract the maximum value from each and exponentiate to get the lower
      # bound for the likelihood
      lik_l <- exp(elbo_l - max(elbo_l))

      # normalize the weights
      weight <- lik_l / sum(lik_l)

      # save the weights
      if (p < 10) weights[l, ] <- weight

      # calculate the average elbo for this individual
      elbo_avg[l] <- sum(elbo_l * weight)

      # multiply the l-th row of each of the alpha and mu matrices by the
      # corresponding weight
      for (k in 1:length(alpha_theta)){
        alpha_theta[[k]][l, ] <- weight[k] * alpha_theta[[k]][l, ]
        mu_theta[[k]][l, ] <- weight[k] * mu_theta[[k]][l, ]
      }
    }

    # sum across the alpha and mu matrices to obtain the final alpha and
    # mu matrices
    alpha <- Reduce("+", alpha_theta)
    mu <- Reduce("+", mu_theta)

    # sum across the ELBO to obtain the final ELBO
    elbo <- sum(elbo_avg)

    return(list(alpha_matrix = alpha, mu_matrix = mu,
                cavi_details = list(converged = T, ELBO = elbo),
                hyperparameters = list(hyperparameters = hp, weights = weights),
                progress = list()))
  }else{

    # apply hybrid method

    # re-define the hyperparameters so that the grid is cartesian product between
    # the variance hyperparameters
    pip <- unique(hp$pip)
    ssq <- unique(hp$ssq)
    sbsq <- unique(hp$ssq)
    hp <- expand.grid(ssq = ssq, sbsq = sbsq)

    # create lists for storing the alpha matrices, mu matrices and elbos
    elbo_theta <- alpha_theta <- mu_theta <- vector("list", length(pip))

    # data.frame for storing the final hyperparameters
    hyp <- data.frame(ssq = NA, sbsq = NA, pip = pip, elbo = NA)

    # iterate over each of the pip
    for (j in 1:length(pip)){

      # fix the j-th value of the pip; repeat for each entry in the hyperparameter
      # grid
      pip_j <- rep(pip[j], nrow(hp))

      # grid search over the ssq and sbsq grid
      out_grid <- grid_search_c(y, D, X, mu, alpha, hp$ssq, hp$sbsq, pip_j,
                                elbo_tol, alpha_tol, max_iter, T)

      # use the best hyperparameters from the grid search to perform a proper
      # CAVI until max_iter is reached or alpha converges
      out <- cavi_c(y, D, X, out_grid$mu, out_grid$alpha, out_grid$ssq,
                    out_grid$sbsq, pip[j], elbo_tol, alpha_tol, max_iter, F)

      # save the final hyperparameters and elbo
      hyp[j, c("ssq", "sbsq", "elbo")] <- c(out_grid$ssq, out_grid$sbsq, out$elbo)

      # save the alpha and mu matrices
      alpha_theta[[j]] <- out$alpha
      mu_theta[[j]] <- out$mu

      # calculate the elbo for each individual under the current hyperparameter
      # setting and save
      elbo_l <- rep(NA, n)
      for (l in 1:n){
        elbo_l[l] <- ELBO_calculator_c(y, D[ , l], X, t(out$ssq_var[l, ]),
                                       t(out$mu[l, ]), t(out$alpha[l, ]),
                                       out_grid$ssq, out_grid$sbsq, pip[j])
      }

      # save the ELBO
      elbo_theta[[j]] <- elbo_l
    }

    # vector for saving the average elbo for each individual
    elbo_avg <- rep(NA, n)

    # matrix for saving importance weights
    weights <- NA
    if (p < 10) weights <- matrix(NA, n, length(pip))

    # calculate the weights for each individual
    for (l in 1:n){

      # retrieve the elbo for the l-th individual for each hyperparameter
      # setting
      elbo_l <- sapply(elbo_theta, `[`, l)

      # subtract the maximum value from each and exponentiate to get the lower
      # bound for the likelihood
      lik_l <- exp(elbo_l - max(elbo_l))

      # normalize the weights
      weight <- lik_l / sum(lik_l)

      # save the weights
      if (p < 10) weights[l, ] <- weight

      # calculate the average elbo for this individual
      elbo_avg[l] <- sum(elbo_l * weight)

      # multiply the l-th row of each of the alpha and mu matrices by the
      # corresponding weight
      for (k in 1:length(alpha_theta)){
        alpha_theta[[k]][l, ] <- weight[k] * alpha_theta[[k]][l, ]
        mu_theta[[k]][l, ] <- weight[k] * mu_theta[[k]][l, ]
      }
    }

    # sum across the alpha and mu matrices to obtain the final alpha and
    # mu matrices
    alpha <- Reduce("+", alpha_theta)
    mu <- Reduce("+", mu_theta)

    # sum across the ELBO to obtain the final ELBO
    elbo <- sum(elbo_avg)

    return(list(alpha_matrix = alpha, mu_matrix = mu,
                cavi_details = list(converged = T, ELBO = elbo),
                hyperparameters = list(hyperparameters = hyp, weights = weights),
                progress = list()))

  }
}
