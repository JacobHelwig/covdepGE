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
cavi_search <- function(X, Z, D, y, hp_method, ssq, sbsq, pip, nssq, nsbsq,
                        npip, ssq_upper_mult, ssq_lower, snr_upper, sbsq_lower,
                        pip_lower, elbo_tol, alpha_tol, max_iter, resp_index){

  # get the dimensions of the data
  n <- nrow(X)
  p <- ncol(X)

  # instantiate initial values of variational parameters; the l, j entry is
  # the variational approximation to the j-th parameter in a regression with
  # the resp_index predictor fixed as the response with weightings taken with
  # respect to the l-th individual
  alpha <- matrix(0.1, n, p)
  mu <- matrix(0, n, p)

  # if the hyperparameter grid has not been fully supplied, create it
  if (any(is.null(c(ssq, sbsq, pip)))){

    # if either sbsq or pip has not been supplied, then use LASSO to estimate
    # the proportion of non-zero coefficients
    if (any(is.null(c(sbsq, pip)))){
      lasso <- glmnet::cv.glmnet(X, y)

      # find the number of non-zero coefficients estimated by LASSO
      # ensure that non0 is an integer in (0, p)
      non0 <- sum(coef(lasso, s = "lambda.1se")[-1] != 0)
      non0 <- min(max(non0, 1), p - 1)

      # an upper bound for pi is the proportion of non-zero coefficients
      pi_upper <- non0 / p
    }

    # if ssq has not been supplied, create the grid; otherwise, find the unique
    # values of the supplied ssq
    if (is.null(ssq)){

      # find the upper bound for the grid of ssq
      ssq_upper <- ssq_upper_mult * var(y)

      # create the grid candidates for ssq
      # ssq <- exp(seq(log(ssq_lower), log(ssq_upper), length.out = nssq))
      ssq <- seq(ssq_lower, ssq_upper, length.out = nssq)

    }else{
      ssq <- unique(ssq)
    }

    # if sbsq has not been supplied, create the grid; otherwise, find the unique
    # values of the supplied sbsq
    if (is.null(sbsq)){

      # find the sum of the variances for each of the columns of X_mat
      s2_sum <- sum(apply(X, 2, var))

      # find the upper bound for the grid of sbsq
      sbsq_upper <- snr_upper / (pi_upper * s2_sum)

      # create the grid candidates for sbsq
      # sbsq <- exp(seq(log(sbsq_lower), log(sbsq_upper), length.out = nsbsq))
      sbsq <- seq(sbsq_lower, sbsq_upper, length.out = nsbsq)

    }else{
      sbsq <- unique(sbsq)
    }

    # if pip has not been supplied, create the grid; otherwise, find the unique
    # values of the supplied pip
    if (is.null(pip)){

      # create posterior inclusion probability grid
      pip <- exp(seq(log(pip_lower), log(pi_upper), length.out = npip))
      # pip <- seq(pip_lower, pi_upper, length.out = npip)

    }else{
      pip <- unique(pip)
    }

    # create the grid
    hp <- expand.grid(pip = pip, ssq = ssq, sbsq = sbsq, elbo = NA)

  }

  if (hp_method == "grid_search"){

    # perform grid search to select the best hyperparameters

    # perform CAVI for each of the hyperparameter settings
    out_grid <- grid_search_c(y, D, X, mu, alpha, hp$ssq, hp$sbsq, hp$pip,
                              elbo_tol, alpha_tol, max_iter)

    # add the ELBO for each of the grid points to the hyperparameter grid
    hp$elbo <- out_grid$elbo_vec

    # use the best hyperparameters from the grid search to perform CAVI until
    # max_iter is reached or alpha converges (no early stopping for ELBO)
    out <- cavi_c(y, D, X, out_grid$mu, out_grid$alpha, out_grid$ssq,
                  out_grid$sbsq, out_grid$pip, -1, alpha_tol, max_iter)

    # save hyperparameter details
    final <- c(ssq = out_grid$ssq, sbsq = out_grid$sbsq, pip = out_grid$pip)
    hp <- list(grid = hp, grid_sz = nrow(hp), final = final)

    # save final variational parameters and elbo
    alpha  <- out$alpha
    mu <- out$mu
    ssq_var <- out$ssq_var
    elbo <- out$elbo

  }else{

    # otherwise, hp_method is model_average or hybrid

    # if hybrid, re-define the grid so that it is the cartesian product between
    # the variance hyperparameters
    if (hp_method == "hybrid"){
      pip <- unique(hp$pip)
      ssq <- unique(hp$ssq)
      sbsq <- unique(hp$sbsq)
      hp <- expand.grid(ssq = ssq, sbsq = sbsq)

      # data.frame for storing final hyperparameters from grid search
      hyp <- data.frame(ssq = NA, sbsq = NA, pip = pip, elbo = NA)
    }

    # get number of parameters to average over; if hybrid, only averaging over
    # pi; otherwise, average over entire grid
    n_hp <- ifelse(hp_method == "hybrid", length(pip), nrow(hp))

    # create lists for storing the variational parameters and elbos for each
    # hyperparameter setting
    elbo_theta <- alpha_theta <- mu_theta <- ssqv_theta <- vector("list", n_hp)

    # iterate over each of the pi if hybrid
    # iterate over each of the hyperparameter settings otherwise
    for (j in 1:n_hp){

      # hybrid CAVI
      if (hp_method == "hybrid"){

        # fix the j-th value of the pip; repeat for each entry in the
        # hyperparameter grid
        pip_j <- rep(pip[j], nrow(hp))

        # grid search over the ssq and sbsq grid
        out_grid <- grid_search_c(y, D, X, mu, alpha, hp$ssq, hp$sbsq, pip_j,
                                  elbo_tol, alpha_tol, max_iter)

        # use the best hyperparameters from the grid search to perform CAVI until
        # max_iter is reached or alpha converges (no early stopping for ELBO)
        out <- cavi_c(y, D, X, out_grid$mu, out_grid$alpha, out_grid$ssq,
                      out_grid$sbsq, pip[j], -1, alpha_tol, max_iter)

        # save the final hyperparameters and elbo
        hyp[j, c("ssq", "sbsq", "elbo")] <- c(out_grid$ssq, out_grid$sbsq, out$elbo)

        # fix the final hyperparameter setting
        hp_j <- data.frame(ssq = out_grid$ssq, sbsq = out_grid$sbsq, pip = pip[j])

      }else{

        # otherwise, model averaging CAVI

        # fix hyperparameter setting
        hp_j <- hp[j, ]

        # perform CAVI for the hyperparameter setting
        out <- cavi_c(y, D, X, mu, alpha, hp_j$ssq, hp_j$sbsq, hp_j$pip, elbo_tol,
                      alpha_tol, max_iter)
      }

      # save the variational parameters
      alpha_theta[[j]] <- out$alpha
      mu_theta[[j]] <- out$mu
      ssqv_theta[[j]] <- out$ssq_var

      # calculate the elbo for each individual under the current hyperparameter
      # setting and save
      elbo_l <- rep(NA, n)
      for (l in 1:n){
        elbo_l[l] <- ELBO_calculator_c(y, D[ , l], X, t(out$ssq_var[l, ]),
                                       t(out$mu[l, ]), t(out$alpha[l, ]),
                                       hp_j$ssq, hp_j$sbsq, hp_j$pip)
      }

      # save the ELBO for all individuals
      elbo_theta[[j]] <- elbo_l

      # sum the ELBO across the individuals to get the elbo for the
      # hyperparameter setting
      hp$elbo[j] <- sum(elbo_l)
    }

    # calculate weights for averaging and average

    # vector for saving the average elbo for each individual
    elbo_avg <- rep(NA, n)

    # calculate the weights for each individual
    for (l in 1:n){

      # retrieve the elbo for the l-th individual for each hyperparameter
      # setting
      elbo_l <- sapply(elbo_theta, `[`, l)

      # subtract the maximum value from each and exponentiate to get the lower
      # bound for the likelihood (subtracting max for numerical stability)
      lik_l <- exp(elbo_l - max(elbo_l))

      # normalize the weights
      weight <- lik_l / sum(lik_l)

      # calculate the average elbo for this individual
      elbo_avg[l] <- sum(elbo_l * weight)

      # multiply the l-th row of each of the variational parameter matrices by
      # the corresponding weight
      for (k in 1:length(alpha_theta)){
        alpha_theta[[k]][l, ] <- weight[k] * alpha_theta[[k]][l, ]
        mu_theta[[k]][l, ] <- weight[k] * mu_theta[[k]][l, ]
        ssqv_theta[[k]][l, ] <- weight[k] * ssqv_theta[[k]][l, ]
      }
    }

    # sum across the weighted variational parameter matrices to obtain the final
    # variational parameter matrices
    alpha <- Reduce("+", alpha_theta)
    mu <- Reduce("+", mu_theta)
    ssq_var <- Reduce("+", ssqv_theta)

    # sum across the ELBO to obtain the final ELBO
    elbo <- sum(elbo_avg)

    # save hyperparameter details
    grid_sz <- nrow(hp) * ifelse(hp_method == "hybrid", length(pip), 1)
    hp <- list(grid = hp, grid_sz = grid_sz, final = "<NA> if hp_method == 'model_average'")

    # if hybrid, add the final hyperparameters chosen by grid search to
    # the hyperparameter details
    if(hp_method == "hybrid") hp$final <- hyp
  }

  # return the final variational parameters, hyperparameter details, and elbo
  return(list(alpha_matrix = alpha, mu_matrix = mu, ssqv_matrix = ssq_var,
              hyperparameters = hp, elbo = elbo))
}
