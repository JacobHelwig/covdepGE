## -----------------------------------------------------------------------------
## -----------------------------cavi--------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## Performs CAVI and hyperparameter selection for n linear regressions, where
## the l-th regression is weighted with respect to the l-th individual
## -----------------------------ARGUMENTS---------------------------------------
## X: n x p numeric matrix; data with the j-th column removed
##
## Z: n x q numeric matrix; extraneous covariates
##
## D: n x n numeric matrix; i, j entry is the weighting of the i-th individual
## with respect to the j-th individual using the j-th individual's bandwidth
##
## y: numeric vector of length n; j-th column of the data that is fixed as the
## reponse
##
## hp_method character in c("grid_search", "model_average", "hybrid"); method
## for setting hyperparameter values based on the hyperparameter grid.
##
## ssq: NULL OR numeric vector with positive entries; candidate values
## of the hyperparameter sigma^2 (prior residual variance)
##
## sbsq: NULL OR numeric vector with positive entries; candidate values
## of the hyperparameter sigma^2_beta (prior slab variance)
##
## pip: NULL OR numeric vector with entries in (0, 1); candidate
## values of the hyperparameter pi (prior inclusion probability)
##
## nssq: positive integer; number of points in ssq if ssq is NULL
##
## nsbsq: positive integer; number of points in sbsq if sbsq is NULL
##
## npip: positive integer; number of points in pip if pip is NULL
##
## ssq_mult: positive numeric; used to obtain an upper bound for ssq
##
## ssq_lower: positive numeric; if ssq is NULL, then ssq_lower will
## be the least value in ssq
##
## snr_upper: positive numeric; used to obtain an upper bound for sbsq
##
## sbsq_lower: positive numeric; if sbsq is NULL, then sbsq_lower will be the
## least value in sbsq
##
## pip_lower: numeric in (0, 1); if pip is NULL, then pip_lower will be the
## least value in pip
##
## pip_upper: NULL OR  numeric in(0, 1); if pip is NULL, then pip_upper will be
## the greatest value in pip
##
## alpha_tol: positive numeric; end CAVI when the Frobenius norm of the
## change in the alpha matrix is within alpha_tol
##
## max_iter: positive integer; if a tolerance criteria has not been met by
## max_iter_grid iterations, end CAVI
## -----------------------------------------------------------------------------
cavi <- function(X, Z, D, y, hp_method, ssq, sbsq, pip, nssq, nsbsq, npip,
                 ssq_mult, ssq_lower, snr_upper, sbsq_lower, pip_lower,
                 pip_upper, alpha_tol, max_iter){

  # get the dimensions of the data
  n <- nrow(X)
  p <- ncol(X)

  # instantiate initial values of variational parameters; the l, j entry is
  # the variational approximation to the j-th parameter in a regression with
  # the resp_index predictor fixed as the response with weightings taken with
  # respect to the l-th individual
  alpha <- matrix(0.2, n, p)
  mu <- matrix(0, n, p)

  # if the hyperparameter grid has not been fully supplied, create it
  if (is.null(ssq) | is.null(sbsq) | is.null(pip)){

    # if either sbsq or pip has not been supplied, then use LASSO to estimate
    # the proportion of non-zero coefficients
    if ((is.null(sbsq) | is.null(pip)) & is.null(pip_upper)){
      lasso <- glmnet::cv.glmnet(X, y)

      # find the number of non-zero coefficients estimated by LASSO
      # ensure that non0 is an integer in [1, p - 1]
      non0 <- sum(coef(lasso, s = "lambda.1se")[-1] != 0)
      non0 <- min(max(non0, 1), p - 1)

      # an upper bound for pi is the proportion of non-zero coefficients
      pip_upper <- non0 / p
    }

    # if ssq has not been supplied, create the grid
    if (is.null(ssq)){

      # find the upper bound for the grid of ssq
      ssq_upper <- ssq_mult * var(y)

      # create the grid candidates for ssq
      ssq <- seq(ssq_lower, ssq_upper, length.out = nssq)

    }

    # if sbsq has not been supplied, create the grid
    if (is.null(sbsq)){

      # find the sum of the variances for each of the columns of X_mat
      s2_sum <- sum(apply(X, 2, var))

      # find the upper bound for the grid of sbsq
      sbsq_upper <- snr_upper / (pip_upper * s2_sum)

      # create the grid candidates for sbsq
      sbsq <- seq(sbsq_lower, sbsq_upper, length.out = nsbsq)

    }

    # if pip has not been supplied, create the grid
    if (is.null(pip)){

      # create posterior inclusion probability grid
      pip <- seq(pip_lower, pip_upper, length.out = npip)

    }
  }

  # create the grid
  hp <- expand.grid(pip = pip, ssq = ssq, sbsq = sbsq, elbo = NA, iter = NA)

  if (hp_method == "grid_search"){

    # perform grid search to select the best hyperparameters

    # perform CAVI for each of the hyperparameter settings
    out_grid <- grid_search_c(y, D, X, mu, alpha, hp$ssq, hp$sbsq, hp$pip,
                              alpha_tol, max_iter)

    # add the ELBO and converged iter for each of the grid points to hp grid
    hp$elbo <- out_grid$elbo_vec
    hp$iter <- out_grid$conv_iter

    # use the best hyperparameters from the grid search to perform CAVI until
    # max_iter is reached or alpha converges
    out <- cavi_c(y, D, X, out_grid$mu, out_grid$alpha, out_grid$ssq,
                  out_grid$sbsq, out_grid$pip, alpha_tol, max_iter)

    # save hyperparameter details
    final <- c(ssq = out_grid$ssq, sbsq = out_grid$sbsq, pip = out_grid$pip)
    hp <- list(grid = hp, final = final)

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
                                  alpha_tol, max_iter)

        # use the best hyperparameters from the grid search to perform CAVI until
        # max_iter is reached or alpha converges
        out <- cavi_c(y, D, X, out_grid$mu, out_grid$alpha, out_grid$ssq,
                      out_grid$sbsq, pip[j], alpha_tol, max_iter)

        # save the final hyperparameters and elbo
        hyp[j, c("ssq", "sbsq", "elbo")] <- c(out_grid$ssq, out_grid$sbsq, out$elbo)

        # fix the final hyperparameter setting
        hp_j <- data.frame(ssq = out_grid$ssq, sbsq = out_grid$sbsq, pip = pip[j])

      }else{

        # otherwise, model averaging CAVI

        # fix hyperparameter setting
        hp_j <- hp[j, ]

        # perform CAVI for the hyperparameter setting
        out <- cavi_c(y, D, X, mu, alpha, hp_j$ssq, hp_j$sbsq, hp_j$pip,
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

      # save the number of iterations to converge
      hp$iter[j] <- out$conv_iter
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
    hp <- list(grid = hp, final = "<NA> if hp_method == 'model_average'")

    # if hybrid, add the final hyperparameters chosen by grid search to
    # the hyperparameter details
    if(hp_method == "hybrid") hp$final <- hyp
  }

  # return the final variational parameters, hyperparameter details, and elbo
  return(list(alpha_matrix = alpha, mu_matrix = mu, ssqv_matrix = ssq_var,
              hyperparameters = hp, elbo = elbo))
}
