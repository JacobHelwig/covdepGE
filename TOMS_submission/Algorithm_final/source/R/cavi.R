## -----------------------------------------------------------------------------
## -----------------------------cavi--------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## Performs CAVI and hyperparameter selection for n linear regressions, where
## the l-th regression is weighted with respect to the l-th observation
## -----------------------------ARGUMENTS---------------------------------------
## X: n x p numeric matrix; data with the j-th column removed
##
## Z: n x q numeric matrix; extraneous covariates
##
## D: n x n numeric matrix; i, j entry is the weighting of the i-th observation
## with respect to the j-th observation using the j-th observation's bandwidth
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
## max_iter: positive integer; if tolerance criteria has not been met by
## max_iter iterations, end CAVI
##
## max_iter_grid: positive integer; if tolerance criteria has not been met by
## max_iter_grid iterations during grid search, end CAVI
## -----------------------------------------------------------------------------
cavi <- function(X, Z, D, y, hp_method, ssq, sbsq, pip, nssq, nsbsq, npip,
                 ssq_mult, ssq_lower, snr_upper, sbsq_lower, pip_lower,
                 pip_upper, alpha_tol, max_iter, max_iter_grid) {

  # get the dimensions of the data
  n <- nrow(X)
  p <- ncol(X)

  # instantiate initial values of variational parameters
  alpha <- matrix(0.2, n, p)
  mu <- matrix(0, n, p)

  # if the hyperparameter grid has not been fully supplied, create it
  if (is.null(ssq) | is.null(sbsq) | is.null(pip)) {

    # if either sbsq or pip has not been supplied and pip_upper has not been
    # supplied, use LASSO to estimate the proportion of non-zero coefficients
    if ((is.null(sbsq) | is.null(pip)) & is.null(pip_upper)) {

      # fit the lasso
      lasso <- glmnet::cv.glmnet(X, y)

      # find the number of non-zero coefficients estimated by LASSO and ensure
      # that non0 is an integer in [1, p - 1]
      non0 <- sum(glmnet::coef.glmnet(lasso, s = "lambda.1se")[-1] != 0)
      non0 <- min(max(non0, 1), p - 1)

      # calculate the proportion of non-zero coefficients
      pip_upper <- non0 / p
    }

    # if ssq has not been supplied, create the grid
    if (is.null(ssq)) {

      # find the upper bound for the grid of ssq
      ssq_upper <- ssq_mult * stats::var(y)

      # create the grid candidates for ssq
      ssq <- seq(ssq_lower, ssq_upper, length.out = nssq)
    }

    # if sbsq has not been supplied, create the grid
    if (is.null(sbsq)) {

      # find the sum of the variances for each of the columns of X
      s2_sum <- sum(apply(X, 2, stats::var))

      # find the upper bound for the grid of sbsq
      sbsq_upper <- snr_upper / (pip_upper * s2_sum)

      # create the grid candidates for sbsq
      sbsq <- seq(sbsq_lower, sbsq_upper, length.out = nsbsq)
    }

    # if pip has not been supplied, create the grid
    if (is.null(pip)) {

      # create posterior inclusion probability grid
      pip <- seq(pip_lower, pip_upper, length.out = npip)
    }
  }

  # create the grid
  hp <- expand.grid(pip = pip, ssq = ssq, sbsq = sbsq, elbo = NA, iter = NA)
  hp <- unique(hp)

  if (hp_method == "grid_search") {

    # perform grid search to select the best hyperparameters

    # perform CAVI for each of the hyperparameter settings
    out_grid <- grid_search_c(
      y, D, X, mu, alpha, hp$ssq, hp$sbsq, hp$pip,
      alpha_tol, max_iter_grid
    )

    # add the elbo and the converged iter to hp
    hp[, c("elbo", "iter")] <- data.frame(out_grid[c("elbo_vec", "iter")])

    # use the best hyperparameters from the grid search to perform CAVI again
    out <- cavi_c(
      y, D, X, out_grid$mu, out_grid$alpha, out_grid$ssq,
      out_grid$sbsq, out_grid$pip, alpha_tol, max_iter
    )

    # save hyperparameter details
    final <- data.frame(
      ssq = out_grid$ssq, sbsq = out_grid$sbsq,
      pip = out_grid$pip, elbo = out$elbo, iter = out$iter
    )
    hp <- list(grid = hp, final = final)

    # save final variational parameters, elbo, and iterations to converge
    alpha <- out$alpha
    mu <- out$mu
    ssq_var <- out$ssq_var
    elbo <- out$elbo
  } else {

    # otherwise, hp_method is model_average or hybrid

    # create a data.frame for storing final hyperparameters from grid search
    if (hp_method == "hybrid") {
      hyp <- data.frame(ssq = NA, sbsq = NA, pip = pip, elbo = NA, iter = NA)
    }

    # get number of parameters to average over; if hybrid, only averaging over
    # pip; otherwise, average over entire grid
    n_hp <- ifelse(hp_method == "hybrid", length(pip), nrow(hp))

    # create lists for storing the variational parameters and elbos for each
    # hyperparameter setting
    elbo_theta <- alpha_theta <- mu_theta <- ssqv_theta <- vector("list", n_hp)

    # iterate over each of the pip if hybrid
    # iterate over each of the hyperparameter settings otherwise
    for (j in 1:n_hp) {

      # hybrid CAVI
      if (hp_method == "hybrid") {

        # fix the values in hp corresponding to the j-th value of pip
        hp_j <- hp[hp$pip == pip[j], ]

        # perform CAVI for each of the hyperparameter settings
        out_grid <- grid_search_c(
          y, D, X, mu, alpha, hp_j$ssq, hp_j$sbsq,
          hp_j$pip, alpha_tol, max_iter_grid
        )

        # add the elbo and the converged iter to hp
        hp[hp$pip == pip[j], c("elbo", "iter")] <- data.frame(
          out_grid[c("elbo_vec", "iter")]
        )

        # use the best hyperparameters from grid search to perform CAVI again
        out <- cavi_c(
          y, D, X, out_grid$mu, out_grid$alpha, out_grid$ssq,
          out_grid$sbsq, pip[j], alpha_tol, max_iter
        )

        # save the final hyperparameters, elbo, and iterations to converge
        hyp[j, ] <- c(
          out_grid$ssq, out_grid$sbsq, out_grid$pip, out$elbo,
          out$iter
        )

        # fix the final hyperparameter setting
        hp_j <- data.frame(
          ssq = out_grid$ssq, sbsq = out_grid$sbsq,
          pip = pip[j]
        )
      } else {

        # otherwise, model averaging CAVI

        # fix hyperparameter setting
        hp_j <- hp[j, ]

        # perform CAVI for the hyperparameter setting
        out <- cavi_c(
          y, D, X, mu, alpha, hp_j$ssq, hp_j$sbsq, hp_j$pip,
          alpha_tol, max_iter
        )

        # add the elbo and the converged iter to hp
        hp[j, c("elbo", "iter")] <- unlist(out[c("elbo", "iter")])
      }

      # save the variational parameters and elbo
      alpha_theta[[j]] <- out$alpha
      mu_theta[[j]] <- out$mu
      ssqv_theta[[j]] <- out$ssq_var

      # calculate the elbo for each observation under the current hyperparameter
      # setting and save
      elbo_l <- rep(NA, n)
      for (l in 1:n) {
        elbo_l[l] <- ELBO_calculator_c(
          y, D[, l], X, t(out$ssq_var[l, ]),
          t(out$mu[l, ]), t(out$alpha[l, ]),
          hp_j$ssq, hp_j$sbsq, hp_j$pip
        )
      }


      # save the ELBO for all observations
      elbo_theta[[j]] <- elbo_l
    }

    # calculate weights for averaging and average

    # vector for saving the average elbo for each observation
    elbo_avg <- rep(NA, n)

    # calculate the weights for each observation
    for (l in 1:n) {

      # retrieve the elbo for the l-th observation for each hyperparameter
      # setting
      elbo_l <- sapply(elbo_theta, `[`, l)

      # subtract the maximum value from each and exponentiate to get the lower
      # bound for the likelihood (subtracting max for numerical stability)
      lik_l <- exp(elbo_l - max(elbo_l))

      # normalize the weights
      weight <- lik_l / sum(lik_l)

      # calculate the average elbo for this observation
      elbo_avg[l] <- sum(elbo_l * weight)

      # multiply the l-th row of each of the variational parameter matrices by
      # the corresponding weight
      for (k in 1:n_hp) {
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
    hp <- list(grid = hp, final = NA)

    # if hybrid, add the final hyperparameters chosen by grid search to
    # the hyperparameter details
    if (hp_method == "hybrid") hp$final <- hyp
  }

  # return the final variational parameters, hyperparameter details, and elbo
  return(list(
    alpha_matrix = alpha, mu_matrix = mu, ssqv_matrix = ssq_var,
    hyperparameters = hp, elbo = elbo
  ))
}
