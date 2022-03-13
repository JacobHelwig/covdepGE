## -----------------------------------------------------------------------------
## -----------------------------ELBO_calculator_R-------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## Calculates ELBO for a fixed response j and individual l
## -----------------------------ARGUMENTS---------------------------------------
## y: n x 1 vector; responses (j-th column of the data)
## D: n x 1 vector; weights (i-th entry is the weight of the i-th
## individual with respect to the l-th individual using the l-th individual's
## bandwidth)
## X_mat: n x p matrix; data_mat with the j-th column removed
## S_sq, mu, alpha: p x 1 vectors; variational parameters. the k-th entry is the
## k-th parameter for the l-th individual
## sigmasq, sigmabeta_sq: double; spipke and slab variance hyperparameters
## pip_est: double; spipke and slab probability of inclusion
## -----------------------------RETURNS-----------------------------------------
## Returns: double; ELBO for the l-th individual and j-th column fixed as the
## response
## -----------------------------------------------------------------------------
ELBO_calculator_R <- function(y, D, X, ssq_var, mu, alpha, ssq, sbsq, pip) {

  # square of: y minus X matrix-multiplied with element-wise product of mu and
  # alpha
  y_Xmualpha_sq <- (y - X %*% (mu * alpha))^2

  # square of mu vector
  mu_sq <- mu^2

  # calculate each of the terms that sum to ELBO
  t1 <- -sum(D * y_Xmualpha_sq) / (2 * ssq)
  t2 <- -sum(D * ((X^2) %*% (alpha * (mu_sq + ssq_var) - alpha^2 * mu_sq))) / (2 * ssq)
  t3 <- sum(alpha * (1 + log(ssq_var))) / 2
  t4 <- -sum(alpha * log((alpha + 0.000001) / pip) + (1 - alpha) * log(
    (1 - alpha + 0.000001) / (1 - pip)))
  t5 <- -sum(alpha * ((mu_sq + ssq_var) / (2 * ssq * sbsq) + log(ssq * sbsq) / 2))
  t6 <- 0.5 * log(1 / (2 * + pi * ssq))

  return(t1 + t2 + t3 + t4 + t5 + t6)
}

## -----------------------------------------------------------------------------
## -----------------------------total_elbo_R------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## for a fixed response, calculate and sum the ELBO for all individuals in the
## data
## -----------------------------ARGUMENTS---------------------------------------
## y: n x 1 vector; responses (j-th column of the data)
## D: n x n matrix; weights (k,l entry is the weight of the k-th individual
## with respect to the l-th individual using the l-th individual's bandwidth)
## X_mat: n x p matrix; data_mat with the j-th column removed
## S_sq, mu, alpha: n x p matrices; variational parameters. the l, k
## entry is the k-th parameter for the l-th individual
## sigmasq, sigmabeta_sq: doubles; spipke and slab variance hyperparameters
## pip_est: double; spipke and slab probability of inclusion
## -----------------------------RETURNS-----------------------------------------
## elbo_tot: double; ELBO for the l-th individual with j-th column fixed as the
## response
## -----------------------------------------------------------------------------
total_ELBO_R <- function(y, D, X, ssq_var, mu, alpha, ssq, sbsq, pip){

  # get sample size
  n <- nrow(X)

  # instantiate variable to store the total ELBO
  elbo_tot <- 0

  # loop over all individuals
  for (l in 1:n){

    ## calculate the ELBO for l-th individual and add it to the total ELBO
    elbo_tot <- elbo_tot + ELBO_calculator_R(y, D[ , l], X, ssq_var[l, ],
                                             mu[l, ], alpha[l, ], ssq, sbsq, pip)
    # elbo_tot <- elbo_tot + ELBO_calculator_c(
    #   y, D[ , l], X_mat, S_sq[l, ], mu_mat[l, ], alpha_mat[l, ],
    #   sigmasq[l], sigmabeta_sq[l], pip_est)
  }

  return(elbo_tot)

}

## -----------------------------------------------------------------------------
## -----------------------------mu_update_R-------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## for a fixed response, update mu matrix
## -----------------------------ARGUMENTS---------------------------------------
## y: n x 1 vector; responses (j-th column of the data)
## D: n x n matrix; weights (k,l entry is the weight of the k-th individual
## with respect to the l-th individual using the l-th individual's bandwidth)
## X_mat: n x p matrix; data_mat with the j-th column removed
## S_sq, mu, alpha: n x p matrices; variational parameters. the l, k entry is
## the k-th parameter for the l-th individual
## sigmasq: double; spipke and slab variance hyperparameter
## -----------------------------------------------------------------------------
mu_update_R <- function(y, D, X, ssq_var, mu, alpha, ssq){

  # get sample size and number of parameters
  n <- nrow(X)
  p <- ncol(X)

  # instantiate matrices for the update loop
  mu_stack <- matrix(NA, n, p)
  alpha_stack <- matrix(NA, n, p)
  X_mu_alpha <- matrix(NA, n, p)
  X_mu_alpha_k <- matrix(NA, n, p)
  y_k <- matrix(NA, n, p)
  d_x_y <- matrix(NA, n, p)

  # loop over the individuals to update mu row by row
  for (l in 1:n){

    # l-th row of mu_mat, alpha_mat stacked n times
    mu_stack <- matrix(mu[l, ], n, p, T)
    alpha_stack <- matrix(alpha[l, ], n ,p, T)

    # take the element-wise product of X_mat, mu_stack, and alpha stack
    X_mu_alpha <- X * mu_stack * alpha_stack

    # the k-th column is the rowSums of X_mu_alpha minus the k-th column of
    # X_mu_alpha (accounts for m \neq k in summation)
    X_mu_alpha_k <- matrix(rowSums(X_mu_alpha), n, p) - X_mu_alpha

    # the k-th column is y minus the k-th column of X_mu_alpha_k
    y_k <- matrix(y, n, p) - X_mu_alpha_k

    # the k-th column is d_:,l * x_:,k * y_k_:,k
    d_x_y <- matrix(D[ , l], n, p) * X * y_k

    # the update of the l-th row of mu
    mu[l, ] <- ssq_var[l, ] * colSums(d_x_y) / ssq
  }

  return(mu)
}

## -----------------------------------------------------------------------------
## -----------------------------alpha_update_R----------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## for a fixed response, update alpha matrix
## -----------------------------ARGUMENTS---------------------------------------
## mu_mat, alpha_mat: n x p matrices; variational parameters. the l, k entry is
## the k-th parameter for the l-th individual
## alpha_logit_term 1, 2, 3: double (1), n x p matrices (2, 3): terms used to
## calculate the logit of alpha
## upper_limit: double; during the alpha update, values of logit(alpha) greater
## than upper_limit will be assigned a probability of 1; this avoids issues
## with exponentiation of large numbers creating Infinity divided by Infinity
## -----------------------------------------------------------------------------
alpha_update_R <- function(ssq_var, mu, alpha, alpha1, alpha2_denom, alpha3,
                           ssq, sbsq, pip){

  # calculate the logit of alpha
  alpha_logit <- alpha1 + ((mu^2) / alpha2_denom) + alpha3

  # transform from logit to probabilities of inclusion; update alpha_mat
  exp_logit <- exp(alpha_logit)
  alpha <- exp_logit / (1 + exp_logit)

  # handle NA's due to division by infinity resulting from exponentiation of
  # large values; these probabilities are indescernible from 1

  # find large values
  upper_limit <- 9
  index1 <- which(alpha_logit > upper_limit)

  # replace them
  alpha[index1] <- 1

  return(alpha)
}

## -----------------------------------------------------------------------------
## -----------------------------cavi_R------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## The core function that performs CAVI to calculate and return the final
## variational estimates of a single regression for each of the n individuals
## -----------------------------ARGUMENTS---------------------------------------
## y: n x 1 vector; responses (j-th column of the data)
## D: n x n matrix; weights (k,l entry is the weight of the k-th individual
## with respect to the l-th individual using the l-th individual's bandwidth)
## X_mat: n x p matrix; data_mat with the j-th column removed
## mu_mat, alpha_mat: n x p matrices; variational parameters. the l, k entry is
## the k-th parameter for the l-th individual
## sigmasq, sigmabeta_sq: doubles; spipke and slab variance hyperparameters
## pip_est: double; spipke and slab probability of inclusion
## tolerance: double; when the square root of the sum of squared changes in
## the elements of alpha are within tolerance, stop iterating
## max_iter: integer; maximum number of iterations
## upper_limit: double; during the alpha update, values of logit(alpha) greater
## than upper_limit will be assigned a probability of 1
## -----------------------------RETURNS-----------------------------------------
## var_alpha: n x p matrix; final alpha values
## var_ELB: double; final value of ELBO summed across all individuals
## converged_iter: integer; number of iterations to reach convergence
## sigmasq: n x 1 vector; fitted error term variance for each individual
## sigmabeta_sq: n x 1 vector; fitted slab variance for each individual
## -----------------------------------------------------------------------------
## [[Rcpp::export]]
cavi_R <- function(y, D, X, mu_, alpha_, ssq, sbsq, pip, tol, max_iter) {

  n <- nrow(X)
  p <- ncol(X)

  # instantiate matrices for updated variational parameters with starting
  # values dictated by the matrices passed as arguments
  mu <- matrix(mu_, n, p)
  alpha <- matrix(alpha_, n, p)

  # matrices for tracking the convergence of alpha parameters
  alpha_last <- matrix(NA, n, p)
  change_alpha <- matrix(NA, n, p)

  # integer for tracking the iteration at which convergence is reached
  converged_iter <- max_iter

  # variational variance of regression coefficients update
  ssq_var <- ssq / (t(t(X^2) %*% D) + 1 / sbsq)

  # 1st and 3rd term of the alpha update, denominator of the second term
  alpha1 <- log(pip / (1 - pip))
  alpha3 <- log(sqrt(ssq_var / (ssq * sbsq)))
  alpha2_denom <- 2 * ssq_var

  # CAVI loop (optimize variational parameters)
  for (k in 1:max_iter){

    # mu update
    mu <- mu_update_R(y, D, X, ssq_var, mu, alpha, ssq)
    #mu_update_c(y, D, X_mat, S_sq, mu, alpha, sigmasq);

    # alpha update

    # save the last value of alpha and update it
    alpha_last <- matrix(alpha, n, p)
    alpha <- alpha_update_R(ssq_var, mu, alpha, alpha1, alpha2_denom, alpha3,
                            ssq, sbsq, pip)
    #alpha_update_c(S_sq, mu, alpha, sigmasq, sigmabeta_sq, pip_est)

    # calculate change in alpha
    change_alpha <- alpha - alpha_last;

    # if the square root of the sum of squared changes in alpha is within the
    # tolerance, break from the for loop
    if (sqrt(sum((change_alpha^2))) < tol){
      converged_iter <- k
      break;
    }
  }

  # calculate ELBO across n individuals
  ELBO <- total_ELBO_R(y, D, X, ssq_var, mu, alpha, ssq, sbsq, pip)
  #ELBO <- total_ELBO_c(y, D, X_mat, S_sq, mu, alpha, sigmasq, sigmabeta_sq, pip_est)

  return(list(mu = mu, alpha = alpha, elbo = ELBO,
              converged_iter = converged_iter, ssq = ssq, sbsq = sbsq))
}

## -----------------------------------------------------------------------------
## -----------------------------grid_search_R-----------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## for a fixed response, run CAVI for each individual across a grid of
## hyperparameters and return the ELBO for each of the grid points
## -----------------------------ARGUMENTS---------------------------------------
## y: n x 1 vector; responses (j-th column of the data)
## D: n x n matrix; weights (k,l entry is the weight of the k-th individual
## with respect to the l-th individual using the l-th individual's bandwidth)
## X_mat: n x p matrix; data_mat with the j-th column removed
## mu_mat, alpha_mat: n x p; matrices of variational parameters. the l, k entry
## is the k-th parameter for the l-th individual
## sigmasq_vec: n_param x 1 vector; error term variance candidates
## sigmabeta_sq_vec: n_param x 1 vector; slab variance candidates
## pip_vec: n_param x 1 vector; candidate spipke and slab probabilities of inclusion
## tolerance: double; when the square root of the sum of squared changes in
## the elements of alpha are within tolerance, stop iterating
## -----------------------------RETURNS-----------------------------------------
## returns list with 2 values:
## elbo: n_param x 1 vector; ELBO for each hyperparameter candidate
## num_converged: integer; the number of hyperparameter candidates for which
## the CAVI converged
## -----------------------------------------------------------------------------
## [[Rcpp::export]]
grid_search_R <- function(y, D, X, mu, alpha, ssq, sbsq, pip, tol, max_iter){

  # get dimensions of the data
  n <- nrow(X)
  p <- ncol(X)

  # get the number of grid points
  n_param <- length(pip)

  # storage for the best ELBO
  elbo_best <- -Inf

  # storage for the hyperparameters and variational parameters corresponding to
  # the current best ELBO
  mu_best <- matrix(NA, n, p)
  alpha_best <- matrix(NA, n, p)
  ssq_best <- rep(NA, n)
  sbsq_best <- rep(NA, n)
  pip_best <- NA

  # instantiate a list to store the result of cavi_R
  out <- list()

  # perform CAVI for each grid point
  for (j in 1:n_param){

    # run CAVI
    out <- cavi_R(y, D, X, mu, alpha, ssq[j], sbsq[j], pip[j], tol, max_iter)
    # out <- cavi_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq_vec[ , j],
    #               update_sigmasq, sigmabetasq_vec[ , j], update_sigmabetasq,
    #               pip_vec[j], tolerance, max_iter, upper_limit)

    # if the new ELBO is greater than the current best, update the best ELBO and
    # corresponding return values
    if (elbo_best < out$elbo){

      elbo_best <- out$elbo
      mu_best <- out$mu
      alpha_best <- out$alpha
      ssq_best <- ssq[j]
      sbsq_best <- sbsq[j]
      pip_best <- pip[j]
      iter_best <- out$converged_iter
    }
  }

  return(list(elbo = elbo_best, mu = mu_best, alpha = alpha_best, ssq = ssq_best,
              sbsq = sbsq_best, pip = pip_best, iterations = iter_best))
}
