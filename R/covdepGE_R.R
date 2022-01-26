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
## sigmasq, sigmabeta_sq: double; spike and slab variance hyperparameters
## pi_est: double; spike and slab probability of inclusion
## -----------------------------RETURNS-----------------------------------------
## Returns: double; ELBO for the l-th individual and j-th column fixed as the
## response
## -----------------------------------------------------------------------------
ELBO_calculator_R <- function(y, D, X_mat, S_sq, mu, alpha, sigmasq,
                              sigmabeta_sq, pi_est) {

  # square of: y minus X matrix-multiplied with element-wise product of mu and
  # alpha
  y_Xmualpha_sq <- (y - X_mat %*% (mu * alpha))^2

  # square of mu vector
  mu_sq <- mu^2

  # calculate each of the terms that sum to ELBO
  t1 <- -sum(D * y_Xmualpha_sq) / (2 * sigmasq)
  t2 <- -sum(D * ((X_mat^2) %*% (alpha * (mu_sq + S_sq) - alpha^2 * mu_sq))) / (
    2 * sigmasq)
  t3 <- sum(alpha * (1 + log(S_sq))) / 2
  t4 <- -sum(alpha * log((alpha + 0.000001) / pi_est) + (1 - alpha) * log(
   (1 - alpha + 0.000001) / (1 - pi_est)))
  t5 <- -sum(alpha * ((mu_sq + S_sq) / (2 * sigmasq * sigmabeta_sq) + log(
   sigmasq * sigmabeta_sq) / 2))
  t6 <- 0.5 * log(1 / (2 * + pi * sigmasq))

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
## sigmasq, sigmabeta_sq: doubles; spike and slab variance hyperparameters
## pi_est: double; spike and slab probability of inclusion
## -----------------------------RETURNS-----------------------------------------
## elbo_tot: double; ELBO for the l-th individual with j-th column fixed as the
## response
## -----------------------------------------------------------------------------
total_ELBO_R <- function(y, D, X_mat, S_sq, mu_mat, alpha_mat, sigmasq,
                         sigmabeta_sq, pi_est){

  # get sample size
  n <- nrow(X_mat)

  # instantiate variable to store the total ELBO
  elbo_tot <- 0

  # loop over all individuals
  for (l in 1:n){

   ## calculate the ELBO for l-th individual and add it to the total ELBO
   elbo_tot <- elbo_tot + ELBO_calculator_R(
     y, D[ , l], X_mat, S_sq[l, ], mu_mat[l, ], alpha_mat[l, ],
     sigmasq[l], sigmabeta_sq[l], pi_est)
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
## sigmasq: double; spike and slab variance hyperparameter
## -----------------------------------------------------------------------------
mu_update_R <- function(y, D, X_mat, S_sq, mu, alpha, sigmasq){

 # get sample size and number of parameters
 n <- nrow(X_mat)
 p <- ncol(X_mat)

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
  X_mu_alpha <- X_mat * mu_stack * alpha_stack

  # the k-th column is the rowSums of X_mu_alpha minus the k-th column of
  # X_mu_alpha (accounts for m \neq k in summation)
  X_mu_alpha_k <- matrix(rowSums(X_mu_alpha), n, p) - X_mu_alpha

  # the k-th column is y minus the k-th column of X_mu_alpha_k
  y_k <- matrix(y, n, p) - X_mu_alpha_k

  # the k-th column is d_:,l * x_:,k * y_k_:,k
  d_x_y <- matrix(D[ , l], n, p) * X_mat * y_k

  # the update of the l-th row of mu
  mu[l, ] <- S_sq[l, ] * colSums(d_x_y) / sigmasq[l]
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
alpha_update_R <- function(S_sq, mu, alpha, sigmasq, sigmabeta_sq, pi_est,
                           upper_limit = 9){

  # get dimensions of the data
  n <- nrow(S_sq)
  p <- ncol(S_sq)

  # 1st and 3rd term of the alpha update, denominator of the second term
  alpha_logit_term1 <- log(pi_est / (1 - pi_est))
  alpha_logit_term3 <- log(sqrt(S_sq / matrix(sigmasq * sigmabeta_sq, n, p)))
  alpha_logit_term2_denom <- 2 * S_sq

  # calculate the logit of alpha
  alpha_logit <- (alpha_logit_term1 + ((mu^2) / alpha_logit_term2_denom) +
                    alpha_logit_term3)

  # transform from logit to probabilities of inclusion; update alpha_mat
  exp_logit <- exp(alpha_logit)
  alpha <- exp_logit / (1 + exp_logit)

  # handle NA's due to division by infinity resulting from exponentiation of
  # large values; these probabilities are indescernible from 1

  # find large values
  index1 <- which(alpha_logit > upper_limit)

  # replace them
  alpha[index1] <- 1

  return(alpha)
}

## -----------------------------------------------------------------------------
## -----------------------------sigmasq_update_R--------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## Update the residual variance term using MAPE
## -----------------------------ARGUMENTS---------------------------------------
## y: n x 1 vector; responses (j-th column of the data)
## D: n x n matrix; weights (k,l entry is the weight of the k-th individual
## with respect to the l-th individual using the l-th individual's bandwidth)
## X_mat: n x p matrix; data_mat with the j-th column removed
## mu_mat, alpha_mat, S_sq: n x p matrices; variational parameters. the l, k
## entry is the k-th parameter for the l-th individual
## sigmasq, sigmabeta_sq: doubles; spike and slab variance hyperparameters
## -----------------------------RETURNS-----------------------------------------
## sigmasq: n x 1 vector; updated residual variance
## -----------------------------------------------------------------------------
## [[Rcpp::export]]
sigmasq_update_R <- function(y, D, X_mat, S_sq, mu_mat, alpha_mat, sigmasq,
                             sigmabeta_sq) {

  # get dimensions of the data
  n <- nrow(X_mat)
  p <- ncol(X_mat)

  # calulate alpha^2 and mu^2
  alpha_sq <- alpha_mat^2
  mu_sq <- mu_mat^2

  # calculate expected value of beta for each individual; l-th row is
  # E(beta) for individual l
  rho <- mu_mat * alpha_mat

  # find fitted values using expected value of beta for each individual; l-th
  # column is fitted values for individual l
  fitted <- X_mat %*% t(rho)

  # calculate the squared residuals for each of the fitted values for each
  # individual; l-th column is residuals for individual l
  resid2 <- (matrix(y, n, n) - fitted)^2

  # calculate the sum of the weighted squared residuals for each individual;
  # l-th value is the SWSR for individual l
  resid_w <- colSums(resid2 * D)

  # calculate the second values in the numerator for each individual; the l-th
  # row is for individual l
  num_term2 <- alpha_mat * S_sq + alpha_mat * mu_sq - alpha_sq * mu_sq

  # calculate the third values in the numerator for each individual; the l-th
  # value is for individual l
  num_term3 <- rowSums(alpha_mat * (S_sq + mu_sq)) / sigmabeta_sq

  # calculate the denominator for each individual; l-th value is for individual l
  denom <- rowSums(alpha_mat) + n

  # iterate over the individuals to update each error variance
  for (l in 1:n){

    # calculate weighted versions of y and X
    weights <- sqrt(D[ , l])
    X_w <- X_mat * matrix(weights, n, p)

    # diagonal elements of X transpose X weighted
    XtX_w <- diag(t(X_w) %*% X_w)

    # second numerator term
    num_term2_l <- t(XtX_w) %*% num_term2[l , ]

    # apply update
    sigmasq[l] <- (resid_w[l] + num_term2_l + num_term3[l]) / denom[l]
  }

  return(sigmasq)
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
## sigmasq, sigmabeta_sq: doubles; spike and slab variance hyperparameters
## pi_est: double; spike and slab probability of inclusion
## tolerance: double; when the square root of the sum of squared changes in
## the elements of alpha are within tolerance, stop iterating
## max_iter: integer; maximum number of iterations
## upper_limit: double; during the alpha update, values of logit(alpha) greater
## than upper_limit will be assigned a probability of 1
## -----------------------------RETURNS-----------------------------------------
## var_alpha: n x p matrix; final alpha values
## var_ELB: double; final value of ELBO summed across all individuals
## converged_iter: integer; number of iterations to reach convergence
## -----------------------------------------------------------------------------
## [[Rcpp::export]]
cavi_R <- function(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq,
                   pi_est, tolerance, max_iter, upper_limit = 9) {

  n <- nrow(X_mat)
  p <- ncol(X_mat)

  # instantiate initial values of variance hyperparameters for all individuals
  sigmasq <- rep(sigmasq, n)
  sigmabeta_sq <- rep(sigmabeta_sq, n)

  # instantiate matrices for updated variational parameters with starting
  # values dictated by the matrices passed as arguments
  mu <- mu_mat
  alpha <- alpha_mat

  # matrices for tracking the convergence of alpha parameters
  alpha_last <- matrix(NA, n, p)
  change_alpha <- matrix(NA, n, p)

  # integer for tracking the iteration at which convergence is reached
  converged_iter <- max_iter

  # CAVI loop (optimize variational parameters)
  for (k in 1:max_iter){

    # S_sq update
    S_sq <- matrix(sigmasq, n, p) / (t(t(X_mat^2) %*% D) + 1 /
                                       matrix(sigmabeta_sq, n, p))

    # mu update
    mu <- mu_update_R(y, D, X_mat, S_sq, mu, alpha, sigmasq);

    # alpha update

    # save the last value of alpha and update it
    alpha_last <- alpha
    alpha <- alpha_update_R(S_sq, mu, alpha, sigmasq, sigmabeta_sq, pi_est)

    # calculate change in alpha
    change_alpha <- alpha - alpha_last;

    # if the square root of the sum of squared changes in alpha is within the
    # tolerance, break from the for loop
    if (sqrt(sum((change_alpha^2))) < tolerance){
      converged_iter <- k
      break;
    }

    # update the error term variance using MAPE
    sigmasq <- sigmasq_update_R(y, D, X_mat, S_sq, mu, alpha, sigmasq, sigmabeta_sq)
 }

 # calculate ELBO across n individuals
 ELBO <- total_ELBO_R(y, D, X_mat, S_sq, mu, alpha, sigmasq, sigmabeta_sq, pi_est)

 # return final alpha matrix, the final ELBO, the number of iterations to
 # converge, and the elbo history matrix
 return(list(var_alpha = alpha, var_elbo = ELBO, converged_iter = converged_iter))
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
## pi_vec: n_param x 1 vector; candidate spike and slab probabilities of inclusion
## tolerance: double; when the square root of the sum of squared changes in
## the elements of alpha are within tolerance, stop iterating
## -----------------------------RETURNS-----------------------------------------
## returns list with 2 values:
## elbo: n_param x 1 vector; ELBO for each hyperparameter candidate
## num_converged: integer; the number of hyperparameter candidates for which
## the CAVI converged
## -----------------------------------------------------------------------------
## [[Rcpp::export]]
grid_search_R <- function(y, D, X_mat, mu_mat, alpha_mat, sigmasq_vec,
                          sigmabeta_sq_vec, pi_vec, tolerance, max_iter,
                          upper_limit = 9){

 # get the number of grid points
 n_param <- length(sigmabeta_sq_vec)

 # vector for storing the ELBO corresponding to each grid point
 elbo <- rep(0, n_param)

 # instantiate a list to store the result of cavi_R
 out <- list()

 # double to store the ELBO of the estimated graphs
 elbo_graph <- 0

 # integer to store the number of iterations it took for cavi_R to converge
 converged_iter <- 0

 # count the number of grid points for which convergence was attained
 converged_ct <- 0

 # perform CAVI for each grid point
 for (j in 1:n_param){

  # run CAVI
  out <- cavi_R(y, D, X_mat, mu_mat, alpha_mat, sigmasq_vec[j],
                sigmabeta_sq_vec[j], pi_vec[j], tolerance, max_iter, upper_limit)

  # if cavi_R has converged, increment the convergence count
  converged_iter <- out$converged_iter
  if (converged_iter != max_iter){
    converged_ct <- converged_ct + 1
  }

  # add the ELBO to the vector of ELBOs
  elbo_graph <- out$var_elbo
  elbo[j] <- elbo_graph
 }


 return(list(elbo = elbo, num_converged = converged_ct))
}
