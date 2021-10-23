set.seed(1)
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/Development")
rm(list = ls())
start_time <- Sys.time()


## This function returns the ELBO for a specific variational parameter (corresponding to a fixed individual in study)
## y: response
## X_mat: data matrix except the response variable
## S_sq: variance parameter for the individual
## mu: mean parameter for that individual
## sigmasq and sigmabeta_sq: hyperparameter settings
## pi_est: hyperparameter of spike and slab mixture proportion
## W: the weights associated with the n individuals wrt the covariate of the fixed individual we are studying.
## n: no. of samples
## p: no. of variables in predictor (no. of variables - 1).
ELBO_calculator <- function(y, X_mat, S_sq, mu, alpha, sigmasq, sigmabeta_sq, pi_est, W, n, p) {
  mu <- matrix(mu, p, 1)
  alpha <- matrix(alpha, p, 1)
  s <- matrix(S_sq, p, 1)
  mu_alpha <- matrix(mu * alpha, p, 1)
  W <- matrix(W, n, 1)
  t1 <- -sum(W * (y - X_mat %*% mu_alpha)^2) / (2 * sigmasq)
  t2 <- -sum(W * ((X_mat)^2 %*% (alpha * (mu^2 + s) - alpha^2 * mu^2))) / (2 * sigmasq)
  t3 <- sum(alpha * ((1 + log(s)))) / 2
  t4 <- -sum(alpha * log((alpha + 0.000001) / pi_est) + (1 - alpha) * log((1 - alpha + 0.000001) / (1 - pi_est)))
  t5 <- -sum(alpha * ((mu^2 + s) / (2 * sigmasq * sigmabeta_sq) + log(sigmasq * sigmabeta_sq) / 2))
  t6 <- sum(0.5 * log(1 / (2 * pi * sigmasq)))
  t <- t1 + t2 + t3 + t4 + t5 + t6
  return(t)
}


## The core function that calculates the variational parameter updates and
## returns the final variational estimates for a single regression.
## Arguments:
## y: n by 1 response, i.e the j th variable whose conditional distribution given the
## remaining variables is being calculated.
## X_mat: n by p data matrix with the j-th variable removed
## sigmasq: variance of response given the parameters (homoscedastic part,
## actual variance sigma_sq/w_i)
## sigmabeta_sq: prior variance of coefficient parameter
## pi_est: estimate of spike and slab mixture proportion.
## S_sq, mu_mat, alpha_mat: n by p matrices of variational parameters; the
## i,j-th entry corresponds to the j-th paramter for the i-th individual
cov_vsvb <- function(y, Z, X_mat, sigmasq, sigmabeta_sq, S_sq, pi_est, mu_mat,
                     alpha_mat) {

  n <- nrow(X_mat); p <- ncol(X_mat)

  # threshold for calculating the reverse logit of alpha
  upper_limit <- 9

  # exit condition tolerance
  tol <- 1e-9

  change_alpha <- matrix(Inf, n, p)#rep(0.001, n * p) # alpha_new - alpha_int

  max_iter <- 100
  iter <- 1

  iter <- 1

  # S_sq update
  S_sq <- t(sigmasq * (t(X_mat^2) %*% D + 1 / sigmabeta_sq)^(-1))
  S_sq_vec <- matrix(t(S_sq), n * p, 1)

  # 1st and 3rd term of the alpha update, denominator of the second term
  alpha_logit_term1 <- logit(pi_est)
  alpha_logit_term3 <- log(sqrt(S_sq / (sigmasq * sigmabeta_sq)))
  alpha_logit_term2_denom <- (2 * S_sq)

  # loop to optimize variational parameters
  while (sqrt(sum(change_alpha^2)) > tol & iter < max_iter) {

    # mu update

    for (l in 1:n){

      # the l-th row of mu_mat, alpha_mat stacked n times
      mu_stack <- matrix(mu_mat[l, ], n, p, T)
      alpha_stack <- matrix(alpha_mat[l, ], n, p, T)

      # the element-wise product of X_mat, mu_stack, and alpha stack;
      # the i,j entry is x_i,j * mu_l,j * alpha_l,j
      X_mu_alpha <- X_mat * mu_stack * alpha_stack

      # the k-th column is the rowSums of X_mu_alpha minus the k-th column of
      # X_mu_alpha (accounts for m \neq k in summation)
      X_mu_alpha_k <- matrix(rowSums(X_mu_alpha), n, p) - X_mu_alpha

      # the k-th column is y minus the k-th column of X_mu_alpha_k
      y_k <- matrix(y, n, p) - X_mu_alpha_k

      # the k-th column is d_:,l * x_:,k * y_k_:,k
      d_x_y <- D[ , l] * X_mat * y_k

      # the update of the l-th row of mu
      mu_mat[l, ] <- S_sq[l, ] / sigmasq * colSums(d_x_y)

    }

    # alpha update

    # save the last value of alpha
    alpha_last <- alpha_mat

    # calculate the logit of alpha
    alpha_logit <- (alpha_logit_term1 +
                      (mu_mat^2 / alpha_logit_term2_denom) +
                      alpha_logit_term3)

    # transform from logit to probabilities of inclusion; update alpha_mat
    exp_logit <- exp(alpha_logit)
    alpha_mat <- exp_logit / (1 + exp_logit)

    # handle NA's due to division by infinity resulting from exponentiation of
    # large values; these probabilities are indescernible from 1
    #alpha_mat[is.infinite(exp_logit)] <- 1
    alpha_mat[alpha_logit > upper_limit] <- 1

    # calculate change in alpha
    change_alpha <- alpha_mat - alpha_last

    # calculate ELBO across n individuals
    e <- 0

    for (i in 1:n) { ## calculates ELBO for the j th variable by adding the contribution of the parameter
      ## corresponding to every individual in study. i th iteration takes the contribution of the variational
      ## parameters corresponding to  the i th individual in study, but the information is borrowed from
      ## all the n individuals depending on the weights coded in D[,i]
      e <- e + ELBO_calculator(y, X_mat, S_sq[i, ], mu_mat[i, ], alpha_mat[i, ], sigmasq, sigmabeta_sq, pi_est, D[, i], n, p)
    }

    # want to maximize this by optimizing sigma beta
    ELBO_LB <- e
    iter <- iter + 1
  }

  # return n times p - 1 of each of the variational parameters
  list(var.alpha = alpha_mat, var.mu = mu_mat, var.S_sq = S_sq, var.elbo = ELBO_LB)
}

source("generate_data.R")
library(MASS)
library(varbvs)

logit <- function(x) {
  if ((x == 0) | (x == 1)) {
    return(0)
  } else {
    return(log(x / (1 - x)))
  }
}

# generate data and covariates
discrete_data <- F # true if discrete example is desired
if (discrete_data) {
  dat <- generate_discrete()
  n <- 100
  p <- 10
  tau <- 0.1 # the bandwidth parameter
}else{
  dat <- generate_continuous()
  n <- 180
  p <- 4
  tau <- 0.56
}

data_mat <- dat$data
Z <- dat$covts

# D is an n by n matrix of weights; the i, j entry is the similarity between
# individuals i and j
D <- matrix(1, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    D[j, i] <- dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, tau)
  }
}

# Scale weights to sum to n
for (i in 1:n) {
  D[, i] <- n * (D[, i] / sum(D[, i]))
}

# List for the variable-specific inclusion probability matrix; the i-th element
# in the list is a n by p matrix corresponding to the i-th predictor;
# in this matrix, the j-th row corresponds to the dependence structure for the
# j-th subject, with the i-th predictor fixed as the response
graph_list <- vector("list", p + 1)

# main loop over the predictors
for (resp_index in 1:(p + 1)) {

  # Set variable number `resp_index` as the response
  y <- data_mat[, resp_index]

  # Set the remaining p variables as predictor
  X_mat <- data_mat[, -resp_index]

  # instantiate initial values of variational parameters
  alpha_mat <- matrix(0.2, n, p)
  sigmabeta_sq <- 3
  sigmasq <- 1
  E <- rnorm(n, 0, sigmasq) # removing this causes discrepency in discrete case
  # dont delete S_sq <- matrix(sigmasq * (DXtX + 1 / sigmabeta_sq)^(-1), n, p) # should be byrow = T?
  S_sq <- sigmasq * (t(X_mat^2) + 1 / sigmabeta_sq)^(-1)
  mu_mat <- matrix(0, n, p, byrow = TRUE)

  # Setting hyperparameter values for sigmasq and the probability of inclusion
  # according to the Carbonetto Stephens model
  idmod <- varbvs(X_mat, y, Z = Z[, 1], verbose = FALSE)
  sigmasq <- mean(idmod$sigma)
  pi_est <- mean(1 / (1 + exp(-idmod$logodds)))

  # values of sigmabeta_sq to optimize according to ELBO
  sigmavec <- c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10)

  # vector for storing the ELBO for each value of sigma in sigmavec
  elb1 <- matrix(0, length(sigmavec), 1)

  # loop to optimize sigma
  for (j in 1:length(sigmavec)) {
    res <- cov_vsvb(y, Z, X_mat, sigmasq, sigmavec[j], S_sq, pi_est, mu_mat,
                    alpha_mat)
    elb1[j] <- res$var.elbo
  }

  # Select the value of sigma_beta that maximizes the elbo
  sigmabeta_sq <- sigmavec[which.max(elb1)]

  # fit another model using this value of sigma_beta
  result <- cov_vsvb(y, Z, X_mat, sigmasq, sigmabeta_sq, S_sq, pi_est, mu_mat,
                     alpha_mat)

  # n by p matrix; the i,j-th entry is the probability of inclusion for the
  # i-th individual for the j-th variable according to the regression on y
  graph_list[[resp_index]] <- result$var.alpha
}


# check to see that this modified code produces the same results as the original code
if (discrete_data){
  load("original_discrete_alpha_matrices.Rdata")
}else{
  load("original_continuous_alpha_matrices.Rdata")
}

same <- T
for (j in 1:length(graph_list)) {
  if (all.equal(graph_list[[j]], mylist[[j]]) != T) {
    same <- F
    break
  }
}
same

# stop timer and see how much time has elapsed
end_time <- Sys.time()
end_time - start_time

## History (continuous data):
  # Original:
    # Time difference of 23.49625 secs
    # Time difference of 23.76863 secs
  # Modified the s_update:
    # Time difference of 16.98189 secs
    # Time difference of 15.27206 secs
  # Modified the mu update:
    # Time difference of 14.16393 secs
    # Time difference of 14.3205 secs
  # Modified the alpha_update:
    # Time difference of 11.38467 secs
    # Time difference of 11.68092 secs
  # Re-organization (removing unnecessary variables)
    # Time difference of 10.42546 secs
    # Time difference of 10.87398 secs
