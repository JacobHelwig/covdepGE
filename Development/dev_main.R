set.seed(1)
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/Development")
rm(list = ls())
library(Rcpp)
source("generate_data.R")
sourceCpp("c_dev.cpp")
start_time <- Sys.time()

source("generate_data.R")
library(MASS)
library(varbvs)

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

  # instantiate initial values of variational parameters and hyperparmeters
  alpha_mat <- matrix(0.2, n, p)
  sigmabeta_sq <- 3
  sigmasq <- 1
  E <- rnorm(n, 0, sigmasq) # removing this causes discrepency in discrete case
  # dont delete S_sq <- matrix(sigmasq * (DXtX + 1 / sigmabeta_sq)^(-1), n, p) # should be byrow = T?
  S_sq <- t(sigmasq * (t(X_mat^2) + 1 / sigmabeta_sq)^(-1))
  mu_mat <- matrix(0, n, p, byrow = TRUE)

  # Setting hyperparameter values for sigmasq and the probability of inclusion
  # according to the Carbonetto Stephens model
  idmod <- varbvs(X_mat, y, Z = Z[, 1], verbose = FALSE)
  sigmasq <- mean(idmod$sigma)
  pi_est <- mean(1 / (1 + exp(-idmod$logodds)))

  # values of sigmabeta_sq to optimize according to ELBO
  sigmavec <- c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10)

  # loop to optimize sigma; for each value of sigma in sigmavec, store the
  # resulting ELBO from fitting n graphs using that value as sigmabeta_sq
  elbo_sigma <- sigma_loop_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmavec, pi_est)

  # Select the value of sigma_beta that maximizes the elbo
  sigmabeta_sq <- sigmavec[which.max(elbo_sigma)]

  # fit another model using this value of sigma_beta
  result <- cov_vsvb_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq, pi_est)

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
  # ELBO calculation in C++
    # Time difference of 9.042475 secs
    # Time difference of 9.159554 secs
  # Mu update in C++
    # Time difference of 4.533927 secs
    # Time difference of 4.517517 secs
  # Alpha update in C++
    # Time difference of 4.118175 secs
    # Time difference of 4.188726 secs
  # Move ELBO calculation outside the variational update loop
    # Time difference of 2.080118 secs
    # Time difference of 1.93863 secs
  # Variational update loop (cov_vsvb function) to C++
    # Time difference of 1.907178 secs
    # Time difference of 1.922939 secs
  # Sigma loop to C++
    # Time difference of 1.895749 secs
    # Time difference of 1.865739 secs
