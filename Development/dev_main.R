set.seed(1)
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/Development")
rm(list = ls())
library(Rcpp)
library(MASS)
library(varbvs)
source("generate_data.R")
sourceCpp("c_dev.cpp")

## _____________________________________________________________________________
## _____________________________covdepGE________________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## Function to estimate the covariance structure for n individuals using
## variational Bayes
## Optimizes the spike and slab paramter sigma_beta_sq by choosing the value
## for each predictor that maximizes the ELBO
## -----------------------------ARGUMENTS---------------------------------------
## data_mat: n by p data matrix
## Z: n by p' matrix of extraneous covariates
## tau: bandwidth parameter; greater values allow for more information to be
## shared from other individuals when estimating the graph for a fixed
## individual. 0.1 by default.
## default
## sigmavec: candidate values of sigmabeta_sq
## -----------------------------RETURNS-----------------------------------------
## TBD
## _____________________________________________________________________________
covdepGE <- function(data_mat, Z, tau = 0.1,
                     sigmavec = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10)){

  start_time <- Sys.time()

  # get sample size and number of parameters
  n <- nrow(data_mat); p <- ncol(data_mat) - 1

  # D is a symmetric n by n matrix of weights; the i, j entry is the similarity
  # between individuals i and j
  D <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      D[j, i] <- dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, tau)
      D[i, j] <- D[j, i]
    }
  }

  # Scale weights to sum to n
  D <- n * (D) * matrix(1 / colSums(D), n, n, T)

  # List for the variable-specific inclusion probability matrix; the j-th element
  # in the list is a n by p matrix corresponding to the k-th predictor;
  # in this matrix, the l-th row corresponds to the dependence structure for the
  # l-th individual, with the j-th predictor fixed as the response
  graph_list <- vector("list", p + 1)

  # main loop over the predictors
  for (resp_index in 1:(p + 1)) {

    # Set variable number `resp_index` as the response
    y <- data_mat[, resp_index]

    # Set the remaining p variables as predictors
    X_mat <- data_mat[, -resp_index]

    # instantiate initial values of variational parameters and hyperparmeters
    alpha_mat <- matrix(0.2, n, p) # argument?
    sigmabeta_sq <- 3 # argument?
    sigmasq <- 1 # argument?
    E <- rnorm(n, 0, sigmasq) # removing this causes discrepency in discrete case
    # dont delete S_sq <- matrix(sigmasq * (DXtX + 1 / sigmabeta_sq)^(-1), n, p) # should be byrow = T?
    S_sq <- t(sigmasq * (t(X_mat^2) + 1 / sigmabeta_sq)^(-1))
    mu_mat <- matrix(0, n, p, byrow = TRUE) # argument?

    # Setting hyperparameter values for sigmasq and the probability of inclusion
    # according to the Carbonetto Stephens model
    idmod <- varbvs(X_mat, y, Z = Z[, 1], verbose = FALSE)
    sigmasq <- mean(idmod$sigma)
    pi_est <- mean(1 / (1 + exp(-idmod$logodds)))

    # loop to optimize sigma; for each value of sigma in sigmavec, store the
    # resulting ELBO
    elbo_sigma <- sigma_loop_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmavec, pi_est)

    # Select the value of sigma_beta that maximizes the ELBO
    sigmabeta_sq <- sigmavec[which.max(elbo_sigma)]

    # fit another model using this value of sigma_beta
    result <- cov_vsvb_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq, pi_est)

    # n by p matrix; the i,j-th entry is the probability of inclusion for the
    # i-th individual for the j-th variable according to the regression on y
    graph_list[[resp_index]] <- result$var.alpha
  }

  # stop timer and see how much time has elapsed
  end_time <- Sys.time()
  print(end_time - start_time)

  return(graph_list)
}

# generate data and covariates
discrete_data <- F # true if discrete example is desired
if (discrete_data) {
  dat <- generate_discrete()
  tau_ <- 0.1 # the bandwidth parameter
}else{
  dat <- generate_continuous()
  tau_ <- 0.56
}

data_mat <- dat$data
Z.cov <- dat$covts

graphs <- covdepGE(data_mat, Z.cov, tau_)

# check to see that this modified code produces the same results as the original code
if (discrete_data){
  load("original_discrete_alpha_matrices.Rdata")
}else{
  load("original_continuous_alpha_matrices.Rdata")
}

same <- T
for (j in 1:length(graphs)) {
  if (all.equal(graphs[[j]], mylist[[j]]) != T) {
    same <- F
    break
  }
}
same

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
  # Modified calculation of weights
    # Time difference of 1.521307 secs
    # Time difference of 1.511709 secs
