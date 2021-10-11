rm(list = ls())
source("cov_vsvb.R")
source("ELBO_calculator.R")
source("generate_data.R")
library(reshape2)
library(MASS)
library(varbvs)
library(ks)

logit <- function(x) {
  if ((x == 0) | (x == 1)) {
    return(0)
  } else {
    return(log(x / (1 - x)))
  }
}

# generate data and covariates
dat <- generate_discrete()
data_mat <- dat$data
Z <- dat$covts

# Data generation
n <- 100
p <- 10

# D is an n by n matrix of weights
D <- matrix(1, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    D[i, j] <- dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, .1)
  }
}

# Scale weights to sum to n
for (i in 1:n) {
  D[, i] <- n * (D[, i] / sum(D[, i]))
}

# the i-th column of D_long is the i-th column of D with the elements repeated
# p times
D_long <- matrix(rep(D, each = p), n * p)

# The variable specific inclusion probability matrix:
# i-th row corresponds to the dependence structure for the i-th subject,
# j-th matrix corresponds to the j th variable as response and the remaining as
# predictors.
mylist <- vector("list", p + 1)

Adj_Mat_vb <- array(0, dim = c(p + 1, p + 1))

# big ind matrix is a matrix of n stacked I_p identities
Big_ind <- matrix(rep(diag(p), n), n * p, p, T)

# main loop
for (resp_index in 1:(p + 1)) {

  # Set variable number `resp_index` as the response
  y <- data_mat[, resp_index]

  # Set the remaining p variables as predictor
  X_mat <- data_mat[, -resp_index]

  # X_vec is a vector of length n*p that is the rows of X_mat "unravelled" by row;
  # that is, the first element of X_vec is the 1,1 of X_mat; the second element
  # is the 1,2; third is 1,3, ect.
  X_vec <- matrix(0, n * p, 1)

  # X is a n by n*p matrix; it consists of rbinding n n by p matrices together
  # the j-th matrix is the j-th row of X_mat in the j-th row, and 0's o.w.
  X <- matrix(0, nrow = n, ncol = n * p)

  for (i in 1:n) {
    for (j in 1:p) {
      k <- p * (i - 1) + 1
      X[i, k + j - 1] <- X_mat[i, j]
      X_vec[k + j - 1] <- X[i, k + j - 1]
    }
  }

  ELBO_LBit <- rep(0, 10000)

  # big diag mat looks exactly the same as X, except it replaces all non-zero
  # entries with 1
  Big_diag_mat <- (X != 0) * 1

  q <- matrix(2, n, 1)

  sigmasq <- 1
  E <- rnorm(n, 0, sigmasq) # removing this causes discrepency

  # a block diagonal matrix; the j-th block is the transpose of the j-th row of
  # X times the j-th row of X; it is an n*p by n*p matrix
  XtX <- t(X) %*% X

  DXtX <- diag(XtX)
  DXtX_rep <- rep(DXtX, p)
  DXtX_mat <- matrix(DXtX_rep, n * p, p)

  # XtX with its diagonal removed and replaced with 0
  Diff_mat <- XtX - diag(DXtX)

  # Initialization of the inclusion probability matrix for a fixed variable
  # with i-th row corresponding to i th subject.
  alpha <- rep(0.2, n * p)

  sigmabeta_sq <- 3
  mu <- rep(0, p)
  true_pi <- 0.5

  # y_long_vec is each element of y repeated p times
  y_long_vec <- rep(y, each = p)

  Xty <- t(X) %*% y
  beta_mat <- matrix(0, n, p, byrow = TRUE)
  mu_mat <- beta_mat

  S_sq <- matrix(sigmasq * (DXtX + 1 / sigmabeta_sq)^(-1), n, p)

  iter <- 1

  DXtX_Big_ind <- DXtX_mat * Big_ind

  candL <- seq(0.1, 0.9, .2) # Different values of hyperparameter true_pi
  elb <- rep(0, length(candL))

  est_q <- rep(0, n)
  beta_matr <- matrix(0, n, p)

  #################### tuning hyperparameters ##################################

  # Setting hyperparameter value as in Carbonetto Stephens model
  idmod <- varbvs(X_mat, y, Z = Z[, 1], verbose = FALSE)
  #rest_index_set <- setdiff(c(1:(p + 1)), resp_index)
  sigmasq <- mean(idmod$sigma)
  pi_est <- mean(1 / (1 + exp(-idmod$logodds)))
  sigmavec <- c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10)

  # vector for storing the ELBO for each value of sigma in sigmavec
  elb1 <- matrix(0, length(sigmavec), 1)

  # loop to optimize sigma
  for (j in 1:length(sigmavec)) {
    res <- cov_vsvb(y, X, Z, XtX, DXtX, Diff_mat, Xty, sigmasq, sigmavec[j], pi_est)
    elb1[j] <- res$var.elbo
  }

  # Select the value of sigma_beta that maximizes the elbo
  sigmabeta_sq <- sigmavec[which.max(elb1)]

  # fit another model using this value of sigma_beta
  result <- cov_vsvb(y, X, Z, XtX, DXtX, Diff_mat, Xty, sigmasq, sigmabeta_sq, pi_est)

  # vector of length n * p of inclusion probabilities
  incl_prob <- result$var.alpha

  # n by p matrix; the i,j-th entry is the probability of inclusion for the
  # i-th individual for the j-th variable according to the regression on y
  heat_alpha <- matrix(incl_prob, n, p, byrow = TRUE)
  mylist[[resp_index]] <- heat_alpha
}

# check to see that this modified code produces the same results as the original code
mylist2 <- mylist
load("original_discrete_alpha_matrices.Rdata")
same <- T
for (j in 1:length(mylist)){
  if (all.equal(mylist[[j]], mylist2[[j]]) != T){
    same <- F
    break
  }
}
same
