library(microbenchmark)

#-------------------------------------------------------------------------------
#--------------------------------DATA CREATION----------------------------------
#-------------------------------------------------------------------------------
rm(list = ls())

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

# Data generation
n <- 180
p <- 4

# function to create sigma matrix for a p-dimensional gaussian given the value
# of an extraneous covariate
Var_cont <- function(z) {
  STR <- 1
  pr <- matrix(0, p + 1, p + 1)
  diag(pr) <- 2
  pr[2, 3] <- STR
  pr[1, 2] <- STR * ((z > -1) && (z < -.33)) + (STR - STR * ((z + .23) / .56)) * ((z > -0.23) && (z < 0.33)) + (0) * ((z > 0.43) && (z < 1))
  pr[1, 3] <- 0 * ((z > -1) && (z < -.33)) + (STR * ((z + .23) / .56)) * ((z > -0.23) && (z < 0.33)) + (STR) * ((z > 0.43) && (z < 1))
  pr[2, 1] <- pr[1, 2]
  pr[3, 1] <- pr[1, 3]
  pr[3, 2] <- pr[2, 3]
  Var <- solve(pr)
  return(Var)
}

# creation of covariate
Z <- c(
  seq(-0.99, -0.331, (-.331 + .99) / 59),
  seq(-0.229, 0.329, (.329 + .229) / 59),
  seq(0.431, .99, (.99 - .431) / 59)
)
Z <- matrix(Z, n, 1)

# creation the data matrix; each individual is generated from a MVN with 0 mean
# and covariance matrix determined by their corresponding extraneous covariate
data_mat <- matrix(0, n, p + 1)
for (i in 1:n) {
  data_mat[i, ] <- mvrnorm(1, rep(0, p + 1), Var_cont(Z[i]))
}

# D is an n by n matrix of weights
D <- matrix(1, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    D[j, i] <- dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, 0.56)
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

# big ind matrix is a matrix of p stacked I_p identities
Big_ind <- matrix(rep(diag(p), p), n * p, p, T)

# main loop
resp_index <- 1
  
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
X <- matrix(rep(0, n^2 * p), nrow = n, ncol = n * p)

for (i in 1:n) {
  for (j in 1:p) {
    k <- p * (i - 1) + 1
    X[i, k + j - 1] <- X_mat[i, j]
    X_vec[k + j - 1] <- X[i, k + j - 1]
  }
}

Xty <- t(X) %*% y

# a block diagonal matrix; the j-th block is the transpose of the j-th row of
# X times the j-th row of X; it is an n*p by n*p matrix
XtX <- t(X) %*% X

DXtX <- diag(XtX)
DXtX_rep <- rep(DXtX, p)
DXtX_mat <- matrix(DXtX_rep, n * p, p, byrow = FALSE)
DXtX_Big_ind <- matrix(DXtX_rep, n * p, p, byrow = FALSE) * Big_ind

# XtX with its diagonal removed and replaced with 0
Diff_mat <- XtX - diag(DXtX)

alpha <- rep(0.2, n * p)
sigmabeta_sq <- 3
sigmasq <- 1
S_sq <- matrix(sigmasq * (DXtX + 1 / sigmabeta_sq)^(-1), n, p)
mu <- rep(0, p)
mu_mat <- matrix(0, n, p, byrow = TRUE)

# y_long_vec is each element of y repeated p times
y_long_vec <- rep(y, each = p)


#################### tuning hyperparameters ##################################

# Setting hyperparameter value as in Carbonetto Stephens model
idmod <- varbvs(X_mat, y, Z = Z[, 1], verbose = FALSE)
sigmasq <- mean(idmod$sigma)

pi_est <- mean(1 / (1 + exp(-idmod$logodds)))
sigmavec <- c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10)

# vector for storing the ELBO for each value of sigma in sigmavec
elb1 <- matrix(0, length(sigmavec), 1)
ELBO_LBit <- rep(0, 10000)


#-------------------------------------------------------------------------------
#--------------------------s_sq UPDATE OPTIMIZATION-----------------------------
#-------------------------------------------------------------------------------

## S_sq update - original
# sigmasq - scalar
# DXtX_Big_ind - n*p rows by p columns
S_sq_orig <- function(){
  for (i in 1:n) {
    S_sq[i, ] <- sigmasq * (t(DXtX_Big_ind) %*% D_long[, i] + 1 / sigmabeta_sq)^(-1)
    # variance parameter
  }
  return(S_sq)
}
S_sq <- S_sq_orig()
S_sq_vec <- matrix(t(S_sq), n * p, 1)

## S_sq update - modified 1 (DXtX_Big_ind) %*% D_long)
S_sq2_update <- function() {
  return(t(sigmasq * (t(DXtX_Big_ind) %*% D_long + 1 / sigmabeta_sq)^(-1)))
}
S_sq2 <- S_sq2_update()
all(S_sq2 == S_sq)

## S_sq update - modified 2 (t(X_mat^2) %*% D)
S_sq3_update <- function(){
  return(t(sigmasq * (t(X_mat^2) %*% D + 1 / sigmabeta_sq)^(-1)))
}

S_sq3 <- S_sq3_update()
all(S_sq3 == S_sq)

microbenchmark(S_sq_orig(),
               S_sq2_update(),
               S_sq3_update())

#-------------------------------------------------------------------------------
#--------------------------mu UPDATE OPTIMIZATION-----------------------------
#-------------------------------------------------------------------------------

## mu update - original 
for (i in 1:n) {
  y_XW <- y_long_vec * X_vec * D_long[, i]
  y_XW_mat <- matrix(y_XW, n, p, byrow = TRUE)
  
  X_mu_alpha <- X_vec * Mu_vec * alpha_vec
  xmualpha_mat <- t(matrix(X_mu_alpha, p, n)) %*% (matrix(1, p, p) - diag(rep(1, p)))
  XW_mat <- matrix(X_vec * D_long[, i], n, p, byrow = TRUE) * xmualpha_mat
  
  mu_mat[i, ] <- (t(y_XW_mat) %*% matrix(1, n, 1) - (t(XW_mat) %*% matrix(1, n, 1))) * (S_sq[i, ] / sigmasq) ### ### CAVI updation of mean variational parameter mu
}
Mu_vec <- matrix(t(mu_mat), n * p, 1)