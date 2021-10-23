setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/Development")

#-------------------------------------------------------------------------------
#--------------------------------DATA CREATION----------------------------------
#-------------------------------------------------------------------------------
rm(list = ls())

source("generate_data.R")
library(microbenchmark)
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

# if a toy example is desired, data_mat and Z are replaced here
toy <- F
if (toy) {
  set.seed(3)
  n <- 5
  p <- 3
  data_mat <- matrix(sample(-3:3, n * (p + 1), T), n, p + 1)
  Z <- matrix(sample(1:n, n), n, 1)
}

# D is an n by n matrix of weights
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

# if a toy example is desired, make D nicer numbers
if (toy){
  set.seed(1)
  D <- matrix(sample(1:5, n * n, T), n, n)
  D <- round(n * D / rowSums(D), 1)
}

# the i-th column of D_long is the i-th column of D with the elements repeated
# p times
D_long <- matrix(rep(D, each = p), n * p)

# The variable specific inclusion probability matrix:
# i-th row corresponds to the dependence structure for the i-th subject,
# j-th matrix corresponds to the j th variable as response and the remaining as
# predictors.
mylist <- vector("list", p + 1)

# big ind matrix is a matrix of n stacked I_p identities
Big_ind <- matrix(rep(diag(p), n), n * p, p, T)

# main loop
resp_index <- 2

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
#--------------------------mu UPDATE OPTIMIZATION-------------------------------
#-------------------------------------------------------------------------------

## mu update - original
original_mupdate <- function(){
  Mu_vec <- matrix(alpha, n * p, 1)
  alpha_vec <- matrix(alpha, n * p, 1, byrow = TRUE)

  for (i in 1:n) {
    y_XW <- y_long_vec * X_vec * D_long[, i]
    y_XW_mat <- matrix(y_XW, n, p, byrow = TRUE)

    X_mu_alpha <- X_vec * Mu_vec * alpha_vec
    xmualpha_mat <- t(matrix(X_mu_alpha, p, n)) %*% (matrix(1, p, p) - diag(rep(1, p)))
    XW_mat <- matrix(X_vec * D_long[, i], n, p, byrow = TRUE) * xmualpha_mat

    mu_mat[i, ] <- (t(y_XW_mat) %*% matrix(1, n, 1) - (t(XW_mat) %*% matrix(1, n, 1))) * (S_sq[i, ] / sigmasq) ### ### CAVI updation of mean variational parameter mu
  }
  Mu_vec <- matrix(t(mu_mat), n * p, 1)
  mu_mat
}

## Mu update - modified 1

modified_mupdate <- function(){
  alpha_mat <- matrix(alpha, n, p, byrow = TRUE)
  mu_mat2 <- matrix(alpha, n, p)

  mu_mat2.copy <- mu_mat2

  for (l in 1:n){

    # vector to ensure that the k-th entry of the l-th row of mu_mat and alpha
    # mat are not summed
    ind0 <- rep(1, p)

    # update the l, k entry of mu
    for (k in 1:p){

      # put a 0 in the k-th position of ind0 to 0 out the k-th entry of mu_stack
      # and alpha_stack
      ind0.k <- ind0; ind0.k[k] <- 0

      # the l-th row of mu_mat, alpha_mat with the k-th entry 0'd, stacked n times
      mu_stack <- matrix(mu_mat2.copy[l, ] * ind0.k, n, p, T)
      alpha_stack <- matrix(alpha_mat[l, ] * ind0.k, n, p, T)

      # update mu
      mu_mat2[l, k] <- ((S_sq[l, k] / sigmasq) *
                          sum((D[ , l] * X_mat[ , k]) *
                                (y - rowSums(X_mat * mu_stack * alpha_stack))))
    }
  }
  mu_mat2
}

## Mu update - modified 2

#modified_mupdate2 <- function(){
alpha_mat <- matrix(1:(n*p), n, p, byrow = TRUE)
mu_mat3 <- matrix(1:(n*p), n, p)

mu_mat3.copy <- mu_mat3

for (l in 1:n){

  # the l-th row of mu_mat, alpha_mat stacked n times
  mu_stack <- matrix(mu_mat3.copy[l, ], n, p, T)
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
  mu_mat3[l, ] <- S_sq[l, ] / sigmasq * colSums(d_x_y)

}
mu_mat3
#}

all.equal(original_mupdate(), modified_mupdate2())
all.equal(modified_mupdate2(), modified_mupdate())


microbenchmark::microbenchmark(original_mupdate(),
                               modified_mupdate(),
                               modified_mupdate2(), times = 500)

#-------------------------------------------------------------------------------
#--------------------------alpha UPDATE OPTIMIZATION----------------------------
#-------------------------------------------------------------------------------


# alpha update - original

# create necessary variables
mu_mat <- matrix(1:(n * p), n, p)
Mu_vec <- matrix(t(mu_mat), n * p, 1)
S_sq_vec <- matrix(t(S_sq), n * p, 1)
thres <- 1e-7
lthres <- logit(thres)
uthres <- logit(1 - thres)

vec_1 <- log(pi_est / (1 - pi_est)) # first term of alpha update
vec_2 <- as.matrix(0.5 * log(S_sq_vec / (sigmasq * sigmabeta_sq))) # second term of alpha update
vec_3 <- as.matrix(Mu_vec^2 / (2 * S_sq_vec)) # third term of alpha update

# logit of alpha is sum of 3 terms for alpha update
unlogitalpha <- vec_1 + vec_2 + vec_3

# find values of alpha that are too large/ small and will pose an issue for
# the reverse logit formula; set them to the upperthreshold/ lower
# threshold values
indlarge <- which(unlogitalpha > uthres)
indsmall <- which(unlogitalpha < lthres)
unlogitalpha[indlarge] <- uthres
unlogitalpha[indsmall] <- lthres

alpha[which(unlogitalpha > 9)] <- 1 # thresholding very large values to 1 for computational stability
alpha[which(unlogitalpha <= 9)] <- 1 / (1 + exp(-unlogitalpha[which(unlogitalpha <= 9)])) ### ### CAVI updation of variational parameter alpha


# modified alpha update
alpha_mt <- matrix(alpha, n, p, T)


logit_ <- logit(pi_est) + (mu_mat^2 / (2 * S_sq)) + 0.5 * log(S_sq / (sigmasq * sigmabeta_sq))
logit2 <- logit(pi_est) + (mu_mat^2 / (2 * S_sq)) + log(sqrt(S_sq / (sigmasq * sigmabeta_sq)))
all(logit_ == logit2)
exp_logit <- exp(logit)
alpha_mat2 <- exp_logit / (1 + exp_logit)
alpha_mat2[is.infinite(exp_logit)] <- 1
all.equal(alpha_mat2, alpha_mt)


#-------------------------------------------------------------------------------
#--------------------------------ELBO CALCULATION-------------------------------
#-------------------------------------------------------------------------------

ELBO_calculator(y, X_mat, S_sq[i, ], mu_mat[i, ], alpha_mat[i, ], sigmasq, sigmabeta_sq, true_pi, D[, i], n, p)
