set.seed(1)

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
## returns the final variational estimates
## for a single regression. So in the arguments, y plays the role of the
## response, i.e the j th variable whose
## conditional distribution given the remaining variables is being calculated.

## This calls the ELBO_calculator function

## X: the data matrix except the j th variable
## XtX is x transpose times x.
## DXtX:diagonal elements of XtX
## Diff_mat: XtX-diag(DXtX)
## Xty: x transpose times y
## sigmasq: variance of response given the parameters (homoscedastic part,
## actual variance sigma_sq/w_i)
## sigmabeta_sq: prior variance of coefficient parameter
## pi_est: estimate of spike and slab mixture proportion.

# fixed x_j; estimate variational parameters for each individual and
# calculate ELBO for sigma_beta optimization
cov_vsvb <- function(y, y_long_vec, Z, X, X_mat, XtX, X_vec, Xty, DXtX, 
                     DXtX_Big_ind, D_long, Diff_mat, sigmasq, sigmabeta_sq, 
                     S_sq, pi_est, mu, mu_mat, alpha, Big_ind, ELBO_LBit) {
  
  n <- nrow(X_mat); p <- ncol(X_mat)
  
  # threshold for calculating the reverse logit of alpha
  thres <- 1e-7
  lthres <- logit(thres)
  uthres <- logit(1 - thres)
  
  # exit condition tolerance
  tol <- 1e-9
  
  
  change_alpha <- rep(0.001, n * p) # alpha_new - alpha_int
  
  max_iter <- 100
  iter <- 1
  Mu_vec <- matrix(rep(mu, n), n * p, 1)
  
  iter <- 1
  
  # S_sq update (WHY IS THERE A NEED FOR ITERATION FOR THIS PARAMETER?)
  S_sq <- t(sigmasq * (t(X_mat^2) %*% D + 1 / sigmabeta_sq)^(-1))
  S_sq_vec <- matrix(t(S_sq), n * p, 1)
  
  # loop to optimize variational parameters
  while (sqrt(sum(change_alpha^2)) > tol & iter < max_iter) {
    alpha_int <- alpha
    
    alpha_mat <- matrix(alpha, n, p, byrow = TRUE)
    
    alpha_vec <- matrix(alpha, n * p, 1, byrow = TRUE)
    
    # S_sq update
    # for (i in 1:n) {
    #   S_sq[i, ] <- sigmasq * (t(DXtX_Big_ind) %*% D_long[, i] + 1 / sigmabeta_sq)^(-1) ## variance parameter
    # }
    # S_sq <- t(sigmasq * (t(X_mat^2) %*% D + 1 / sigmabeta_sq)^(-1))
    # S_sq_vec <- matrix(t(S_sq), n * p, 1)
    
    # mu update
    for (i in 1:n) {
      y_XW <- y_long_vec * X_vec * D_long[, i]
      y_XW_mat <- matrix(y_XW, n, p, byrow = TRUE)
      
      X_mu_alpha <- X_vec * Mu_vec * alpha_vec
      xmualpha_mat <- t(matrix(X_mu_alpha, p, n)) %*% (matrix(1, p, p) - diag(rep(1, p)))
      XW_mat <- matrix(X_vec * D_long[, i], n, p, byrow = TRUE) * xmualpha_mat
      
      mu_mat[i, ] <- (t(y_XW_mat) %*% matrix(1, n, 1) - (t(XW_mat) %*% matrix(1, n, 1))) * (S_sq[i, ] / sigmasq) ### ### CAVI updation of mean variational parameter mu
    }
    Mu_vec <- matrix(t(mu_mat), n * p, 1)
    
    # alpha update
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
    
    #
    alpha[which(unlogitalpha > 9)] <- 1 # thresholding very large values to 1 for computational stability
    alpha[which(unlogitalpha <= 9)] <- 1 / (1 + exp(-unlogitalpha[which(unlogitalpha <= 9)])) ### ### CAVI updation of variational parameter alpha
    
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
    
    alpha_new <- alpha
    change_alpha <- alpha_new - alpha_int
    
    ELBO_LBit[iter] <- ELBO_LB
    iter <- iter + 1
  }
  
  # how has ELBO evolved?
  ELBO_LBit <- ELBO_LBit[1:(iter - 1)]
  
  # return n times p - 1 of each of the variational parameters
  list(var.alpha = alpha, var.mu = mu_mat, var.S_sq = S_sq, var.elbo = ELBO_LB, var.elboit = ELBO_LBit)
}

# source("cov_vsvb.R")
# source("ELBO_calculator.R")
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
  #pi_est <- 0.5
  
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
  
  # loop to optimize sigma
  for (j in 1:length(sigmavec)) {
    res <- cov_vsvb(y, y_long_vec, Z, X, X_mat, XtX, X_vec, Xty, DXtX, 
                    DXtX_Big_ind, D_long, Diff_mat, sigmasq, sigmavec[j], 
                    S_sq, pi_est, mu, mu_mat, alpha, Big_ind, ELBO_LBit)
    elb1[j] <- res$var.elbo
  }
  
  # Select the value of sigma_beta that maximizes the elbo
  sigmabeta_sq <- sigmavec[which.max(elb1)]
  
  # fit another model using this value of sigma_beta
  result <- cov_vsvb(y, y_long_vec, Z, X, X_mat, XtX, X_vec, Xty, DXtX, 
                     DXtX_Big_ind, D_long, Diff_mat, sigmasq, sigmabeta_sq, 
                     S_sq, pi_est, mu, mu_mat, alpha, Big_ind, ELBO_LBit)

  # vector of length n * p of inclusion probabilities
  incl_prob <- result$var.alpha
  
  # n by p matrix; the i,j-th entry is the probability of inclusion for the
  # i-th individual for the j-th variable according to the regression on y
  heat_alpha <- matrix(incl_prob, n, p, byrow = TRUE)
  mylist[[resp_index]] <- heat_alpha
}

# check to see that this modified code produces the same results as the original code
mylist2 <- mylist
load("original_continuous_alpha_matrices.Rdata")
same <- T
for (j in 1:length(mylist)) {
  if (all.equal(mylist[[j]], mylist2[[j]]) != T) {
    same <- F
    break
  }
}
same

# stop timer and see how much time has elapsed
end_time <- Sys.time()
end_time - start_time

## History:
  # Original:
    # Time difference of 23.49625 secs
    # Time difference of 23.76863 secs
  # Modified the s_update:
    # Time difference of 16.98189 secs
    # Time difference of 15.27206 secs