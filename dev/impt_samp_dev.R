setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev")
rm(list = ls())
source("generate_data.R")

dat <- generate_continuous()
data_mat <- dat$data
Z <- dat$covts

source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/weights.R")
Rcpp::sourceCpp("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/src/covdepGE_c.cpp")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/covdepGE_R.R")

# get weights
D <- get_weights(Z, 2, T, NA)
bandwidths <- D$bandwidths
D <- D$D

y <- data_mat[ , 1]
X <- data_mat[ , -1]
n <- nrow(X)
p <- ncol(X)
mu <- matrix(0, n, p)
alpha <- matrix(0.2, n, p)
ssq <- 0.5
sbsq <- 0.01
pip <- 0.2
elbo_tol <- 0.01
alpha_tol <- 1e-5
max_iter <- 100

# generate the grid of hyperparameters
ssq_grid <- exp(seq(log(0.1), log(4), length.out = 5))
sbsq_grid <- exp(seq(log(1e-5), log(0.2), length.out = 5))
pi_grid <- exp(seq(log(1e-3), log(0.5), length.out = 5))
hp <- expand.grid(ssq = ssq_grid, sbsq = sbsq_grid, pip = pi_grid)

# find the prior density for hyperparameters
hp$ssq_p <- dgamma(hp$ssq, 1, 1)
hp$sbsq_p <- dgamma(hp$sbsq, 9e-1, 1)
hp$pi_p <- dbeta(hp$pip, .4, .9)
hp$prior <- hp$ssq_p * hp$sbsq_p * hp$pi_p

# find the importance density for hyperparameters
hp$ssq_q <- 1 / length(ssq_grid)
hp$sbsq_q <- 1 / length(sbsq_grid)
hp$pi_q <- 1 / length(pi_grid)
hp$impt <- hp$ssq_q * hp$sbsq_q * hp$pi_q

# find the ratio of the prior density to the importance density for hyperparameters
hp$p_q <- hp$prior / hp$impt

# create lists for storing the alpha matrices and elbos
elbo_theta <- alpha_theta <- vector("list", nrow(hp))

# iterate over each of the hyperparameters
for (j in 1:nrow(hp)){

  # fix hyperparameter setting
  hp_j <- hp[j , c("ssq", "sbsq", "pip")]

  # perform CAVI for the hyperparameter setting
  out <- cavi_c(y, D, X, mu, alpha, hp_j$ssq, hp_j$sbsq, hp_j$pip, elbo_tol, alpha_tol, max_iter, F)

  # save the alpha matrices
  alpha_theta[[j]] <- out$alpha

  # calculate the elbo for each individual under the current hyperparameter
  # setting and save
  elbo_l <- rep(NA, n)
  for (l in 1:n){
    elbo_l[l] <- ELBO_calculator_c(y, D[ , l], X, t(out$ssq_var[l, ]), t(out$mu[l, ]), t(out$alpha[l, ]), hp_j$ssq, hp_j$sbsq, hp_j$pip)
  }

  # save the ELBO
  elbo_theta[[j]] <- elbo_l
}

# calculate the weights for each individual
for (l in 1:n){

  # retrieve the elbo for the l-th individual for each hyperparameter setting
  elbo_l <- sapply(elbo_theta, `[`, l)

  # subtract the maximum value from each and exponentiate to get the log-likelihood
  lik_l <- exp(elbo_l - max(elbo_l))

  # multiply by the ratio of prior and importance density to find the
  # unnormalized weight for each hyperparameter setting
  weight <- lik_l * hp$impt

  # normalize the weights
  weight <- weight / sum(weight)

  # multiply the l-th row of each of the alpha matrices by the corresponding weight
  for (k in 1:length(alpha_theta)){
    alpha_theta[[k]][l, ] <- weight[k] * alpha_theta[[k]][l, ]
  }
}

# sum across the alpha matrices to obtain the final alpha matrix
alpha <- Reduce("+", alpha_theta)

