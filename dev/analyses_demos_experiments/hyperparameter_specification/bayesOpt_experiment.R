setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev")
rm(list = ls())
source("generate_data.R")
Rcpp::sourceCpp("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/src/covdepGE_c.cpp")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/weights.R")
library(GPfit)

dat <- generate_continuous()

data <- dat$data
y <- data[ , 1]
X <- data[ , -1]
n <- nrow(X)
p <- ncol(X)
mu <- matrix(0, n, p)
alpha <- matrix(0.2, n, p)
Z <- dat$covts
D <- get_weights(Z, 2, T)$D

# get an estimate to the upper bound for variance hyperparameters
lasso <- glmnet::cv.glmnet(X, y)
non0 <- sum(coef(lasso, s = "lambda.min")[-1] != 0)
pi_hat <- non0 / p
var_sum <- sum(apply(X, 2, var))
sbsq_upper <- 25 / (pi_hat * var_sum)
ssq_upper <- 10 * var(y)

# create an initial grid to get ELBO on
var_lower <- 1e-3
sbsq0 <- seq(var_lower, sbsq_upper / 10, length.out = 3)
ssq0 <- seq(var_lower, ssq_upper / 10, length.out = 3)
pip_lower <- 0.01
pip_upper <- 0.99
pip0 <- seq(pip_lower, pip_upper / 2, length.out = 3)
grid0 <- expand.grid(ssq = ssq0, sbsq = sbsq0, pip = pip0)

# find ELBO at each point on the initial grid
elbos <- rep(NA, nrow(grid0))
for (j in 1:nrow(grid0)){
  print(grid0[j, ])
  ssq_j <- grid0$ssq[j]
  sbsq_j <- grid0$sbsq[j]
  pip_j <- grid0$sbsq[j]
  out_j <- cavi_c(y, D, X, mu, alpha, ssq_j, sbsq_j, pip_j, 1, 1e-12, 1e3, T)
  elbos[j] <- out_j$elbo
}

grid0 <- cbind.data.frame(grid0, elbos)

# filter out the elbo that are NaN
grid0 <- grid0[!is.na(grid0$elbos), ]

# create a grid of candidate points
d <- 10
ssq_cand <- exp(seq(log(var_lower), log(ssq_upper), length.out = d))
sbsq_cand <- exp(seq(log(var_lower), log(sbsq_upper), length.out = d))
pip_cand <- exp(seq(log(pip_lower), log(pip_upper), length.out = d))
grid_cand <- expand.grid(ssq = ssq_cand, sbsq = sbsq_cand, pip = pip_cand)

# normalize the candidates to [0, 1]
ssq_cand01 <- (grid_cand$ssq - var_lower) / (ssq_upper - var_lower)
sbsq_cand01 <- (grid_cand$sbsq - var_lower) / (sbsq_upper - var_lower)
pip_cand01 <- (grid_cand$pip - pip_lower) / (pip_upper - pip_lower)
grid_cand01 <- cbind(ssq_cand01, sbsq_cand01, pip_cand01)

# normalize the initial grid values to 01
ssq01 <- (grid0$ssq - var_lower) / (ssq_upper - var_lower)
sbsq01 <- (grid0$sbsq - var_lower) / (sbsq_upper - var_lower)
pip01 <- (grid0$pip - pip_lower) / (pip_upper - pip_lower)
grid0 <- cbind(grid0, ssq01, sbsq01, pip01)

# fit a GP regression to obtain a surrogate to the ELBO as a function of the
# hyperparameters
gp_fit <- GP_fit(grid0[ , c("ssq01", "sbsq01", "pip01")], grid0$elbos)

# get expected value and sd at each of the candidate grid points
gp_pred <- predict(gp_fit, grid_cand01)
mu_gp <- gp_pred$Y_hat
sigma_gp <- sqrt(gp_pred$MSE)

# get the expected improvement at each of the grid points
y_best <- max(grid0$elbos)
gamma <- ifelse(sigma_gp == 0, 0, (mu_gp - y_best) / sigma_gp)
exp_imp <- sigma_gp * (gamma * pnorm(gamma) + dnorm(gamma))

# find the next value of the hyperparamters
next_hp01 <- as.numeric(grid_cand01[which.max(exp_imp), ])
next_hp <- as.numeric(grid_cand[which.max(exp_imp), ])
next_elbo <- cavi_c(y, D, X, mu, alpha, next_hp01[1], next_hp01[2], next_hp01[3], 1, 1e-12, 1e3, T)
next_elbo <- next_elbo$elbo
grid0 <- rbind(grid0, c(next_hp, next_elbo, next_hp01))

bayesOpt2D <- function(max_iter = 10){

  # create a grid of values to evaluate the surrogate to the loss
  surr_grid <- seq(0, 1, length.out = 100)
  surr_grid <- expand.grid(surr_grid, surr_grid)

  # create the initial grid and normalize the X values to [0, 1]
  eval <- data.frame(beta1 = -3:3, beta2 = -3:3)
  eval <- rbind.data.frame(eval, MLE)
  beta_lim <- c(-3, 3)
  eval$beta1_norm <- with(eval, (beta1 - beta_lim[1]) / (beta_lim[2] - beta_lim[1]))
  eval$beta2_norm <- with(eval, (beta2 - beta_lim[1]) / (beta_lim[2] - beta_lim[1]))
  eval$y <- apply(eval, 1, function(par) log_lik(c(par[1], par[2]), mod_df))

  for (j in 1:max_iter){

    # fit a GP to obtain a surrogate to the loss and plot
    gp_fit <- GP_fit(X = eval[ , c("beta1_norm", "beta2_norm")], Y = eval$y,
                     corr = list(type = "exponential", power = 1.95))
    plot(gp_fit)

    # find the expected value of the surrogate at each of the grid points
    gp_preds <- predict(gp_fit, surr_grid)
    mu_gp <- gp_preds$Y_hat
    sig_gp <- sqrt(gp_preds$MSE)

    # find the expected improvement for each of the grid points
    y_best <- max(eval$y)
    gamma <- ifelse(sig_gp == 0, 0, (mu_gp - y_best) / sig_gp)
    expc_imp <- sig_gp * (gamma * pnorm(gamma) + dnorm(gamma))

    # choose the next point to evaluate as the point corresponding to the greatest
    # expected improvement; unnormalize from [0, 1] and find the loss at this point
    beta_norm_next <- surr_grid[which.max(expc_imp), ]
    while(any(apply(eval, 1, function(beta) all(c(beta[1], beta[2]) == beta_norm_next)))){

      # if beta_norm_next has already been evaluated, choose the point with the next
      # greatest value
      k <- 2
      beta_norm_next <- surr_grid[order(expc_imp, decreasing = T)[k], ]
      k <- k + 1
    }
    beta_norm_next <- as.numeric(beta_norm_next)
    beta_next <- beta_norm_next * (beta_lim[2] - beta_lim[1]) + beta_lim[1]
    y_next <- log_lik(beta_next, mod_df)
    eval <- rbind.data.frame(eval, c(beta_next, beta_norm_next, y_next))
  }
  return(eval)
}

bayesOpt2D(10)
