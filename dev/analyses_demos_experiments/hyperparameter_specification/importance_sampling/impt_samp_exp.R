start <- Sys.time()
library(covdepGE)
library(ggplot2)
library(doParallel)
library(doRNG)

# function for generating the data and the covariates
generate_continuous <- function(n1 = 60, n2 = 60, n3 = 60, p = 4){

  # create covariate for individuals in each of the three intervals

  # define the dimensions of the data
  n <- sum(n1, n2, n3)

  # define the limits of the intervals
  limits1 <- c(-3, -1)
  limits2 <- c(-1, 1)
  limits3 <- c(1, 3)

  # define the covariate values within each interval
  z1 <- runif(n1, limits1[1], limits1[2])
  z2 <- runif(n2, limits2[1], limits2[2])
  z3 <- runif(n3, limits3[1], limits3[2])
  Z <- matrix(sort(c(z1, z2, z3)), n, 1)

  # create precision matrices

  # the shared part of the structure for all three intervals is a 2 on the
  # diagonal and a 1 in the (2, 3) position
  common_str <- diag(p + 1)
  common_str[2, 3] <- 1

  # define constants for the structure of interval 2
  beta1 <- diff(limits2)^-1
  beta0 <- -limits2[1] * beta1

  # interval 2 has two different linear functions of Z in the (1, 2) position
  # and (1, 3) positions; define structures for each of these components
  int2_str12 <- int2_str13 <- matrix(0, p + 1, p + 1)
  int2_str12[1, 2] <- int2_str13[1, 3] <- 1

  # define the precision matrices for each of the individuals in interval 2
  int2_prec <- lapply(z2, function(z) common_str +
                        ((1 - beta0 - beta1*z)*int2_str12) +
                        ((beta0 + beta1*z)*int2_str13))

  # interval 1 has a 1 in the (1, 2) and interval 3 has a 1 in the (1, 3) position;
  # define structures for each of these components
  int1_str12 <- int3_str13 <- matrix(0, p + 1, p + 1)
  int1_str12[1, 2] <- int3_str13[1, 3] <- 1

  # define the precision matrices for each of the individuals in interval 1 and interval 3
  int1_prec <- rep(list(common_str + int1_str12), n1)
  int3_prec <- rep(list(common_str + int3_str13), n3)

  # put all of the precision matrices into one list
  prec_mats <- c(int1_prec, int2_prec, int3_prec)

  # symmetrize the precision matrices
  prec_mats <- lapply(prec_mats, function(mat) t(mat) + mat)

  # invert the precision matrices to get the covariance matrices
  cov_mats <- lapply(prec_mats, solve)

  # generate the data using the covariance matrices
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))

  return(list(data = data_mat, covts = Z, true_precision = prec_mats))
}

# function for calculating the log-odds
log_odds <- function(p) log(p / (1 - p))

# function for calculating the prior density of sbsq given log-odds pi and the
# explanatory data
p_sbsq <- function(t, lo_pi, dat){

  # calculate the sum of the variances of the columns of dat
  var_sum <- sum(diag(var(dat)))

  # calculate the value of pi corresponding to the log-odds
  mu_pi <- exp(lo_pi) / (1 + exp(lo_pi))

  # calculate the density
  (pi * var_sum) / (1 + t * pi * var_sum)^2 * (t > 0)
}

# list for storing the results of each of the trials
n_trials <- 10
set.seed(1)

# spin up parallel backend
# find the available number of cores
tot_cores <- detectCores() - 10

# number of cores needed per trial
trial_cores <- 6

# how many parallel executions can be run in parallel?
(outer_cores <- floor(tot_cores / trial_cores))

registerDoParallel(outer_cores)

res <- foreach(trial_ind = 1:n_trials,
               .packages = c("covdepGE", "doParallel")) %dorng%{

                 # spin up the cores for the current trial
                 registerDoParallel(trial_cores)

                 # generate the data
                 dat <- generate_continuous()
                 X <- dat$data
                 Z <- dat$covts

                 # get the dimensions of the data
                 n <- nrow(X)
                 p <- ncol(X) - 1

                 # find the parameters of the prior on the log-odds:
                 mu_lo <- log_odds(1 / p)
                 sd_lo <- abs(log_odds((0.5 * p - 1) / p))

                 # # plot the prior on PIP
                 # ggplot() + geom_function(fun = function(x) dnorm(log_odds(x), mu_lo, sd_lo))
                 #
                 # # plot the prior on sbsq with pi : log_odds pi = mu_lo and response is variable 1
                 # ggplot() + geom_function(fun = p_sbsq, args = list(lo_pi = mu_lo, dat = X[ , -1])) + xlim(c(0, 10))

                 # generate the grid of hyperparameters for each of the variables fixed as the
                 # response
                 # store the grid for each in an m x (p + 1) matrix, where m is the number of
                 # hyperparameter settings
                 # also store the prior densities for each of the settings
                 ssq_upper_mult <- 10
                 nssq <- nsbsq <- npip <- 10
                 m <- prod(nssq, nsbsq, npip)
                 ssq_hp <- ssq_p_mat <- sbsq_hp <- sbsq_p_mat <- pip_hp <- pip_p_mat <- matrix(NA, m, p + 1)
                 for (resp_index in 1:(p + 1)){

                   # fix the resp_index column of the data as the response
                   y <- X[ , resp_index]
                   X_j <- X[ , -resp_index]

                   # find the upper bound for the grid of ssq
                   ssq_upper <- ssq_upper_mult * var(y)

                   # find the lasso estimate to the proportion of non-zero coefficients
                   lasso <- glmnet::cv.glmnet(X_j, y)
                   non0 <- sum(coef(lasso, s = "lambda.min")[-1] != 0)
                   non0 <- max(non0, 1)
                   pi_hat <- non0 / p

                   # find an upper bound for the grid of pi
                   pi_upper <- min(0.5, 2 * pi_hat)

                   # find the sum of the variances for each of the columns of X_mat
                   s2_sum <- sum(apply(X_j, 2, var))

                   # find the upper bound for the grid of sbsq
                   sbsq_upper <- 25 / (pi_hat * s2_sum)

                   # create the grid candidates for ssq and sbsq
                   var_lower <- 1e-5
                   ssq <- exp(seq(log(var_lower), log(ssq_upper), length.out = nssq))
                   sbsq <- exp(seq(log(var_lower), log(sbsq_upper), length.out = nsbsq))

                   # create posterior inclusion probability grid
                   pip <- exp(seq(log(0.01), log(pi_upper), length.out = npip))

                   # create the hyperparameter grid; find the prior density for each
                   # hyperparameter setting
                   hp <- expand.grid(ssq = ssq, sbsq = sbsq, pip = pip)
                   ssq_p <- 1 / length(ssq)
                   pip_p <- dnorm(log_odds(hp$pip), mu_lo, sd_lo)
                   sbsq_p <- sapply(1:nrow(hp), function(hp_ind) p_sbsq(hp[hp_ind, "sbsq"],
                                                                        hp[hp_ind, "pip"], X_j))
                   prior_p <- ssq_p * sbsq_p * pip_p

                   # store the hyperparameters and the prior densities in the resp_index column
                   # of the corresponding storage matrix
                   ssq_hp[ , resp_index] <- hp$ssq
                   sbsq_hp[ , resp_index] <- hp$sbsq
                   pip_hp[ , resp_index] <- hp$pip
                   ssq_p_mat[ , resp_index] <- ssq_p
                   sbsq_p_mat[ , resp_index] <- sbsq_p
                   pip_p_mat[ , resp_index] <- pip_p
                 }

                 # perform inference with the priors specified above
                 out_prior <- covdepGE(X, Z,
                                       ssq = ssq_hp, ssq_p = ssq_p_mat,
                                       sbsq = sbsq_hp, sbsq_p = sbsq_p_mat,
                                       pip = pip_hp, pip_p = pip_p_mat,
                                       grid_search = F, parallel = T, stop_cluster = F, warnings = F)

                 # perform inference with flat priors
                 out_unif_prior <- covdepGE(X, Z,
                                            ssq = ssq_hp,
                                            sbsq = sbsq_hp,
                                            pip = pip_hp,
                                            grid_search = F, parallel = T,
                                            stop_cluster = F, warnings = F)

                 # perform inference with grid search
                 out_grid <- covdepGE(X, Z,
                                      ssq = ssq_hp,
                                      sbsq = sbsq_hp,
                                      pip = pip_hp,
                                      grid_search = T, parallel = T,
                                      stop_cluster = F, warnings = F)

                 # return each of the resulting models and the data
                 list(data = dat, models = list(
                   prior = out_prior, unif_prior = out_unif_prior, grid = out_grid
                 ))
               }

# save res
save(res, file = "impt_samp_models_elbo_sum.Rda")
Sys.time() - start
