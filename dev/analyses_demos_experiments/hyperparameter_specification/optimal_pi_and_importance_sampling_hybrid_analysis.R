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

# how many trials?
n_trials <- 100
set.seed(1)

# create the hyperparameter grid
# pip <- c(1e-5, 1e-4, 1e-3, 0.01, 0.025, seq(0.05, 0.5, 0.05))
pip <- seq(0.1, 0.5, 0.1)
sbsq <- c(1e-5, 1e-3, 1e-2, 0.1, 0.25, 0.5, 1, 3, 5, 10)
ssq <- c(1e-3, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 2, 3)
hp_grid <- expand.grid(ssq = ssq, sbsq = sbsq, pip = pip)

# spin up parallel backend
# find the available number of cores
cores <- detectCores() - 5

registerDoParallel(cores)

res <- foreach(trial_ind = 1:n_trials, .packages = "covdepGE") %dorng%{

                 # generate the data
                 dat <- generate_continuous()
                 X <- dat$data
                 Z <- dat$covts

                 # fit both models
                 out_grid <- tryCatch(covdepGE(X, Z, ssq = hp_grid$ssq,
                                               sbsq = hp_grid$sbsq,
                                               pip = hp_grid$pip, prog_bar = F),
                                      error = function(msg) msg)
                 out_impt <- tryCatch(covdepGE(X, Z, ssq = hp_grid$ssq,
                                               sbsq = hp_grid$sbsq,
                                               pip = hp_grid$pip, grid_search = F, prog_bar = F),
                                      error = function(msg) msg)
                 cat("...", trial_ind, sep = "")
                 # return each of the resulting models and the data
                 list(data = dat, grid = out_grid, impt = out_impt)
               }

# save res
#save(res, file = "opt_pi.Rda")
# save(res, file = "impt_hybrid.Rda")
save(res, file = "impt_hybrid2.Rda")
Sys.time() - start
stopImplicitCluster()
