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

# src the scripts
if ("covdepGE" %in% .packages()) detach("package:covdepGE", unload = TRUE)
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/covdepGE_main.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/cavi_search.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/weights.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/checks.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/gg_covdepGE.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/covdepGE_R.R")
Rcpp::sourceCpp("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/src/covdepGE_c.cpp")

set.seed(2)

# gen the data
cont <- generate_continuous()
X <- cont$data
Z <- cont$covts

# vanilla
out1 <- covdepGE(X, Z) # hybrid
out2 <- covdepGE(X, Z, hp_method = "grid_search")
out3 <- covdepGE(X, Z, hp_method = "model_average")

# increase the number of hyperparameters
out <- covdepGE(X, Z, nssq = 7, nsbsq = 7, npip = 7)

# increase ssq mult
out <- covdepGE(X, Z, ssq_mult = 10)

# increase ssq lower
out <- covdepGE(X, Z, ssq_lower = 0.1)

# increase snr_upper
out <- covdepGE(X, Z, snr_upper = 25)

# increase sbsq lower
out <- covdepGE(X, Z, sbsq_lower = .1)

# increase pip lower
out <- covdepGE(X, Z, pip_lower = 0.2)

# provide tau
out <- covdepGE(X, Z, tau = 100000)
out <- covdepGE(X, Z, tau = 1e-16)
out <- covdepGE(X, Z, tau = seq(1e-16, 1, length.out = 180))

# change the norm
out <- covdepGE(X, Z, norm = 1)

# do not center the data
out <- covdepGE(X, Z, center_data = F)

# increase the ELBO tol
out <- covdepGE(X, Z, elbo_tol = 1000)

