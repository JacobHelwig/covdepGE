#-------------------------------------------------------------------------------
#-------------------FUNCTION TO GENERATE CONTINUOUS DATA------------------------
#-------------------------------------------------------------------------------

# function to generate the continuous data and covariates
# takes a RNG seed and sample size for each interval
# returns the continuous data, covariates, and true precision matrices
generate_continuous <- function(seed = 1, n1 = 60, n2 = 60, n3 = 60, p = 4){

  set.seed(seed)

  # create covariate for individuals in each of the three intervals

  # define the dimensions of the data
  n <- sum(n1, n2, n3)

  # define the limits of the intervals
  limits1 <- c(-.990, -.331)
  limits2 <- c(-.229, 0.329)
  limits3 <- c(0.431, 0.990)

  # define the covariate values within each interval
  z1 <- seq(limits1[1], limits1[2], length = n1)
  z2 <- seq(limits2[1], limits2[2], length = n2)
  z3 <- seq(limits3[1], limits3[2], length = n3)
  Z <- matrix(c(z1, z2, z3), n, 1)

  # create precision matrices

  # the shared part of the structure for all three intervals is a 2 on the
  # diagonal and a 1 in the (2, 3) position
  common_str <- diag(p + 1)
  common_str[2, 3] <- 1

  # define constants for the structure of interval 2
  const1 <- 0.23
  const2 <- 0.56

  # interval 2 has two different linear functions of Z in the (1, 2) position
  # and (1, 3) positions; define structures for each of these components
  int2_str12 <- int2_str13 <- matrix(0, p + 1, p + 1)
  int2_str12[1, 2] <- int2_str13[1, 3] <- 1

  # define the precision matrices for each of the individuals in interval 2
  int2_prec <- lapply(z2, function(z) common_str +
                        ((1 - (z + const1) / const2) * int2_str12) +
                        ((z + const1) / const2 * int2_str13))

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


#-------------------------------------------------------------------------------
#-------------------FUNCTION TO GENERATE DISCRETE DATA--------------------------
#-------------------------------------------------------------------------------

# function to generate the discrete data and covariates
# takes a RNG seed, sample size, number of predictors, lambda (root of the
# non-zero elements in the precision matrix), values to populate the
# extraneous covariate vector with, and boolean for whether the two groups
# should have the same covariance
# returns the continuous data, covariates, and true covariance matrix
generate_discrete <- function(seed = 1, n = 100, p = 10, lambda = 15,
                              cov1 = -0.1, cov2 = 0.1, same = T){

  set.seed(seed)

  # generating the precision matrix: Assume two discrete covariate levels, one
  # for each group
  Lam1 <- c(rep(lambda, 4), rep(0, p - 3))
  Lam2 <- c(rep(0, 4), rep(lambda, p - 3))

  # if same is true, the individuals in both groups will have the same
  # covariance matrix
  if (same) Lam2 <- Lam1

  # create covariance matrix for both groups
  Var1 <- solve(Lam1 %*% t(Lam1) + diag(rep(10, p + 1)))
  Var2 <- solve(Lam2 %*% t(Lam2) + diag(rep(10, p + 1)))

  # create the extraneous covariate; individuals in group j have a covariate
  # vector of length p with cov_j as the only entry, j\in {1,2}
  Z <- matrix(c(rep(cov1, n %/% 2), rep(cov2, n %/% 2)), n, p)

  # create the data matrix; individuals in group j are generated from a MVN with
  # 0 mean vector and covariance matrix Var_j, j\in {1,2}
  X1 <- MASS::mvrnorm(n %/% 2, rep(0, p + 1), Var1)
  X2 <- MASS::mvrnorm(n %/% 2, rep(0, p + 1), Var2)
  data_mat <- rbind(X1, X2)

  return(list(data = data_mat, covts = Z, true_covariance = list(Var1, Var2)))

}
