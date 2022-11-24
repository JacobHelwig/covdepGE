# function for generating discrete covariate dependent data
disc_cov_dep_data <- function(p, n1, n2, lambda){

  # create covariate for observations

  # define number of samples and covariate
  n <- sum(n1, n2)
  z1 <- rep(1, n1)
  z2 <- rep(2, n2)
  Z <- matrix(c(z1, z2), n, 1)

  # create the precision matrices
  lambda1 <- c(rep(lambda, 4), rep(0, p - 4))
  lambda2 <- c(rep(0, p - 4), rep(lambda, 4))
  prec1 <- tcrossprod(lambda1) + 10 * diag(p)
  prec2 <- tcrossprod(lambda2) + 10 * diag(p)
  true_precision <- c(replicate(n1, prec1, simplify = F),
                      replicate(n2, prec2, simplify = F))

  # invert the precision matrices to get the covariance matrices
  cov_mats <- lapply(true_precision, solve)

  # generate the data using the covariance matrices
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p)))

  return(list(X = data_mat, Z = Z, true_precision = true_precision))
}

# function for generating continuous covariate dependent data
cont_cov_dep_data <- function(p, n1, n2, n3){

  # create covariate for observations in each of the three intervals

  # define number of samples
  n <- sum(n1, n2, n3)

  # define the intervals
  limits1 <- c(-3, -1)
  limits2 <- c(-1, 1)
  limits3 <- c(1, 3)

  # define the interval labels
  interval <- c(rep(1, n1), rep(2, n2), rep(3, n3))

  # draw the covariate values within each interval
  z1 <- sort(stats::runif(n1, limits1[1], limits1[2]))
  z2 <- sort(stats::runif(n2, limits2[1], limits2[2]))
  z3 <- sort(stats::runif(n3, limits3[1], limits3[2]))
  Z <- matrix(c(z1, z2, z3), n, 1)

  # the shared part of the structure for all three intervals is a 2 on the
  # diagonal and a 1 in the (2, 3) position
  common_str <- diag(p)
  common_str[2, 3] <- 1

  # define constants for the structure of interval 2
  beta1 <- diff(limits2)^-1
  beta0 <- -limits2[1] * beta1

  # define omega12 and omega 13
  omega12 <- (Z < 1) * pmin(1, 1 - beta0 - beta1 * Z)
  omega13 <- (Z > -1) * pmin(1, beta0 + beta1 * Z)

  # interval 2 has two different linear functions of Z in the (1, 2) position
  # and (1, 3) positions; define structures for each of these components
  str12 <- str13 <- matrix(0, p, p)
  str12[1, 2] <- str13[1, 3] <- 1

  # create the precision matrices
  prec_mats <- vector("list", n)
  for (j in 1:n){
    prec_mats[[j]] <- common_str + omega12[j] * str12 + omega13[j] * str13
  }

  # symmetrize the precision matrices
  true_precision <- lapply(prec_mats, function(mat) t(mat) + mat)

  # invert the precision matrices to get the covariance matrices
  cov_mats <- lapply(true_precision, solve)

  # generate the data using the covariance matrices
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p)))

  return(list(X = data_mat, Z = Z, true_precision = true_precision,
              interval = interval))
}
