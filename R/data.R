## -----------------------------------------------------------------------------
#' @title generate_continuous
#' @export
## -----------------------------DESCRIPTION-------------------------------------
#' @description function to generate a `1`-dimensional extraneous covariate and
#' `p`-dimensional Gaussian data with a precision matrix having entries that
#' vary as a continuous function of the extraneous covariate
## -----------------------------ARGUMENTS---------------------------------------
#' @param p positive integer; number of variables in the data matrix
#'
#' @param n1 positive integer; number of individuals in the first interval
#'
#' @param n2 positive integer; number of individuals in the second interval
#'
#' @param n3 positive integer; number of individuals in the third interval
## -----------------------------RETURNS-----------------------------------------
#' @return Returns list with the following values:
#'
#' 1. `data`: a `(n1 + n2 + n3) x p` numeric matrix, where the `i`-th row is
#' drawn from a `p`-dimensional Gaussian with mean `0` and precision matrix
#' `true_precision`
#'
#' 2. `covts`: a `(n1 + n2 + n3) x 1` numeric matrix, where the `i`-th entry is
#' the extraneous covariate `z_i` for individual `i`.
#'
#' The first `n1` individuals have `z_i` from from a uniform distribution on the
#' interval `(-3, -1)`.
#'
#' Individuals `n1 + 1` to `n1 + n2` have `z_i` from from a uniform distribution
#' on the interval `(-1, 1)`.
#'
#' Individuals `n1 + n2 + 1` to `n1 + n2 + n3` have `z_i` from a uniform
#' distribution on the interval `(1, 3)`
#'
#' 3. `true_precision`: `list` of `n1 + n2 + n3 p x p` matrices; the `i`-th
#' `matrix` is the precision matrix for the `i`-th individual
#'
#' All precision matrices have `2` on the diagonal and `1` in the
#' `(2, 3)/(3, 2)` position
#'
#' Individuals in the first interval (`(-3, -1)`) have a `1` in the
#' `(1, 2)/(2, 1)` position, while individuals in the third interval (`(1, 3)`)
#' have a `1` in the `(1, 3)/(3, 1)` position.
#'
#' Individuals in the second interval (`(-1, 1)`) have 2 entries that vary as
#' a linear function of their extraneous covariate. Let `beta = 1/2`. Then,
#' the `(1, 2)/(2, 1)` position for the `i`-th individual in interval 2 is equal
#' to `beta * (1 - z_i)`, while the `(1, 3)/(3, 1)` entry is equal to
#' `beta * (1 + z)`.
#'
#' Thus, as `z_i` approaches `-1` from the right, the precision matrix becomes
#' more similar to the matrix for individuals in interval 1. Similarly, as `z_i`
#' approaches `1` from the left, the matrix becomes more similar to the matrix
#' for individuals in interval `3`
#'
#' 4. `interval`: `vector` of length `n1 + n2 + n3`; contains the ground truth
#' interval assignments for each of the individuals
## -----------------------------EXAMPLES----------------------------------------
#' @examples
## -----------------------------------------------------------------------------
generate_continuous <- function(p = 5, n1 = 60, n2 = 60, n3 = 60){

  # create covariate for individuals in each of the three intervals

  # define the dimensions of the data
  n <- sum(n1, n2, n3)

  # define the intervals and assignments
  limits1 <- c(-3, -1)
  limits2 <- c(-1, 1)
  limits3 <- c(1, 3)
  interval <- c(rep(1, n1), rep(2, n2), rep(3, n3))

  # define the covariate values within each interval
  z1 <- sort(stats::runif(n1, limits1[1], limits1[2]))
  z2 <- sort(stats::runif(n2, limits2[1], limits2[2]))
  z3 <- sort(stats::runif(n3, limits3[1], limits3[2]))
  Z <- matrix(c(z1, z2, z3), n, 1)

  # create precision matrices

  # the shared part of the structure for all three intervals is a 2 on the
  # diagonal and a 1 in the (2, 3) position
  common_str <- diag(p)
  common_str[2, 3] <- 1

  # define constants for the structure of interval 2
  beta1 <- diff(limits2)^-1
  beta0 <- -limits2[1] * beta1

  # interval 2 has two different linear functions of Z in the (1, 2) position
  # and (1, 3) positions; define structures for each of these components
  int2_str12 <- int2_str13 <- matrix(0, p, p)
  int2_str12[1, 2] <- int2_str13[1, 3] <- 1

  # define the precision matrices for each of the individuals in interval 2
  int2_prec <- lapply(z2, function(z) common_str +
                        ((1 - beta0 - beta1*z)*int2_str12) +
                        ((beta0 + beta1*z)*int2_str13))

  # interval 1 has a 1 in the (1, 2) and interval 3 has a 1 in the (1, 3) position;
  # define structures for each of these components
  int1_str12 <- int3_str13 <- matrix(0, p, p)
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
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p)))

  return(list(data = data_mat, covts = Z, true_precision = prec_mats,
              interval = interval))
}
