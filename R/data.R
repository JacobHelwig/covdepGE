## -----------------------------------------------------------------------------
#' @title Generate Covariate-Dependent Data
#' @export
## -----------------------------DESCRIPTION-------------------------------------
#' @description function to generate a `1`-dimensional extraneous covariate and
#' `p`-dimensional Gaussian data with a precision `matrix` that varies as a
#' continuous function of the extraneous covariate
## -----------------------------ARGUMENTS---------------------------------------
#' @param p positive `integer`; number of variables in the data `matrix`
#'
#' @param n1 positive `integer`; number of observations in the first interval
#'
#' @param n2 positive `integer`; number of observations in the second interval
#'
#' @param n3 positive `integer`; number of observations in the third interval
## -----------------------------RETURNS-----------------------------------------
#' @return Returns list with the following values:
#'
#'  \enumerate{
#'    \item `data`: a `(n1 + n2 + n3)` `x` `p` `numeric` `matrix`, where the
#'    `i`-th row is drawn from a `p`-dimensional Gaussian with mean `0` and
#'    precision `matrix` `true_precision[[i]]`
#'
#'    \item `covts`: a `(n1 + n2 + n3)` `x` `1` `numeric` `matrix`, where the
#'    `i`-th entry is the extraneous covariate `z_i` for observation `i`
#'
#'    The first `n1` observations have `z_i` from from a uniform distribution on
#'    the interval `(-3,-1)` (the first interval)
#'
#'    Observations `n1 + 1` to `n1 + n2` have `z_i` from from a uniform
#'    distribution on the interval `(-1,1)` (the second interval)
#'
#'    observations `n1 + n2 + 1` to `n1 + n2 + n3` have `z_i` from a uniform
#'    distribution on the interval `(1,3)` (the third interval)
#'
#'    \item `true_precision`: `list` of `n1 + n2 + n3` `p` `x` `p` matrices; the
#'    `i`-th `matrix` is the precision `matrix` for the `i`-th observation
#'
#'    All precision matrices have `2` on the diagonal and `1` in the
#'    `(2,3)`/`(3,2)` positions
#'
#'    Observations in the first interval have a `1` in the `(1,2)`/`(2,1)`
#'    positions, while observations in the third interval have a `1` in the
#'    `(1,3)`/`(3,1)` positions
#'
#'    Observations in the second interval (`(-1,1)`) have 2 entries that vary as
#'    a linear function of their extraneous covariate. Let `beta = 1/2`. Then,
#'    the `(1, 2)`/`(2, 1)` positions for the `i`-th observation in interval 2
#'    are `beta * (1 - z_i)`, while the `(1, 3)`/`(3, 1)` entries are
#'    `beta * (1 + z_i)`
#'
#'    Thus, as `z_i` approaches `-1` from the right, the associated precision
#'    `matrix` becomes more similar to the `matrix` for observations in the
#'    first interval. Similarly, as `z_i` approaches `1` from the left, the
#'    `matrix` becomes more similar to the `matrix` for observations in the
#'    third interval
#'
#'    \item `interval`: `vector` of length `n1 + n2 + n3`; contains the ground
#'    truth interval assignments for each of the observations
#' }
## -----------------------------EXAMPLES----------------------------------------
#' @examples
#'
#' library(ggplot2)
#'
#' # get the data
#' set.seed(1)
#' data <- generateData()
#' X <- data$data
#' Z <- data$covts
#' interval <- data$interval
#' prec <- data$true_precision
#'
#' # get overall and within interval sample sizes
#' n <- nrow(X)
#' n1 <- sum(interval == 1)
#' n2 <- sum(interval == 2)
#'
#' # visualize the distribution of the extraneous covariate
#' ggplot(data.frame(Z = Z, interval = as.factor(interval))) +
#' geom_histogram(aes(Z, fill = interval), color = "black", bins = n %/% 5)
#'
#' # visualize the true precision matrices in each of the intervals
#'
#' # interval 1
#' matViz(prec[[1]], incl_val = TRUE) +
#' ggtitle("True precision matrix, interval 1")
#'
#' # interval 2 (varies continuously with Z)
#' int2_mats <- prec[interval == 2]
#' int2_inds <- c(5, n2 %/% 2, n2 - 5)
#' lapply(int2_inds, function(j) matViz(int2_mats[[j]], incl_val = TRUE) +
#' ggtitle(paste("True precision matrix, interval 2, observation", j)))
#'
#' # interval 3
#' matViz(prec[[length(prec)]], incl_val = TRUE) +
#' ggtitle("True precision matrix, interval 3")
#'
#' # fit the model and visualize the estimated precision matrices
#' (out <- covdepGE(X, Z))
#' plot(out)
#'
#' # visualize the inclusion probabilities for variables (1, 3) and (1, 2)
#' inclusionCurve(out, 1, 2)
#' inclusionCurve(out, 1, 3)
## -----------------------------------------------------------------------------
generateData <- function(p = 5, n1 = 60, n2 = 60, n3 = 60){

  # create covariate for observations in each of the three intervals

  # define number of samples
  n <- sum(n1, n2, n3)

  # define the intervals and assignments
  limits1 <- c(-3, -1)
  limits2 <- c(-1, 1)
  limits3 <- c(1, 3)
  interval <- c(rep(1, n1), rep(2, n2), rep(3, n3))

  # draw the covariate values within each interval
  z1 <- sort(stats::runif(n1, limits1[1], limits1[2]))
  z2 <- sort(stats::runif(n2, limits2[1], limits2[2]))
  z3 <- sort(stats::runif(n3, limits3[1], limits3[2]))
  Z <- `matrix`(c(z1, z2, z3), n, 1)

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
  int2_str12 <- int2_str13 <- `matrix`(0, p, p)
  int2_str12[1, 2] <- int2_str13[1, 3] <- 1

  # define the precision matrices for each of the observations in interval 2
  int2_prec <- lapply(z2, function(z) common_str +
                        ((1 - beta0 - beta1*z)*int2_str12) +
                        ((beta0 + beta1*z)*int2_str13))

  # interval 1 has a 1 in the (1, 2) position and interval 3 has a 1 in the
  # (1, 3) position; define structures for each of these components
  int1_str12 <- int3_str13 <- matrix(0, p, p)
  int1_str12[1, 2] <- int3_str13[1, 3] <- 1

  # define the precision matrices for each of the observations in interval 1 and
  # interval 3
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
