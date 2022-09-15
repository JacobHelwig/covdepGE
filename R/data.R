## -----------------------------------------------------------------------------
## Distributed under GPL (≥ 3) license
#'
#' @title Generate Covariate-Dependent Data
#' @export
## -----------------------------DESCRIPTION-------------------------------------
#' @description Generate a \eqn{1}-dimensional extraneous covariate
#' and \eqn{p}-dimensional Gaussian data with a precision matrix that varies as
#' a continuous function of the extraneous covariate. This data is distributed
#' similar to that used in the simulation study from (1)
## -----------------------------ARGUMENTS---------------------------------------
#' @param p positive integer; number of variables in the data matrix. `5` by
#' default
#'
#' @param n1 positive integer; number of observations in the first interval.
#' `60` by default
#'
#' @param n2 positive integer; number of observations in the second interval.
#' `60` by default
#'
#' @param n3 positive integer; number of observations in the third interval.
#' `60` by default
#'
#' @param Z `NULL` or numeric vector; extraneous covariate values for each
#' observation. If `NULL`, `Z` will be generated from a uniform distribution on
#' each of the intervals
#'
#' @param true_precision `NULL` OR list of matrices of dimension
#' \eqn{p \times p}{p x p}; true precision matrix for each observation. If
#' `NULL`, the true precision matrices will be generated dependent on `Z`.
#' `NULL` by default
## -----------------------------RETURNS-----------------------------------------
#' @return Returns list with the following values:
#'
#'    \item{X}{a `(n1 + n2 + n3)` \eqn{\times p}{x p} numeric matrix, where
#'    the \eqn{i}-th row is drawn from a \eqn{p}-dimensional Gaussian with mean
#'    \eqn{0} and precision matrix `true_precision[[i]]`}
#'
#'    \item{Z}{a `(n1 + n2 + n3)` \eqn{\times 1}{x 1} numeric matrix, where
#'    the \eqn{i}-th entry is the extraneous covariate \eqn{z_i}{zi} for
#'    observation \eqn{i}}
#'
#'    \item{true_precision}{list of `n1 + n2 + n3` matrices of dimension
#'    \eqn{p \times p}{p x p}; the \eqn{i}-th matrix is the precision matrix for
#'    the \eqn{i}-th observation}
#'
#'    \item{interval}{vector of length `n1 + n2 + n3`; interval assignments
#'    for each of the observations, where the \eqn{i}-th entry is the interval
#'    assignment for the \eqn{i}-th observation}
## -----------------------------DETAILS-----------------------------------------
#' @details
#' # Extraneous Covariate
#'
#' If `Z = NULL`, then the generation of `Z` is as follows:
#'
#' The first `n1` observations have \eqn{z_i}{zi} from from a uniform
#' distribution on the interval \eqn{(-3, -1)} (the first interval).
#'
#' Observations `n1 + 1` to `n1 + n2` have \eqn{z_i}{zi} from from a uniform
#' distribution on the interval \eqn{(-1, 1)} (the second interval).
#'
#' Observations `n1 + n2 + 1` to `n1 + n2 + n3` have \eqn{z_i}{zi} from a
#' uniform distribution on the interval \eqn{(1, 3)} (the third interval).
#'
#' # Precision Matrices
#'
#' If `true_precision = NULL`, then the generation of the true precision
#' matrices is as follows:
#'
#' All precision matrices have \eqn{2} on the diagonal and \eqn{1} in the
#' \eqn{(2, 3)/ (3, 2)} positions.
#'
#' Observations in the first interval have a \eqn{1} in the
#' \eqn{(1, 2) / (1, 2)} positions, while observations in the third interval
#' have a \eqn{1} in the \eqn{(1, 3)/ (3, 1)} positions.
#'
#' Observations in the second interval have \eqn{2} entries that vary as a
#' linear function of their extraneous covariate. Let
#' \eqn{\beta = 1/2}{beta = 1/2}. Then, the \eqn{(1, 2)/(2, 1)} positions for
#' the \eqn{i}-th observation in the second interval are
#' \eqn{\beta\cdot(1 - z_i)}{beta(1 - zi)}, while the \eqn{(1, 3)/ (3, 1)}
#' entries are \eqn{\beta\cdot(1 + z_i)}{beta(1 + zi)}.
#'
#' Thus, as \eqn{z_i}{zi} approaches \eqn{-1} from the right, the associated
#' precision matrix becomes more similar to the matrix for observations in the
#' first interval. Similarly, as \eqn{z_i}{zi} approaches \eqn{1} from the left,
#' the matrix becomes more similar to the matrix for observations in the third
#' interval.
## -----------------------------EXAMPLES----------------------------------------
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # get the data
#' set.seed(12)
#' data <- generateData()
#' X <- data$X
#' Z <- data$Z
#' interval <- data$interval
#' prec <- data$true_precision
#'
#' # get overall and within interval sample sizes
#' n <- nrow(X)
#' n1 <- sum(interval == 1)
#' n2 <- sum(interval == 2)
#' n3 <- sum(interval == 3)
#'
#' # visualize the distribution of the extraneous covariate
#' ggplot(data.frame(Z = Z, interval = as.factor(interval))) +
#'   geom_histogram(aes(Z, fill = interval), color = "black", bins = n %/% 5)
#'
#' # visualize the true precision matrices in each of the intervals
#'
#' # interval 1
#' matViz(prec[[1]], incl_val = TRUE) +
#'   ggtitle(paste0("True precision matrix, interval 1, observations 1,...,", n1))
#'
#' # interval 2 (varies continuously with Z)
#' cat("\nInterval 2, observations ", n1 + 1, ",...,", n1 + n2, sep = "")
#' int2_mats <- prec[interval == 2]
#' int2_inds <- c(5, n2 %/% 2, n2 - 5)
#' lapply(int2_inds, function(j) matViz(int2_mats[[j]], incl_val = TRUE) +
#'          ggtitle(paste("True precision matrix, interval 2, observation", j + n1)))
#'
#' # interval 3
#' matViz(prec[[length(prec)]], incl_val = TRUE) +
#'   ggtitle(paste0("True precision matrix, interval 3, observations ",
#'                  n1 + n2 + 1, ",...,", n1 + n2 + n3))
#'
#' # fit the model and visualize the estimated graphs
#' (out <- covdepGE(X, Z))
#' plot(out)
#'
#' # visualize the posterior inclusion probabilities for variables (1, 3) and (1, 2)
#' inclusionCurve(out, 1, 2)
#' inclusionCurve(out, 1, 3)
#' }
## -----------------------------------------------------------------------------
generateData <- function(p = 5, n1 = 60, n2 = 60, n3 = 60, Z = NULL,
                         true_precision = NULL){

  # create covariate for observations in each of the three intervals

  # define number of samples
  n <- ifelse(is.null(true_precision), sum(n1, n2, n3), length(true_precision))

  # define the intervals
  limits1 <- c(-3, -1)
  limits2 <- c(-1, 1)
  limits3 <- c(1, 3)

  # if Z and true_precision have not been provided, generate Z
  interval <- NULL
  if (is.null(Z) & is.null(true_precision)){

    # define the interval labels
    interval <- c(rep(1, n1), rep(2, n2), rep(3, n3))

    # draw the covariate values within each interval
    z1 <- sort(stats::runif(n1, limits1[1], limits1[2]))
    z2 <- sort(stats::runif(n2, limits2[1], limits2[2]))
    z3 <- sort(stats::runif(n3, limits3[1], limits3[2]))
    Z <- matrix(c(z1, z2, z3), n, 1)
  }else if(!is.null(Z) & is.null(true_precision)){

    # Z has been provided and true_precision has not
    # divide Z into the 3 intervals
    interval <- as.integer(cut(Z, c(-Inf, -1, 1, Inf), labels = 1:3))
    z1 <- Z[interval == 1]
    z2 <- Z[interval == 2]
    z3 <- Z[interval == 3]

    # get the sample size in each of the intervals
    n1 <- length(z1)
    n2 <- length(z2)
    n3 <- length(z3)
  }else if(!is.null(Z) & !is.null(true_precision)){

    # Z and true_precision have been provided
    stop("Z and true_precision cannot both be provided")
  }

  # if they have not been provided, create the precision matrices
  if (is.null(true_precision)){

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

    # define the precision matrices for each of the observations in interval 2
    int2_prec <- lapply(z2, function(z) common_str +
                          ((1 - beta0 - beta1 * z) * int2_str12) +
                          ((beta0 + beta1 * z) * int2_str13))

    # interval 1 has a 1 in the (1, 2) position and interval 3 has a 1 in the
    # (1, 3) position; define structures for each of these components
    int1_str12 <- int3_str13 <- matrix(0, p, p)
    int1_str12[1, 2] <- int3_str13[1, 3] <- 1

    # define the precision matrices for each of the observations in interval 1
    # and interval 3
    int1_prec <- rep(list(common_str + int1_str12), n1)
    int3_prec <- rep(list(common_str + int3_str13), n3)

    # put all of the precision matrices into one list
    prec_mats <- c(int1_prec, int2_prec, int3_prec)

    # symmetrize the precision matrices
    true_precision <- lapply(prec_mats, function(mat) t(mat) + mat)
  }

  # invert the precision matrices to get the covariance matrices
  cov_mats <- lapply(true_precision, solve)

  # generate the data using the covariance matrices
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p)))

  return(list(X = data_mat, Z = Z, true_precision = true_precision,
              interval = interval))
}
