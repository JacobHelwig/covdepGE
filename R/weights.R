## -----------------------------------------------------------------------------
## -----------------------------silverman---------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## calculate silverman's rule of thumb for a data vector according to:
## https://github.com/statsmodels/statsmodels/blob/main/statsmodels/nonparametric/bandwidths.py
## -----------------------------ARGUMENTS---------------------------------------
## x: vector of length n; data vector for which the bandwidth will be estimated
## -----------------------------RETURNS-----------------------------------------
## returns bandwidth estimate
## -----------------------------------------------------------------------------
silverman <- function(x) {

  # apply and return silverman's rule of thumb
  sigma <- (0.9 * min(stats::sd(x), stats::IQR(x) / 1.35) * length(x)^(-0.2))

  # ensure that sigma is greater than 0
  if (sigma == 0) {
    warning("Cannot calculate Silverman's rule of thumb for constant Z; returning 1")
    sigma <- 1
  }

  return(sigma)
}

## -----------------------------------------------------------------------------
## -----------------------------phi0_k.z----------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## calulate the density of a point under a mixture of n Gaussians with common
## standard deviation
## -----------------------------ARGUMENTS---------------------------------------
## z: numeric; point for which the density will be calculated
##
## mu: numeric vector of length n; j-th entry is the mean of the j-th Gaussian
##
## sigma: positive numeric; common standard deviation of the Gaussians
## -----------------------------RETURNS-----------------------------------------
## returns density
## -----------------------------------------------------------------------------
phi0_k.z <- function(z, mu, sigma) {

  # calculate and return the density of z
  return((1 / length(mu)) * sum(stats::dnorm(z, mu, sigma)))
}

## -----------------------------------------------------------------------------
## -----------------------------phi0_k------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## calulate the densities of a vector of points under a mixture
## of n Gaussians with common standard deviation
## -----------------------------ARGUMENTS---------------------------------------
## Z: numeric vector of length N; points for which the densities will be
## calculated
##
## mu: numeric vector of length n; j-th entry is the mean of the j-th Gaussian
##
## sigma: positive numeric; common standard deviation of the Gaussians
## -----------------------------RETURNS-----------------------------------------
## returns numeric vector of N densities
## -----------------------------------------------------------------------------
phi0_k <- function(Z, mu, sigma) {
  return(sapply(Z, phi0_k.z, mu = mu, sigma = sigma))
}

## -----------------------------------------------------------------------------
## -----------------------------get_bandwidths----------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## calculate observation specific bandwidths for data following the
## methodology described in "A Two-Step Geometric Framework For Density Modeling"
## -----------------------------ARGUMENTS---------------------------------------
## X: n x p matrix; data
## -----------------------------RETURNS-----------------------------------------
## returns vector of n bandwidths
## -----------------------------------------------------------------------------
get_bandwidths <- function(X) {

  # get dimensions of X
  n <- nrow(X)
  p <- ncol(X)

  # find the component-wise density for each of the observations
  # also find silverman's rule of thumb for each of the columns of X
  densities <- matrix(NA, n, p)
  sigma <- rep(NA, p)
  for (j in 1:p) {

    # find the value of sigma corresponding to the j-th predictor
    sigma[j] <- silverman(X[, j])

    # calculate the resulting densities
    densities[, j] <- phi0_k(X[, j], X[, j], sigma[j])
  }

  # calculate the harmonic mean for the sigma^2
  H <- 1 / mean(1 / sigma)

  # calculate the square root of the row-wise product of the densities
  rowProds_sqrt <- rep(NA, n)
  for (l in 1:n) {
    rowProds_sqrt[l] <- prod(sqrt(densities[l, ]))
  }

  # return the final bandwidths
  return(H / rowProds_sqrt)
}

## -----------------------------------------------------------------------------
## -----------------------------get_weights-------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
## function to calculate weight matrix
## -----------------------------ARGUMENTS---------------------------------------
## Z: n by q matrix; extraneous covariates
##
## norm: numeric in [1, Inf]; norm to use when calculating weights
##
## tau: NULL OR vector of length n with positive entries; bandwidth parameters.
## If NULL, use KDE to get bandwidths
## -----------------------------RETURNS-----------------------------------------
## D: n x n matrix of weights; i, j entry is the weighting of the i-th
## observation with respect to the j-th observation using the j-th observation's
## bandwidth
##
## bandwidths: vector of length n; observation-specific bandwidths
## -----------------------------------------------------------------------------
get_weights <- function(Z, norm, tau) {

  # get n
  n <- nrow(Z)

  # if tau is null, use kde to get observation-specific bandwidths
  if (is.null(tau)) {
    tau <- get_bandwidths(Z)
  } else if (length(tau) == 1) {
    tau <- rep(tau, n)
  }

  D <- matrix(NA, n, n)

  # calculate the weighting of the j-th observation with respect to the i-th
  # observation using the i-th observation's bandwidth
  for (i in 1:n) {
    for (j in i:n) {

      # take the p-norm
      diff_vec <- Z[i, ] - Z[j, ]

      if (norm == 2) {

        # take the 2-norm, use crossprod
        diff_norm <- sqrt(as.numeric(crossprod(diff_vec)))
      } else if (is.infinite(norm)) {

        # take the infinity norm
        diff_norm <- max(abs(diff_vec))
      } else {

        # take the p-norm
        diff_norm <- (sum(abs(diff_vec)^norm))^(1 / norm)
      }


      # given the norm, find the weight of j with respect to i
      D[j, i] <- stats::dnorm(diff_norm, 0, tau[i])

      # find the weight of i with respect to j
      D[i, j] <- stats::dnorm(diff_norm, 0, tau[j])
    }
  }

  # Scale weights to sum to n down the columns
  D <- n * (D) * matrix(1 / colSums(D), n, n, TRUE)

  return(list(D = D, bandwidths = tau))
}
