## _____________________________________________________________________________
## _____________________________silverman_______________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## Function to calculate silverman's rule of thumb for a data vector according
## to: https://github.com/statsmodels/statsmodels/blob/main/statsmodels/nonparametric/bandwidths.py
## -----------------------------ARGUMENTS---------------------------------------
## x: n x 1 vector; data vector for which the bandwidth will be estimated
## -----------------------------RETURNS-----------------------------------------
## returns scalar bandwidth estimate
silverman <- function(x){

  # apply and return silverman's rule of thumb
  sigma <- (0.9 * min(sd(x), IQR(x) / 1.35) * length(x)^(-0.2))
  return(sigma)
}

## _____________________________________________________________________________
## _____________________________phi0_k.z________________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## Function to calulate the density of a point under a mixture density of n
## Gaussians
## -----------------------------ARGUMENTS---------------------------------------
## z: scalar; point for which the density will be calculated
## mu: n x 1 vector; j-th entry is the mean of the j-th Gaussian
## sigma: scalar; common standard deviation of the Gaussians
## -----------------------------RETURNS-----------------------------------------
## returns scalar density
phi0_k.z <- function(z, mu, sigma){

  # calculate and return the density of z
  return((1 / length(mu)) * sum(dnorm(z, mu, sigma)))
}

## _____________________________________________________________________________
## _____________________________phi0_k__________________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## Function to calulate the densities of a vector of points under a mixture
## density of n Gaussians
## -----------------------------ARGUMENTS---------------------------------------
## Z: N x 1 vector; points for which the densities will be calculated
## mu: n x 1 vector; j-th entry is the mean of the j-th Gaussian
## sigma: scalar; common standard deviation of the Gaussians
## -----------------------------RETURNS-----------------------------------------
## returns N x 1 vector of densities
phi0_k <- function(Z, mu, sigma){
  return(sapply(Z, phi0_k.z, mu = mu, sigma = sigma))
}

## _____________________________________________________________________________
## _____________________________get_bandwidths__________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## Function to calculate individual specific bandwidths for data following the
## methodology described in "A Two-Step Geometric Framework For Density Modeling"
## -----------------------------ARGUMENTS---------------------------------------
## X: n x p matrix; data
## -----------------------------RETURNS-----------------------------------------
## returns n x 1 vector of bandwidths
get_bandwidths <- function(X){

  # get dimensions of X
  n <- nrow(X)
  p <- ncol(X)

  # find the component-wise density for each of the individuals
  # also find silverman's rule of thumb for each of the columns of X
  densities <- matrix(NA, n, p)
  sigma <- rep(NA, p)
  for (j in 1:p){

    # find the value of sigma corresponding to the j-th predictor
    sigma[j] <- silverman(X[ , j])

    # calculate the resulting densities
    densities[ , j] <- phi0_k(X[ , j], X[ , j], sigma[j])
  }

  # calculate the harmonic mean for the sigma^2
  H <- 1 / mean(1 / sigma)

  # calculate the square root of the row-wise product of the densities
  rowProds_sqrt <- rep(NA, n)
  for (l in 1:n){
    rowProds_sqrt[l] <- prod(sqrt(densities[l, ]))
  }

  # return the final bandwidths
  return(H / rowProds_sqrt)
}

## _____________________________________________________________________________
## _____________________________get_weights_____________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## Function to calculate weight matrix
## -----------------------------ARGUMENTS---------------------------------------
## Z: n by p' matrix; extraneous covariates
## norm: scalar in [1, Inf]; norm to use when calculating weights
## kde: boolean; if T, use 2-step KDE methodology described in (2) to calculate
## individual-specific bandwidths
## tau: n x 1 vector, entries in (0, Inf); bandwidth parameters
## -----------------------------RETURNS-----------------------------------------
## D: n x n matrix of weights; j, i entry is the weighting of the j-th
## individual with respect to the i-th individual using the i-th individual's
## bandwidth
## bandwidths: n x 1 vector; individual-specific bandwidths
get_weights <- function(Z, norm, kde, tau){

  # get n
  n <- nrow(Z)

  # if kde, get individual-specific bandwidths
  if (kde){
    tau <- get_bandwidths(Z)
  }else if (length(tau) == 1){
    tau <- rep(tau, n)
  }

  D <- matrix(NA, n, n)

  # calculate the weighting of the j-th individual with respect to the i-th
  # individual using the i-th individual's bandwidth
  for (i in 1:n) {

    # if kde = F, then weights are symmetric (until column scaling);
    # skip repeating unneccesary calculations
    j_inds <- 1:n
    if (!kde) j_inds <- i:n

    # fix the i-th individual's bandwidth
    tau_i <- tau[i]

    for (j in j_inds) {

      # take the p-norm
      diff_vec <- Z[i, ] - Z[j, ]

      if (norm == 2){

        # take the 2-norm, use crossprod
        diff_norm <- sqrt(as.numeric(crossprod(diff_vec)))
      }else if (is.infinite(norm)){

        # take the infinity norm
        diff_norm <- max(abs(diff_vec))
      }else{

        # take the p-norm
        diff_norm <- (sum(abs(diff_vec)^norm))^(1 / norm)
      }


      # given the norm, find the weight using the normal density
      D[j, i] <- stats::dnorm(diff_norm, 0, tau_i)

      # if kde = F, then the weight matrix is symmetric up until scaling
      if (!kde) D[i, j] <- D[j, i]
    }
  }

  # Scale weights to sum to n down the columns
  D <- n * (D) * matrix(1 / colSums(D), n, n, T)

  return(list(D = D, bandwidths = tau))
}
