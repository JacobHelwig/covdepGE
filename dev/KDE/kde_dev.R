# takes a n x 1 vector of data
# calculate silverman's rule of thumb for sigma
# https://github.com/statsmodels/statsmodels/blob/main/statsmodels/nonparametric/bandwidths.py
silverman <- function(x){

  # apply and return silverman's rule of thumb
  sigma <- (0.9 * min(sd(x), IQR(x) / 1.35) * length(x)^(-0.2))
  return(sigma)
}

# takes a point z, a vector of n means, and scalar sigma
# returns the density of z under the density of a mixture of Gaussians, each
# centered at mu with standard deviation sigma
phi0_k.z <- function(z, mu, sigma){

  # calculate and return the density of z
  return((1 / length(mu)) * sum(dnorm(z, mu, sigma)))
}

# vectorized version of phi0_k - takes a vector Z of length N, vector of n means
# and a scalar sigma
# Returns a vector of N densities
phi0_k <- function(Z, mu, sigma){
  return(sapply(Z, phi0_k.z, mu = mu, sigma = sigma))
}

# function that takes a datamatrix X and returns a bandwidth for each of the
# subjects (vector of length n)
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

# generate some data and try it out
set.seed(1)
X <- mvtnorm::rmvnorm(100, -2:2)
(bw <- get_bandwidths(X[ , 1:2, drop = F]))
silverman(X[ , 1])
silverman(X[ , 2])
integrate(phi0_k, mu = X[1, ], sigma = bw[1], -Inf, Inf)
library(ggplot2)
ggplot()+ stat_function(fun = function(x) phi0_k(x, mu = X[ , ], sigma = bw[1])) + xlim(-100, 100)

# smaller sample for explicit calculations
set.seed(1)
X <- cbind(runif(3, 2, 5), runif(3, -3, 1))
X
(bw <- get_bandwidths(X))

# manual calculations

# calculate Silverman for both columns
(silv1 <- 0.9 * min(sd(X[ , 1]), IQR(X[ , 1]) / 1.35) * nrow(X)^-0.2)
(silv2 <- 0.9 * min(sd(X[ , 2]), IQR(X[ , 2]) / 1.35) * nrow(X)^-0.2)

# calculate the harmonic mean
(h <- ncol(X) / ((1 / silv1) + (1 / silv2)))

# calculate the square root of the product of the step 1 densities for
# each individual

# individual 1

# sum and normalize the densities for both columns
(dens11 <- 1 / nrow(X) * sum(dnorm(X[1, 1], X[ , 1], silv1)))
(dens12 <- 1 / nrow(X) * sum(dnorm(X[1, 2], X[ , 2], silv2)))

# get the root of the product of the densities
(dens1 <- sqrt(dens11 * dens12))

# individual 2

# sum and normalize the densities for both columns
(dens21 <- 1 / nrow(X) * sum(dnorm(X[2, 1], X[ , 1], silv1)))
(dens22 <- 1 / nrow(X) * sum(dnorm(X[2, 2], X[ , 2], silv2)))

# get the root of the product of the densities
(dens2 <- sqrt(dens21 * dens22))

# individual 3

# sum and normalize the densities for both columns
(dens31 <- 1 / nrow(X) * sum(dnorm(X[3, 1], X[ , 1], silv1)))
(dens32 <- 1 / nrow(X) * sum(dnorm(X[3, 2], X[ , 2], silv2)))

# get the root of the product of the densities
(dens3 <- sqrt(dens31 * dens32))

# final bandwidths:
dens <- c(dens1, dens2, dens3)
(bw.manual <- h / dens)

bw.manual - bw
