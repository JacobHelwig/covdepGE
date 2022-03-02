set.seed(1)
n <- 10; p <- 5
X <- matrix(rnorm(n * p), n, p)
W <- matrix(runif(n * n), n, n)

# will result in n * p sigmas
# the l,j entry of LS_sigma is the residual error for a regression on the j-th
# variable weighted with respect to the l-th individual
j <- 1

LS_sigma2 <- rep(NA, n)

# fix the j-th variable as the response
y <- X[ , j]
X_j <- X[ , -j]

for (l in 1:n){

  # perform the weighted regression with respect to individual l
  w <- W[ , l]

  lm_lj <- lm(sqrt(w) * y~0 + I(cbind(1, X_j) * sqrt(w)))
  resid <- lm_lj$residuals / sqrt(w)
  rss <- sum(w * resid^2)
  resvar <- rss / lm_lj$df.residual
  LS_sigma2[l] <- resvar

  # same as:
  # lm_lj <- lm(y~X_j, weights = w)
  # LS_sigma2[l] <- summary(lm_lj)$sigma^2
}

LS_sigma2

# the m, l entry of LS_sig2_l is the least squares estimate to the error term
# variance for indvidual m in a regression weighted with respect to indvidual l
LS_sig2_l <- matrix(LS_sigma2, n, n, T) / W
