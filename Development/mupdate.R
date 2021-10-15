
# data generation
rm(list = ls())

n <- 4; p <- 3

(X <- matrix(1:(n * p), n, p)) # data

(D <- matrix(1:(n * n) + max(X), n, n)) # weights

(y <- 1:n + max(D)) # response

(mu <- matrix(1:(n * p) + max(y), n, p)) # mu matrix
#mu <- matrix(0, n, p)
(alpha <- matrix(1:(n * p) + max(mu), n, p)) # alpha matrix
#alpha <- matrix(0, n, p)

# METHOD 1: manual update of mu_1^2 (the 2, 1 entry of the mu matrix)
(D[1, 2] * X[1, 1] * (y[1] - X[1, 2] * mu[2, 2] * alpha[2, 2] - X[1, 3] * mu[2, 3] * alpha[2, 3])) +
  (D[2, 2] * X[2, 1] * (y[2] - X[2, 2] * mu[2, 2] * alpha[2, 2] - X[2, 3] * mu[2, 3] * alpha[2, 3])) +
  (D[3, 2] * X[3, 1] * (y[3] - X[3, 2] * mu[2, 2] * alpha[2, 2] - X[3, 3] * mu[2, 3] * alpha[2, 3])) +
  (D[4, 2] * X[4, 1] * (y[4] - X[4, 2] * mu[2, 2] * alpha[2, 2] - X[4, 3] * mu[2, 3] * alpha[2, 3]))

# METHOD 1.2

total <- (D[1, 2] * X[1, 1] * (y[1] -
                                 X[1, 1] * mu[2, 1] * alpha[2, 1] -
                                 X[1, 2] * mu[2, 2] * alpha[2, 2] -
                                 X[1, 3] * mu[2, 3] * alpha[2, 3])) +
  (D[2, 2] * X[2, 1] * (y[2] -
                          X[2, 1] * mu[2, 1] * alpha[2, 1] -
                          X[2, 2] * mu[2, 2] * alpha[2, 2] -
                          X[2, 3] * mu[2, 3] * alpha[2, 3])) +
  (D[3, 2] * X[3, 1] * (y[3] -
                          X[3, 1] * mu[2, 1] * alpha[2, 1] -
                          X[3, 2] * mu[2, 2] * alpha[2, 2] -
                          X[3, 3] * mu[2, 3] * alpha[2, 3])) +
  (D[4, 2] * X[4, 1] * (y[4] -
                          X[4, 1] * mu[2, 1] * alpha[2, 1] -
                          X[4, 2] * mu[2, 2] * alpha[2, 2] -
                          X[4, 3] * mu[2, 3] * alpha[2, 3]))

residual <- ((D[1, 2] * X[1, 1] * (X[1, 1] * mu[2, 1] * alpha[2, 1])) +
               (D[2, 2] * X[2, 1] * (X[2, 1] * mu[2, 1] * alpha[2, 1])) +
               (D[3, 2] * X[3, 1] * (X[3, 1] * mu[2, 1] * alpha[2, 1])) +
               (D[4, 2] * X[4, 1] * (X[4, 1] * mu[2, 1] * alpha[2, 1])))
total + residual

# METHOD 2: update using for loops
mu_mat <- matrix(NA, n, p)
for (l in 1:n){

  # vector to ensure that the k-th entry of the l-th row of mu_mat and alpha
  # mat are not summed
  ind0 <- rep(1, p)

  # update the l, k entry of mu
  for (k in 1:p){

    # put a 0 in the k-th position of ind0 to 0 out the k-th entry of mu_stack
    # and alpha_stack
    ind0.k <- ind0; ind0.k[k] <- 0

    # the l-th row of mu_mat, alpha_mat with the k-th entry 0'd, stacked n times
    mu_stack <- matrix(mu[l, ] * ind0.k, n, p, T)
    alpha_stack <- matrix(alpha[l, ] * ind0.k, n, p, T)

    # update mu
    mu_mat[l, k] <- sum((D[ , l] * X[ , k]) *
                           (y - rowSums(X * mu_stack * alpha_stack)))
  }
}
mu_mat


# METHOD 3: update using old method

# initialize the necessary variables
Mu_vec <- matrix(t(mu), n * p, 1)
alpha_vec <- matrix(t(alpha), n * p, 1)
y_long_vec <- as.vector(t(y%*%matrix(1,1,p)))
X_vec <- matrix(0, n * p, 1)
X2 <- matrix(rep(0, n^2 * p), nrow = n, ncol = n * p)
D_long <- matrix(0,n*p,n)

for (i in 1:n){
  D_long[ , i] <- matrix(t(D[, i] %*% matrix(1, 1, p)),n * p, 1)
}

for (i in 1:n) {
  for (j in 1:p) {
    k <- p * (i - 1) + 1
    X2[i, k + j - 1] <- X[i, j]
    X_vec[k + j - 1] <- X2[i, k + j - 1]
  }
}

mu_mat3 <- matrix(NA, n, p)

# update loop
for (i in 1:n) {
  y_XW <- y_long_vec * X_vec * D_long[, i]
  y_XW_mat <- matrix(y_XW, n, p, byrow = TRUE)

  X_mu_alpha <- X_vec * Mu_vec * alpha_vec
  xmualpha_mat <- t(matrix(X_mu_alpha, p, n)) %*% (matrix(1, p, p) - diag(rep(1, p)))
  XW_mat <- matrix(X_vec * D_long[, i], n, p, byrow = TRUE) * xmualpha_mat

  mu_mat3[i, ] <- (t(y_XW_mat) %*% matrix(1, n, 1) - (t(XW_mat) %*% matrix(1, n, 1)))
}

mu_mat3


# METHOD 4: manual update of mu_1^2 using original coding method
# (the 2, 1 entry of the mu matrix)
(y[1] * X[1, 1] * D[1, 2] + y[2] * X[2, 1] * D[2, 2] + y[3] * X[3, 1] * D[3, 2] +
  y[4] * X[4, 1] * D[4, 2] -
  X[1, 1] * D[1, 2] * (X[1, 2] * mu[1, 2] * alpha[1, 2] +
                         X[1, 3] * mu[1, 3] * alpha[1, 3]) -
  X[2, 1] * D[2, 2] * (X[2, 2] * mu[2, 2] * alpha[2, 2] +
                         X[2, 3] * mu[2, 3] * alpha[2, 3]) -
  X[3, 1] * D[3, 2] * (X[3, 2] * mu[3, 2] * alpha[3, 2] +
                         X[3, 3] * mu[3, 3] * alpha[3, 3]) -
  X[4, 1] * D[4, 2] * (X[4, 2] * mu[4, 2] * alpha[4, 2] +
                         X[4, 3] * mu[4, 3] * alpha[4, 3]))

