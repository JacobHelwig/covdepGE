# function for calculating ELBO
ELBO <- function(y, X, D, mu, alpha, S_sq, sigmasq, sigmabeta_sq, pi_est){

  n <- nrow(X); p <- ncol(X)

  t1 <- t2 <- 0
  for (i in 1:n){
    sum_t1 <- 0
    sum_t2 <- 0
    for (j in 1:p){
      sum_t1 <- sum_t1 + X[i, j] * mu[j] * alpha[j]
      sum_t2 <- sum_t2 + X[i, j]^2 * (alpha[j] * (mu[j]^2 + S_sq[j]) - alpha[j]^2 * mu[j]^2)
    }
    t1 <- t1 + D[i] * (y[i] - sum_t1)^2
    t2 <- t2 + D[i] * sum_t2
  }
  t1 <- -t1 / (2 * sigmasq)
  t2 <- -t2 / (2 * sigmasq)

  t3 <- t4 <- t5 <- 0
  for (j in 1:p){
    t3 <- t3 + alpha[j] * (1 + log(S_sq[j]))
    t4 <- t4 + alpha[j] * log((alpha[j] + 0.000001) / pi_est) + (1 - alpha[j]) * log((1 - alpha[j] + 0.000001) / (1 - pi_est))
    t5 <- t5 + alpha[j] * ((mu[j]^2 + S_sq[j]) / (2 * sigmasq * sigmabeta_sq) + log(sigmasq * sigmabeta_sq) / 2)
  }
  t3 <- t3 / 2
  t4 <- -t4
  t5 <- -t5

  t6 <- log(1 / (2 * pi * sigmasq)) / 2

  return(sum(t1, t2, t3, t4, t5, t6))
}

# function to calculate the ELBO for all of the individuals
ELBO_tot <- function(y, X, D, mu, alpha, S_sq, sigmasq, sigmabeta_sq, pi_est){

  n <- nrow(X)

  ELBO_tot <- 0
  for (l in 1:n){
    ELBO_tot <- ELBO_tot + ELBO(y, X, D[ , l], as.numeric(mu[l , ]), as.numeric(alpha[l , ]), as.numeric(S_sq[l , ]), sigmasq[l], sigmabeta_sq[l], pi_est)
  }
  ELBO_tot
}

# S_sq update
Ssq_update <- function(X, D, sigmasq, sigmabeta_sq){

  n <- nrow(X); p <- ncol(X)
  new_Ssq <- matrix(NA, n, p)

  # update S_sq for each individual
  for (l in 1:n){

    Ssq_l <- rep(NA, p)
    # update the k-th coordinate of S_sq
    for (k in 1:p){
      sum_lk <- 0
      for (i in 1:n){
        sum_lk <- sum_lk + X[i, k]^2 * D[i, l]
      }
      Ssq_l[k] <- sigmasq[l] / (1 / sigmabeta_sq[l] + sum_lk)
    }
    new_Ssq[l , ] <- Ssq_l
  }

  new_Ssq
}

# variational update for mu
mu_update <- function(y, X, D, mu, alpha, S_sq, sigmasq){

  n <- nrow(X); p <- ncol(X)

  new_mu <- matrix(NA, n, p)

  for (l in 1:n){

    # update mu_l
    mu_l <- rep(NA, p)
    for (k in 1:p){

      # update the k-th coordinate of mu_l
      mu_lk <- 0
      for (i in 1:n){

        sum_lk <- 0
        for (m in 1:p){

          if (m != k){
            sum_lk <- sum_lk + X[i, m] * mu[l, m] * alpha[l, m]
          }
        }
        mu_lk <- mu_lk + D[i, l] * X[i, k] * (y[i] - sum_lk)
      }
      mu_l[k] <- S_sq[l, k] / sigmasq[l] * mu_lk
    }
    new_mu[l, ] <- mu_l
  }
  new_mu
}

# variational update for alpha
alpha_update <- function(mu, S_sq, sigmasq, sigmabeta_sq, pi_est){
  n <- nrow(mu); p <- ncol(mu)

  new_alpha <- matrix(NA, n, p)

  logit_pi <- log(pi_est / (1 - pi_est))

  # update alpha for individual l
  for (l in 1:n){

    alpha_l <- rep(NA, p)

    # update the k-th coordindate of alpha for individual l
    for (k in 1:p){
      alpha_l[k] <- logit_pi + mu[l, k]^2 / (2 * S_sq[l, k]) + log(sqrt(S_sq[l, k]) / (sqrt(sigmasq[l] * sigmabeta_sq[l])))
    }
    new_alpha[l, ] <- alpha_l
  }
  rep1 <- new_alpha > 9
  new_alpha <- exp(new_alpha) / (1 + exp(new_alpha))
  new_alpha[rep1] <- 1
  new_alpha
}

variance_update <- function(y, X, D, mu, alpha, S_sq, sigmabeta_sq){

  n <- nrow(X); p <- ncol(X)
  new_sigmasq <- new_sigmabeta_sq <- rep(NA, n)

  # update the variance for each individual
  for (l in 1:n){

    # calculate X and y weighted with respect to individual l
    X_w <- matrix(NA, n, p)
    y_w <- rep(NA)
    for (i in 1:n){
      X_w[i, ] <- sqrt(D[i, l]) * X[i, ]
      y_w[i] <- sqrt(D[i, l]) * y[i]
    }

    # calculate the expected value of beta for individual l
    e_beta <- mu[l, ] * alpha[l, ]


    t1 <- sum((y_w - X_w %*% e_beta)^2)

    t2 <- t3 <- num <- den <- 0
    t4 <- n
    for (k in 1:p){
      t2 <- t2 + sum((X_w[ , k])^2) * (alpha[l, k] * S_sq[l, k] + alpha[l, k] * mu[l, k]^2 - alpha[l, k]^2 * mu[l, k]^2)
      t3 <- t3 + alpha[l, k] * (S_sq[l, k] + mu[l, k]^2)
      t4 <- t4 + alpha[l, k]
      num <- num + alpha[l, k] * (S_sq[l, k] + mu[l, k]^2)
      den <- den + alpha[l, k]
    }
    t3 <- t3 / sigmabeta_sq[l]
    new_sigmasq[l] <- (t1 + t2 + t3) / t4
    den <- new_sigmasq[l] * den
    new_sigmabeta_sq[l] <- num / den
  }

  list(sigmasq = new_sigmasq, sigmabeta_sq = new_sigmabeta_sq)
}

# CAVI function
CAVI <- function(y, X, D, mu, alpha, sigmasq, sigmabeta_sq, pi_est, max_iter, tol){

  for (i in 1:max_iter){

    alpha_last <- alpha

    # update S_sq
    S_sq <- Ssq_update(X, D, sigmasq, sigmabeta_sq)

    # update mu
    mu <- mu_update(y, X, D, mu, alpha, S_sq, sigmasq)

    # update alpha
    alpha <- alpha_update(mu, S_sq, sigmasq, sigmabeta_sq, pi_est)

    # change in alpha
    delta_alpha <- norm(alpha - alpha_last, "f")
    if (delta_alpha < tol){
      break
    }

    # update variance hyperparameters
    var_upd <- variance_update(y, X, D, mu, alpha, S_sq, sigmabeta_sq)
    sigmasq <- var_upd$sigmasq
    sigmabeta_sq <- var_upd$sigmabeta_sq
  }

  ELBO <- ELBO_tot(y, X, D, mu, alpha, S_sq, sigmasq, sigmabeta_sq, pi_est)

  return(list(alpha = alpha, mu = mu, ELBO = ELBO))
}
