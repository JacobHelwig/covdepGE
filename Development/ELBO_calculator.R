## This function returns the ELBO for a specific variational parameter (corresponding to a fixed individual in study)
## y: response
## X_mat: data matrix except the response variable
## S_sq: variance parameter for the individual
## mu: mean parameter for that individual
## sigmasq and sigmabeta_sq: hyperparameter settings
## true_pi: hyperparameter of spike and slab mixture proportion
## W: the weights associated with the n individuals wrt the covariate of the fixed individual we are studying.
## n: no. of samples
## p: no. of variables in predictor (no. of variables - 1).

## sigmasq, sigmabeta_sq, true_pi (HYPERPAREMTERS)
## S_sq, mu, alpha (variational parameters)
## W is weights 
ELBO_calculator <- function(y, X_mat, S_sq, mu, alpha, sigmasq, sigmabeta_sq, true_pi, W, n, p) {
  mu <- matrix(mu, p, 1)
  alpha <- matrix(alpha, p, 1)
  s <- matrix(S_sq, p, 1)
  mu_alpha <- matrix(mu * alpha, p, 1)
  W <- matrix(W, n, 1)
  t1 <- -sum(W * (y - X_mat %*% mu_alpha)^2) / (2 * sigmasq)
  t2 <- -sum(W * ((X_mat)^2 %*% (alpha * (mu^2 + s) - alpha^2 * mu^2))) / (2 * sigmasq)
  t3 <- sum(alpha * ((1 + log(s)))) / 2
  # 0.0000001 to avoid division by zero since alpha = 1 often 
  t4 <- -sum(alpha * log((alpha + 0.000001) / true_pi) + (1 - alpha) * log((1 - alpha + 0.000001) / (1 - true_pi)))
  t5 <- -sum(alpha * ((mu^2 + s) / (2 * sigmasq * sigmabeta_sq) + log(sigmasq * sigmabeta_sq) / 2))
  t6 <- sum(0.5 * log(1 / (2 * pi * sigmasq)))
  t <- t1 + t2 + t3 + t4 + t5 + t6
  return(t)
}
