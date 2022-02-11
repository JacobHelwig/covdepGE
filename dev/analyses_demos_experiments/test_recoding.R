source("dev/analyses_demos_experiments/recoding_covdepGE.R")
source("R/covdepGE_R.R")

set.seed(1)

n_ <- 10; p_ <- 4
np <- n_ * p_

alpha_ <- matrix(0.2, n_, p_)
S_sq_ <- matrix(0.1, n_, p_)
mu_ <- matrix(0, n_, p_)
X_ <- matrix(rnorm(np), n_, p_)
D_ <- matrix(runif(n_ * n_), n_, n_)
y_ <- as.numeric(X_ %*% c(runif(2, -1, 1), rep(0, p_ - 2)))
sigmasq_ <- rep(var(y_), n_)
sigmabeta_sq_ <- rep(1, n_)
pi_est_ <- 0.2

ELBO(y_, X_, D_[ , 1], mu_[1 , ], alpha_[1 , ], S_sq_[1 , ], sigmasq_[1], sigmabeta_sq_[1], pi_est_)
ELBO_calculator_R(y_, D_[ , 1], X_, S_sq_[1 , ], mu_[1 , ], alpha_[1 , ], sigmasq_[1], sigmabeta_sq_[1], pi_est_)

ELBO_tot(y_, X_, D_, mu_, alpha_, S_sq_, sigmasq_, sigmabeta_sq_, pi_est_)
total_ELBO_R(y_, D_, X_, S_sq_, mu_, alpha_, sigmasq_, sigmabeta_sq_, pi_est_)

Ssq_update(X_, D_, sigmasq_, sigmabeta_sq_)

sum(mu_update(y_, X_, D_, mu_, alpha_, S_sq_, sigmasq_) - mu_update_R(y_, D_, X_, S_sq_, mu_, alpha_, sigmasq_))

sum(alpha_update(mu_, S_sq_, sigmasq_, sigmabeta_sq_, pi_est_) - alpha_update_R(S_sq_, mu_, alpha_, sigmasq_, sigmabeta_sq_, pi_est_))

lapply(1:2, function(j) variance_update(y_, X_, D_, mu_, alpha_, S_sq_, sigmabeta_sq_)[[j]]
- sigma_update_R(y_, D_, X_, S_sq_, mu_, alpha_, sigmasq_, sigmabeta_sq_)[[j]])

out1 <- CAVI(y_, X_, D_, mu_, alpha_, sigmasq_, sigmabeta_sq_, pi_est_, 1000, 1e-100)

out2 <- cavi_R(y_, D_, X_, mu_, alpha_, sigmasq_, T, sigmabeta_sq_, T, pi_est_, 1e-100, 1000)

sum(out1$alpha - out2$var_alpha)
sum(out1$mu - out2$var_mu)
