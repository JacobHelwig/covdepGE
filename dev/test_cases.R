# package
library(covdepGE)
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/data.R")

# src the scripts
if ("covdepGE" %in% .packages()) detach("package:covdepGE", unload = TRUE)
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/main.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/cavi.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/weights.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/plots.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/data.R")
Rcpp::sourceCpp("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/src/covdepGE_c.cpp")

set.seed(1)

# gen the data
cont <- generate_continuous()
X <- cont$data
Z <- cont$covts

# vanilla
out1 <- covdepGE(X, Z) # hybrid
out2 <- covdepGE(X, Z, hp_method = "grid_search")
out3 <- covdepGE(X, Z, hp_method = "model_average")

# increase the number of hyperparameters
out <- covdepGE(X, Z, nssq = 7, nsbsq = 7, npip = 7)

# increase ssq mult
out <- covdepGE(X, Z, ssq_mult = 10)

# increase ssq lower
out <- covdepGE(X, Z, ssq_lower = 0.1)

# increase snr_upper
out <- covdepGE(X, Z, snr_upper = 50)

# increase sbsq lower
out <- covdepGE(X, Z, sbsq_lower = .1)

# increase pip lower
out <- covdepGE(X, Z, pip_lower = 0.2)

# provide tau
out <- covdepGE(X, Z, tau = 100000)
out <- covdepGE(X, Z, tau = 1e-16)
out <- covdepGE(X, Z, tau = seq(1e-16, 1, length.out = 180))

# change the norm
out <- covdepGE(X, Z, norm = 1)

# do not center the data
out <- covdepGE(X, Z, center_data = F)

# increase the ELBO tol
out <- covdepGE(X, Z, elbo_tol = 1000)

# increase alpha tol
out <- covdepGE(X, Z, alpha_tol = 1000)

# decrease max iter
out <- covdepGE(X, Z, max_iter = 1)

# decrease edge threshold
out <- covdepGE(X, Z, edge_threshold = 0.05)

# change symmetrization method
out <- covdepGE(X, Z, sym_method = "max")

# turn off the progress bar
out <- covdepGE(X, Z, prog_bar = F)

# run in parallel
library(covdepGE)
out <- covdepGE(X, Z, parallel = T)
out <- covdepGE(X, Z, parallel = T, num_workers = 5)

doParallel::registerDoParallel(5)
out <- covdepGE(X, Z, parallel = T)
doParallel::stopImplicitCluster()
