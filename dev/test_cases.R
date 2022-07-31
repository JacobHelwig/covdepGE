library(covdepGE)
set.seed(1)

# gen the data
cont <- generateData()
X <- cont$X
Z <- cont$Z

# pass the covts
Z[61, ] <- 0
cont <- generateData(Z = Z)
cont$true_precision[[61]]

# pass the prec_mats
cont$true_precision[[61]] <- diag(5) * 0.00001
cont <- generateData(true_precision = cont$true_precision)
cont$true_precision[[61]]
cont$X[61, ]

# weird dimensions
cont <- generateData(p = 3, 4, 5, 6)

# regen
cont <- generateData()
X <- cont$X
Z <- cont$Z

# vanilla
out1 <- covdepGE(X, Z, parallel = T, num_workers = 5) # hybrid
out1$hyperparameters$variable1$final
out2 <- covdepGE(X, Z, hp_method = "grid_search", parallel = T, num_workers = 5)
out2$hyperparameters$variable1$final
out3 <- covdepGE(X, Z, hp_method = "model_average", parallel = T, num_workers = 5)
out3$hyperparameters$variable1$grid[1:5, ]

# increase the number of hyperparameters
out <- covdepGE(X, Z, nssq = 7, nsbsq = 7, npip = 7, parallel = T, num_workers = 5)
out$hyperparameters$variable5$final

# increase ssq mult
out <- covdepGE(X, Z, ssq_mult = 10, parallel = T, num_workers = 5)
unique(out$hyperparameters$variable1$grid$ssq)

# increase ssq lower
out <- covdepGE(X, Z, ssq_lower = 0.1, parallel = T, num_workers = 5)
unique(out$hyperparameters$variable1$grid$ssq)

# increase snr_upper
out <- covdepGE(X, Z, snr_upper = 50, parallel = T, num_workers = 5)
unique(out$hyperparameters$variable1$grid$sbsq)

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
