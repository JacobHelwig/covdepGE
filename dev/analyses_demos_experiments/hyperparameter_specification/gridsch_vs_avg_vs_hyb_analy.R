setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_demos_experiments/hyperparameter_specification")

library(covdepGE)

load("gridsch_vs_avg_vs_hyb.Rda")

# separate models from the data
mods <- lapply(res, `[[`, 2)

# separate the small grid from the large grid
small <- lapply(mods, lapply, `[[`, "small")
large <- lapply(mods, lapply, `[[`, "large")

# separate the methods
small_grid <- unlist(lapply(small, `[`, "grid_search"), F)
small_avg <- unlist(lapply(small, `[`, "model_average"), F)
small_hyb <- unlist(lapply(small, `[`, "hybrid"), F)
large_grid <- unlist(lapply(large, `[`, "grid_search"), F)
large_avg <- unlist(lapply(large, `[`, "model_average"), F)
large_hyb <- unlist(lapply(large, `[`, "hybrid"), F)
names(small_grid) <- names(small_avg) <- names(small_hyb) <- names(
  large_grid) <- names(large_avg) <- names(large_hyb) <- paste0("trial", 1:length(small_grid))

# extract the model averaging weights
small_avg_weights <- lapply(lapply(small_avg, `[[`, "hyperparameters"), lapply, `[[`, "weights")
for (trial in 1:length(small_avg)){
  small_avg[[trial]]$hyperparameters <- NA
}
large_avg_weights <- lapply(lapply(large_avg, `[[`, "hyperparameters"), lapply, `[[`, "weights")
for (trial in 1:length(large_avg)){
  large_avg[[trial]]$hyperparameters <- NA
}

# save each of the separate models
save(small_grid, file = "small_grid5.Rda")
save(small_avg, file = "small_avg5.Rda")
save(small_hyb, file = "small_hyb5.Rda")
save(small_avg_weights, file = "small_avg_weights5.Rda")
save(large_grid, file = "large_grid5.Rda")
save(large_avg, file = "large_avg5.Rda")
save(large_hyb, file = "large_hyb5.Rda")
save(large_avg_weights, file = "large_avg_weights5.Rda")

# load them all
rm(list = ls())
gc()
load("small_grid5.Rda")
load("small_avg5.Rda")
load("small_hyb5.Rda")
load("large_grid5.Rda")
load("large_avg5.Rda")
load("large_hyb5.Rda")

# get the true graphs

# function for generating the data and the covariates
generate_continuous <- function(n1 = 60, n2 = 60, n3 = 60, p = 4){

  # create covariate for individuals in each of the three intervals

  # define the dimensions of the data
  n <- sum(n1, n2, n3)

  # define the limits of the intervals
  limits1 <- c(-3, -1)
  limits2 <- c(-1, 1)
  limits3 <- c(1, 3)

  # define the covariate values within each interval
  z1 <- runif(n1, limits1[1], limits1[2])
  z2 <- runif(n2, limits2[1], limits2[2])
  z3 <- runif(n3, limits3[1], limits3[2])
  Z <- matrix(sort(c(z1, z2, z3)), n, 1)

  # create precision matrices

  # the shared part of the structure for all three intervals is a 2 on the
  # diagonal and a 1 in the (2, 3) position
  common_str <- diag(p + 1)
  common_str[2, 3] <- 1

  # define constants for the structure of interval 2
  beta1 <- diff(limits2)^-1
  beta0 <- -limits2[1] * beta1

  # interval 2 has two different linear functions of Z in the (1, 2) position
  # and (1, 3) positions; define structures for each of these components
  int2_str12 <- int2_str13 <- matrix(0, p + 1, p + 1)
  int2_str12[1, 2] <- int2_str13[1, 3] <- 1

  # define the precision matrices for each of the individuals in interval 2
  int2_prec <- lapply(z2, function(z) common_str +
                        ((1 - beta0 - beta1*z)*int2_str12) +
                        ((beta0 + beta1*z)*int2_str13))

  # interval 1 has a 1 in the (1, 2) and interval 3 has a 1 in the (1, 3) position;
  # define structures for each of these components
  int1_str12 <- int3_str13 <- matrix(0, p + 1, p + 1)
  int1_str12[1, 2] <- int3_str13[1, 3] <- 1

  # define the precision matrices for each of the individuals in interval 1 and interval 3
  int1_prec <- rep(list(common_str + int1_str12), n1)
  int3_prec <- rep(list(common_str + int3_str13), n3)

  # put all of the precision matrices into one list
  prec_mats <- c(int1_prec, int2_prec, int3_prec)

  # symmetrize the precision matrices
  prec_mats <- lapply(prec_mats, function(mat) t(mat) + mat)

  # invert the precision matrices to get the covariance matrices
  cov_mats <- lapply(prec_mats, solve)

  # generate the data using the covariance matrices
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))

  return(list(data = data_mat, covts = Z, true_precision = prec_mats))
}

# function to get the number of true 1 or 0
true01 <- function(pred_graphs, true_graphs, true = 0){
  sum(sapply(1:length(pred_graphs), function(k) sum(
    pred_graphs[[k]] == true_graphs[[k]] & true_graphs[[k]] == true)))
}

# function that takes a list of models and returns the number of true 1 or 0 for each
true01_mods <- function(mods, true_graphs, true = 0){
  sapply(1:length(mods), function(j) true01(mods[[j]]$graphs, true_graphs, true))
}

cont <- generate_continuous()$true_precision
true_prec <- lapply(cont, function(prec) ((prec - diag(diag(prec))) != 0) * 1)

# true 1/0
true0 <- sum(sapply(true_prec, `==`, 0))
true1 <- sum(sapply(true_prec, `==`, 1))

# find sensitivity/ specificity for each of the methods

# grid small - specificity
spec_small_grid <- true01_mods(small_grid, true_prec) / true0
summary(spec_small_grid)

# average small - specificity
spec_small_avg <- true01_mods(small_avg, true_prec) / true0
summary(spec_small_avg)

# hybrid small - specificity
spec_small_hyb <- true01_mods(small_hyb, true_prec) / true0
summary(spec_small_hyb)

# grid small - sensitivity
sens_small_grid <- true01_mods(small_grid, true_prec, 1) / true1
summary(sens_small_grid)

# average small - sensitivity
sens_small_avg <- true01_mods(small_avg, true_prec, 1) / true1
summary(sens_small_avg)

# hybrid small - sensitivity
sens_small_hyb <- true01_mods(small_hyb, true_prec, 1) / true1
summary(sens_small_hyb)

# grid large - specificity
spec_large_grid <- true01_mods(large_grid, true_prec) / true0
summary(spec_large_grid)

# average large - specificity
spec_large_avg <- true01_mods(large_avg, true_prec) / true0
summary(spec_large_avg)

# hybrid large - specificity
spec_large_hyb <- true01_mods(large_hyb, true_prec) / true0
summary(spec_large_hyb)

# grid large - sensitivity
sens_large_grid <- true01_mods(large_grid, true_prec, 1) / true1
summary(sens_large_grid)

# average large - sensitivity
sens_large_avg <- true01_mods(large_avg, true_prec, 1) / true1
summary(sens_large_avg)

# hybrid large - sensitivity
sens_large_hyb <- true01_mods(large_hyb, true_prec, 1) / true1
summary(sens_large_hyb)

# find relative difference in sensivity and specificity with grid as the baseline

# average small - specificity
rel_spec_small_avg <- (spec_small_avg - spec_small_grid) / spec_small_grid
summary(rel_spec_small_avg)

# hybrid small - specificity
rel_spec_small_hyb <- (spec_small_hyb - spec_small_grid) / spec_small_grid
summary(rel_spec_small_hyb)

# average small - sensitivity
rel_sens_small_avg <- (sens_small_grid - sens_small_avg) / sens_small_grid
summary(rel_sens_small_avg)

# hybrid small - sensitivity
rel_sens_small_hyb <- (sens_small_grid - sens_small_hyb) / sens_small_grid
summary(rel_sens_small_hyb)

# average large - specificity
rel_spec_large_avg <- (spec_large_avg - spec_large_grid) / spec_large_grid
summary(rel_spec_large_avg)

# hybrid large - specificity
rel_spec_large_hyb <- (spec_large_hyb - spec_large_grid) / spec_large_grid
summary(rel_spec_large_hyb)

# average large - sensitivity
rel_sens_large_avg <- (sens_large_grid - sens_large_avg) / sens_large_grid
summary(rel_sens_large_avg)

# hybrid large - sensitivity
rel_sens_large_hyb <- (sens_large_grid - sens_large_hyb) / sens_large_grid
summary(rel_sens_large_hyb)

# find the importance weighted average pi

# do it for small

# get the saved importance weights
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_demos_experiments/hyperparameter_specification/small_avg_weights5.Rda")
pip_small <- small_grid$trial1$hyperparameters$`Variable 1`$grid$pip
dim(small_avg_weights$trial1$`Variable 1`)
small_avg_weights$trial1$`Variable 1` %*% pip_small

# function that takes a

small_hyb$trial1$hyperparameters$`Variable 1`$weights %*% small_hyb$trial1$hyperparameters$`Variable 1`$hyperparameters$pip
