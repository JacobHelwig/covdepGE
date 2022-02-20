setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_demos_experiments/hyperparameter_specification/condition_number_analysis")
start <- Sys.time()

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

library(covdepGE)

set.seed(1)
n <- 180; p <- 24
n_trials <- 100

results <- vector("list", n_trials)

names(results) <- paste0("trial", 1:n_trials)

doParallel::registerDoParallel(14)

pb <- txtProgressBar(min = 0, max = n_trials, style = 3)

for (j in 1:n_trials){

  # generate the data
  cont <- generate_continuous(p = p)
  X <- cont$data
  Z <- cont$covts

  # run the algorithm
  out <- tryCatch(covdepGE(X, Z, max_iter = 1e3, parallel = T, warnings = F,
                           stop_cluster = F),
                  error = function(msg) as.character(msg))

  # save the data and the results
  results[[j]] <- list(data = cont, results = out)

  setTxtProgressBar(pb, j)
}

close(pb)

doParallel::stopImplicitCluster()

Sys.time() - start

#save(results, file = "cond_number_models.Rda")

# get dimensions of the data
dt_dim <- dim(results$trial1$data$data)
n <- dt_dim[1]; p <- dt_dim[2]
n_trials <- length(results)

# get all of the summaries
summs <- lapply(results, `[[`, "results")

# find the sigmasq for each model
ssq <- lapply(lapply(summs, `[[`, "hyperparameters"), `[[`, "sigmasq")

# unfold the sigmasq matrices into a n_trials x p nested list of n-vectors
ssq_nest <- lapply(1:n_trials, function(trial_ind) lapply(1:p, function(var_ind) ssq[[trial_ind]][ , var_ind]))
names(ssq_nest) <- paste0("trial", 1:n_trials)
for (j in 1:n_trials) names(ssq_nest[[j]]) <- paste0("variable", 1:p)

# find ssq that are blown: either NA or in excess of 5
blown_sig <- lapply(ssq_nest, lapply, function(sigma) (sigma > 5 | is.na(sigma)))

# find counts of individuals by variable and trial
blown_ct_var <- lapply(blown_sig_nest, sapply, sum)

# find counts by trials of variables with blown sigmas
blown_ct_tr <- sapply(blown_ct_var, function(ssq) sum(ssq > 0))

# logical for trials with at least one blown variable
blown_trials <- blown_ct_tr > 0

# find blown individual indices for each variable
blown_inds <- lapply(blown_sig[blown_trials], sapply, which)

# filter out the zero length variables from blown_inds
lapply(blown_inds, function(trial) trial[sapply(trial, function(variable) length(variable) > 0)])

# all of the data
datas <- lapply(results, `[[`, "data")

# all of the data mats
data_mats <- lapply(datas, `[[`, "data")

# all data matrices with/ without at least one blown ind
blown_mats <- data_mats[blown_trials]
unblown_mats <- data_mats[!blown_trials]

# find the condition number for both
blown_cond_nums <- sapply(blown_mats, kappa)
unblown_cond_nums <- sapply(unblown_mats, kappa)
summary(blown_cond_nums)
summary(unblown_cond_nums)

# find the weights
wts <- lapply(summs, `[[`, "weights")

# find the weighted data_mats for everyone
wted_mats <- lapply(1:n_trials, function(trial_ind) lapply(1:n, function(ind) wts[[trial_ind]][ , ind] * data_mats[[trial_ind]]))
wted_mats <- lapply(1:n_trials, function(tr_ind)
  lapply(1:p, function(var_ind)
    lapply(1:n, function(indv_ind) sqrt(wts[[tr_ind]][ , indv_ind]) * data_mats[[tr_ind]][ , -var_ind])))

# find the condition
