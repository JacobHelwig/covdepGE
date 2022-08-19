rm(list = ls())
library(covdepGE)
library(foreach)
library(JGL)
library(mclust)
library(mgm)

(now <- format(Sys.time(), "%Y%m%d_%H%M%S"))

# initialize storage for results, time, and progress tracking
set.seed(1)
n_trials <- 100
results <- vector("list", n_trials)
names(results) <- c(paste0("trial", 1:n_trials))
pb <- txtProgressBar(0, n_trials, style = 3)

# define data dimensions
p <- 100
(n <- 2 * 3 * p)
(nj <- n %/% 3)

# p <- 5
# n <- 180
# (nj <- n %/% 3)

prec_arr_dim <- c(p, p, n)

# get number of available workers
(num_workers <- parallel::detectCores() - 1)

# function for evaluating estimated graphs compared to ground truth
eval_est <- function(est, true){

  # get true number of edges and non-edges
  num_edge <- sum(true, na.rm = T)
  num_non <- sum(true == 0, na.rm = T)

  # calculate sensitivity and specificity
  true_edge <- sum(est == true & true == 1, na.rm = T)
  true_non <- sum(est == true & true == 0, na.rm = T)
  sens <- true_edge / num_edge
  spec <- true_non / num_non

  list(sens = sens, spec = spec)
}

# function to approximate the AIC for JGL
aic_JGL <- function(X, prec){

  # iterate over each of the clusters
  aic <- 0
  for (k in 1:length(X)){

    # fix the data for k-th cluster; get n and covariance
    n_k <- nrow(X[[k]])
    cov_k <- cov(X[[k]])

    # 3 terms in AIC
    aic1 <- n_k * sum(diag(cov_k %*% prec[[k]]))
    aic2 <- -n_k * log(det(prec[[k]]))
    aic3 <- 2 * sum(prec[[k]] != 0)
    aic <- aic + aic1 + aic2 + aic3
  }
  aic
}

# function to perform clustering and cross-validation for JGL
JGL.cv <- function(X, Z, lambda1_min = 0.1, lambda2_min = 1e-5,
                   lambda1_max = 0.3, lambda2_max = 0.01, num_workers){

  start <- Sys.time()

  # spin up parallel backend and cluster the data based on Z
  doParallel::registerDoParallel(num_workers)
  clust <- Mclust(Z, verbose = F)
  X_k <- lapply(1:clust$G, function(k) X[clust$classification == k, ])

  # create a grid of lambda1 and lambda2
  lambda1 <- seq(lambda1_min, lambda1_max, length = num_workers)
  lambda2 <- exp(seq(log(lambda2_min), log(lambda2_max), length = num_workers))

  # optimize lambda1 with lambda2 fixed as the smallest value
  aic_lambda1 <- foreach(lambda = lambda1, .combine = c, .export = "aic_JGL",
                         .packages = "JGL")%dopar%
    {

      # fit the model and return the AIC
      out <- JGL(Y = X_k,
                 lambda1 = lambda,
                 lambda2 = lambda2_min,
                 return.whole.theta = T)
      aic_JGL(X_k, out$theta)
  }

  # fix lambda 1 and optimize lambda2
  lambda1_opt <- lambda1[which.min(aic_lambda1)]
  out_lambda2 <- foreach(lambda = lambda2, .export = "aic_JGL",
                         .packages = "JGL")%dopar%
    {

      # fit the model and return the AIC with the model
      out <- JGL(Y = X_k,
                 lambda1 = lambda1_opt,
                 lambda2 = lambda,
                 return.whole.theta = T)
      list(out = out, aic = aic_JGL(X_k, out$theta))
    }

  # stop cluster
  doParallel::stopImplicitCluster()

  # save the model with the lowest AIC as the final model
  opt_ind <- which.min(sapply(out_lambda2, `[[`, "aic"))
  out <- out_lambda2[[opt_ind]]$out

  # save the lambda grid, optimal lambda and classification
  out$lambda1_grid <- lambda1
  out$lambda2_grid <- lambda2
  out$lambda1 <- lambda1_opt
  out$lambda2 <- lambda2[opt_ind]
  out$classification <- clust$classification

  # record time
  out$time <- as.numeric(Sys.time() - start, units = "secs")

  # get the estimated graphs and return
  n <- nrow(X)
  p <- ncol(X)
  out$str <- array(unlist(out$theta[clust$classification]), c(p, p, n))
  out$str <- (out$str != 0) * 1 - replicate(n, diag(p))
  out
}

# function to perform bandwidth selection and inference for tvmgm in parallel
tvmgm.par <- function(X, Z, num_workers){

  start <- Sys.time()

  # spin up parallel backend and re-scale Z to [0, 1]
  doParallel::registerDoParallel(num_workers)
  z01 <- Z - min(Z)
  z01 <- z01 / max(z01)

  # choose optimal bandwidth
  p <- ncol(X)
  bw <- bwSelect(data = X,
                 type = rep("g", p),
                 level = rep(1, p),
                 bwSeq = seq(0.1, 0.4, 0.1),
                 bwFolds = 5,
                 bwFoldsize = 5,
                 modeltype = "mgm",
                 k = 2,
                 pbar = F,
                 timepoints = z01)
  bw <- as.numeric(names(which.min(bw$meanError)))

  # run tvmgm
  out <- tvmgm(data = X,
               type = rep("g", p),
               level = rep(1, p),
               timepoints = z01,
               estpoints = z01,
               bandwidth = bw,
               k = 2,
               pbar = F)

  # save the selected bandwidth
  out$bw <- bw

  # stop the cluster and record time
  doParallel::stopImplicitCluster()
  out$time <- as.numeric(Sys.time() - start, units = "secs")

  # get graphs and return
  out$str <- (out$pairwise$wadj != 0) * 1
  out
}

# perform trials
for (j in 1:n_trials){

  # record the time the trial started
  trial_start <- Sys.time()

  # generate the data and create storage for the models
  data <- generateData(p, nj, nj, nj)
  trial <- vector("list", 4)
  names(trial) <- c("data", "covdepGE", "mgm", "JGL")
  trial$data <- data

  # convert the true precision to an array and then to a graph; mask diagonal
  prec <- array(unlist(data$true_precision), prec_arr_dim)
  graph <- (prec != 0) * 1 + replicate(n, diag(rep(NA, p)) * 1)

  # fit each method, save details about results, time, etc.

  # mgm
  out_mgm <- tvmgm.par(X = data$X,
                       Z = data$Z,
                       num_workers = num_workers)
  out_mgm[c("sens", "spec")] <- eval_est(out_mgm$str, graph)[c("sens", "spec")]
  trial$mgm <- out_mgm
  rm(list = "out_mgm")
  gc()

  # JGL
  out_JGL <- JGL.cv(X = data$X,
                    Z = data$Z,
                    num_workers = num_workers)
  out_JGL[c("sens", "spec")] <- eval_est(out_JGL$str, graph)[c("sens", "spec")]
  trial$JGL <- out_JGL
  rm(list = "out_JGL")
  gc()

  # covdepGE
  out_covdepGE <- covdepGE(X = data$X,
                           Z = data$Z,
                           parallel = T,
                           num_workers = num_workers)
  out_covdepGE$time <- as.numeric(out_covdepGE$model_details$elapsed, units = "secs")
  out_covdepGE$str <- array(unlist(out_covdepGE$graphs$graphs), dim = prec_arr_dim)
  out_covdepGE[c("sens", "spec")] <- eval_est(out_covdepGE$str, graph)[c("sens", "spec")]
  trial$covdepGE <- out_covdepGE
  rm(list = "out_covdepGE")
  gc()

  # record the elapsed time for the trial
  trial$data$trial_time <- as.numeric(Sys.time() - trial_start, units = "secs")

  # save the trial and update the progress bar
  results[[j]] <- trial
  setTxtProgressBar(pb, j)
  save(results, file = paste0("results_", now, ".Rda"))
}
