rm(list = ls())
library(JGL)
library(mclust)
library(covdepGE)
library(foreach)

# define data dimensions
p <- 100
n <- 2 * 3 * p
(nj <- n %/% 3)

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
aic_apx <- function(X, prec){

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

# generate the data
set.seed(1)
data <- generateData(p, nj, nj, nj)

# convert the true precision to an array and then to a graph; mask diagonal
prec_arr_dim <- c(p, p, n)
prec <- array(unlist(data$true_precision), prec_arr_dim)
graph <- (prec != 0) * 1 + replicate(n, diag(rep(NA, p)) * 1)

X <- data$X
Z <- data$Z

# spin up parallel backend and cluster the data based on Z
doParallel::registerDoParallel(num_workers)
clust <- Mclust(Z, verbose = F)
X_k <- lapply(1:clust$G, function(k) X[clust$classification == k, ])

# create a grid of lambda1 and lambda2
lambda1_min <- 0.1
lambda2_min <- 1e-5
lambda1_max <- 0.3
lambda2_max <- 0.01
#lambda1 <- exp(seq(log(lambda1_min), log(lambda1_max), length = num_workers))
lambda1 <- seq(lambda1_min, lambda1_max, length = num_workers)
lambda2 <- exp(seq(log(lambda2_min), log(lambda2_max), length = num_workers))

# optimize lambda1 with lambda2 fixed as the smallest value
aic_lambda1 <- foreach(lambda = lambda1, .combine = c, .packages = "JGL")%dopar%{

  # fit the model and return the AIC
  out <- JGL(Y = X_k,
             lambda1 = lambda,
             lambda2 = lambda2_min,
             return.whole.theta = T)
  aic_apx(X_k, out$theta)
}

# fix lambda 1 and optimize lambda2
lambda1_opt <- lambda1[which.min(aic_lambda1)]
out_lambda2 <- foreach(lambda = lambda2, .packages = "JGL")%dopar%{

  # fit the model and return the AIC with the model
  out <- JGL(Y = X_k,
             lambda1 = lambda1_opt,
             lambda2 = lambda,
             return.whole.theta = T)
  list(out = out, aic = aic_apx(X_k, out$theta))
}

# save the model with the lowest AIC as the final model
opt_ind <- which.min(sapply(out_lambda2, `[[`, "aic"))
out <- out_lambda2[[opt_ind]]$out

# save the optimal lambda
out$lambda1 <- lambda1_opt
out$lambda2 <- lambda2[opt_ind]

# get the estimated graphs
n <- nrow(X)
p <- ncol(X)
out$graphs <- array(unlist(out$theta[clust$classification]), c(p, p, n))
out$graphs <- (out$graphs != 0) * 1 - replicate(n, diag(p))

# get sensitivity and specificity
eval_est(out$graphs, graph)

out$lambda1
out$lambda2
doParallel::stopImplicitCluster()
