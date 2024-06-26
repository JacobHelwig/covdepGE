library(covdepGE)
library(loggle)
library(JGL)
library(mclust)
library(mgm)

# function to fit and evaluate results for covdepGE
covdepGE.eval <- function(X, Z, hp_method, true, n_workers, max_iter_grid){

  start <- Sys.time()

  # get dimensions of the data and fit covdepGE
  n <- nrow(X)
  p <- ncol(X)
  out <- covdepGE(X = X,
                  Z = Z,
                  hp_method = hp_method,
                  parallel = n_workers > 1,
                  num_workers = n_workers,
                  max_iter_grid = max_iter_grid,
                  prog_bar=FALSE)

  # record time and get the array of graphs
  out$time <- as.numeric(Sys.time() - start, units = "secs")
  out$str <- array(unlist(out$graphs$graphs), dim = c(p, p, n))

  # covert the unique graphs to a list of sparse matrices
  out$unique_graphs <- out$graphs$unique_graphs
  for (k in 1:length(out$unique_graphs)){
    out$unique_graphs[[k]]$graph <- Matrix::Matrix(
      out$unique_graphs[[k]]$graph, sparse = T)
  }

  # remove large objects, put the unique graphs back in the graphs sublist
  out$pip <- out$graphs$inclusion_probs_sym
  out$variational_params <- out$graphs <- out$weights <- NULL
  out$graphs$unique_graphs <- out$unique_graphs
  out$unique_graphs <- NULL

  # get performance, convert graphs to a sparse array, and return
  perf <- eval_est(out$str, true)
  out[names(perf)] <- perf
  out$str <- sp.array(out$str, n)
  message("\ncovdepGE complete ", Sys.time(), "\n")
  out
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

  # verify that the aic is valid and return
  aic <- ifelse(is.numeric(aic), aic, Inf)
  aic
}

# function to perform clustering, cross-validation and evaluation for JGL
JGL.eval <- function(X, Z, true){

  start0 <- Sys.time()

  # cluster the data based on Z
  clust <- Mclust(Z, verbose = F, G = 2:12)
  X_k <- lapply(1:clust$G, function(k) X[clust$classification == k, ])

  # create a grid of lambda1 and lambda2
  lambda1_min <- 0.15
  lambda2_min <- 1e-5
  lambda1_max <- 0.4
  lambda2_max <- 0.01
  lambda1 <- seq(lambda1_min, lambda1_max, 0.005)
  lambda2 <- exp(seq(log(lambda2_min), log(lambda2_max),
                     length = length(lambda1) %/% 2))

  # optimize lambda1 with lambda2 fixed as the smallest value
  aic_lambda1 <- vector("list", length(lambda1))
  for(k in 1:length(lambda1)){

    # fit the model and return lambda, AIC, and time to fit
    start <- Sys.time()
    out <- JGL(Y = X_k,
               lambda1 = lambda1[k],
               lambda2 = lambda2_min,
               return.whole.theta = T)
    time <- as.numeric(Sys.time() - start, units = "secs")
    aic_lambda1[[k]] <- list(lambda = lambda1[k], aic = aic_JGL(X_k, out$theta),
                             time = time)
  }

  # fix lambda 1 and optimize lambda2
  lambda1_opt <- sapply(aic_lambda1, `[[`, "aic")
  lambda1_opt <- lambda1[which.min(lambda1_opt)]
  aic_lambda2 <- vector("list", length(lambda2))
  for(k in 1:length(lambda2)){

    # fit the model and return lambda, AIC, and time to fit
    start <- Sys.time()
    out <- JGL(Y = X_k,
               lambda1 = lambda1_opt,
               lambda2 = lambda2[k],
               return.whole.theta = T)
    time <- as.numeric(Sys.time() - start, units = "secs")
    aic_lambda2[[k]] <- list(lambda = lambda2[k], aic = aic_JGL(X_k, out$theta),
                             time = time)
  }

  # select the optimal lambda2 and fit the final model
  lambda2_opt <- sapply(aic_lambda2, `[[`, "aic")
  lambda2_opt <- lambda2[which.min(lambda2_opt)]
  out <- JGL(Y = X_k,
             lambda1 = lambda1_opt,
             lambda2 = lambda2_opt,
             return.whole.theta = T)

  # record time
  out$time <- as.numeric(Sys.time() - start0, units = "secs")

  # save the lambda grid, optimal lambda and classification
  out$lambda1_grid <- aic_lambda1
  out$lambda2_grid <- aic_lambda2
  out$lambda1 <- lambda1_opt
  out$lambda2 <- lambda2_opt
  out$classification <- clust$classification

  # get the estimated graphs
  n <- nrow(X)
  p <- ncol(X)
  out$str <- array(unlist(out$theta[clust$classification]), c(p, p, n))
  out$str <- (out$str != 0) * 1 - replicate(n, diag(p))

  # get performance, convert graphs to a sparse array, and return
  perf <- eval_est(out$str, true)
  out[names(perf)] <- perf
  out$str <- sp.array(out$str, n)
  out
}

# function to perform bandwidth selection, run tvmgm, and evaluate the results
tvmgm.eval <- function(X, Z, true){

  start <- Sys.time()

  # if the covariate is multidimensional, sort observations in X and ground truth
  if (ncol(Z) > 1){
    sort_inds <- sort_Z(Z)
    X <- X[sort_inds, ]
    true <- true[sort_inds]
    Z <- 1:nrow(X)
  }

  # get dimensions of the data
  n <- nrow(X)
  p <- ncol(X)

  # re-scale Z to [0, 1]
  z01 <- Z - min(Z)
  z01_est <- z01 <- z01 / max(z01)

  # choose optimal bandwidth
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
               estpoints = z01_est,
               bandwidth = bw,
               k = 2,
               pbar = F)

  # record the time
  out$time <- as.numeric(Sys.time() - start, units = "secs")

  # save the selected bandwidth and remove large objects
  out$bw <- bw
  out$tvmodels <- out$interactions <- out$intercepts <- NULL

  # get graphs, remove pairwise (it is large)
  out$str <- (out$pairwise$wadj != 0) * 1
  out$pairwise <- NULL

  # get performance, convert graphs to a sparse array, and return
  perf <- eval_est(out$str, true)
  out[names(perf)] <- perf
  out$str <- sp.array(out$str, n)
  out
}

# function to fit and evaluate results for loggle
loggle.eval <- function(X, Z, true, n_workers, cutoff){

  start <- Sys.time()

  # if the covariate is multidimensional, sort observations in X and ground truth
  if (ncol(Z) > 1){
    sort_inds <- sort_Z(Z)
    X <- X[sort_inds, ]
    true <- true[sort_inds]
  }

  # get dimensions of the data
  n <- nrow(X)
  p <- ncol(X)

  # determine if the covariate is discrete
  Z_star <- unique(Z)
  discrete <- length(Z_star) <= 2

  # there are issues with estimating graphs at the end points of the time
  # interval; don't estimate these
  pos <- (cutoff + 1):(n - cutoff)
  # pos <- cutoff:(n - cutoff)

  # fit loggle
  # out <- R.utils::withTimeout(
  #   quiet(loggle.cv(t(X),
  #                   pos = pos,
  #                   d.list = c(0, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2),
  #                   num.thread = n_workers)),
  #   timeout = 2 * 60 * 60 * n_workers)
  out <- quiet(loggle.cv(t(X),
                         pos = pos,
                         d.list = c(0, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2),
                         num.thread = n_workers))
  closeAllConnections()

  # record time and get the array of graphs
  out$time <- as.numeric(Sys.time() - start, units = "secs")
  out$str <- array(NA, dim = c(p, p, n))


  for (j in 1:n){

    # if the observation is in the cutoff region, assign the graph for the last
    # observation outside of the cutoff region
    if (j %in% pos){
      graph_j <- out$cv.select.result$adj.mat.opt[[j - cutoff]]
      # print(paste('in', j, j-cutoff))
    }else if (j < cutoff + 1){
      # print(j)
      graph_j <- out$cv.select.result$adj.mat.opt[[1]]
    }else if (j > n - cutoff){
      # print(paste(j, n - 2 * cutoff, length(pos)))
      graph_j <- out$cv.select.result$adj.mat.opt[[n - 2 * cutoff]]
    }else{
      stop(paste0('Error in resolving cutoff regions for j=', j))
    }
    out$str[,, j] <- as.matrix(graph_j - diag(p))
  }

  # remove large objects
  out$cv.result.h <- NULL

  # get performance, convert graphs to a sparse array, and return
  perf <- eval_est(out$str, true)
  out[names(perf)] <- perf
  out$str <- sp.array(out$str, n)
  message("\nloggle complete ", Sys.time(), "\n")
  out
}
