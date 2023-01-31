library(covdepGE)
library(foreach)
library(loggle)
library(mgm)
library(varbvs)

# function to fit and evaluate results for covdepGE
covdepGE.eval <- function(X, Z, true, n_workers){

  start <- Sys.time()

  # get dimensions of the data and fit covdepGE
  n <- nrow(X)
  p <- ncol(X)
  out <- covdepGE(X = X,
                  Z = Z,
                  parallel = T,
                  num_workers = n_workers)

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

# function to fit and evaluate results for loggle
loggle.eval <- function(X, Z, true, n_workers){

  start <- Sys.time()

  # if the covariate is multidimensional, sort observations in X and ground truth
  if (ncol(Z) == 2){
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

  # if the covariates are discrete, only estimate the graphs in the middle of
  # each of the clusters
  if (discrete){
    pos <- sapply(Z_star, function(z) round(median(which(Z == z))))
  }else{

    # otherwise, the covariates are continuous

    # there are issues with estimating graphs at the end points of the time
    # interval; don't estimate these
    cutoff <- 10
    pos <- cutoff:(n - cutoff)

  }

  # fit loggle
  out <- R.utils::withTimeout(
    quiet(loggle.cv(t(X),
                    pos = pos,
                    d.list = c(0, 0.001, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2),
                    num.thread = n_workers)),
    timeout = 2 * 60 * 60 * n_workers)
  closeAllConnections()

  # record time and get the array of graphs
  out$time <- as.numeric(Sys.time() - start, units = "secs")
  out$str <- array(NA, dim = c(p, p, n))

  # if discrete, assign all observations in the same cluster the corresponding
  # graph
  if (discrete){
    for (j in 1:length(Z_star)){

      # get the graph estimated for the j-th unique covariate value and assign
      # it to all observations with this covariate
      graph_j <- as.matrix(out$cv.select.result$adj.mat.opt[[j]] - diag(p))
      out$str[,, Z == Z_star[j]] <- graph_j
    }
  }else{

    # otherwise, the covariate is continuous

    for (j in 1:n){

      # if the observation is in the cutoff region, assign the graph for the last
      # observation outside of the cutoff region
      if (j < cutoff){
        graph_j <- out$cv.select.result$adj.mat.opt[[1]]
      }else if (j > n - cutoff){
        graph_j <- out$cv.select.result$adj.mat.opt[[n - 2 * cutoff + 1]]
      }else{
        graph_j <- out$cv.select.result$adj.mat.opt[[j - cutoff + 1]]
      }
      out$str[,, j] <- as.matrix(graph_j - diag(p))
    }
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

# function to perform bandwidth selection, run tvmgm, and evaluate the results
tvmgm.eval <- function(X, Z, true){

  start <- Sys.time()

  # if the covariate is multidimensional, sort observations in X and ground truth
  if (ncol(Z) == 2){
    sort_inds <- sort_Z(Z)
    X <- X[sort_inds, ]
    true <- true[sort_inds]
    Z <- 1:nrow(X)
  }

  # get dimensions of the data
  n <- nrow(X)
  p <- ncol(X)

  # determine if the covariate is discete
  discrete <- length(unique(Z)) <= 2

  # re-scale Z to [0, 1]
  z01 <- Z - min(Z)
  z01 <- z01 / max(z01)

  # if the covariate are discrete, only estimate the graphs for the unique
  # timepoints
  z01_est <- z01
  if (discrete){
    z01_est <- unique(z01)
  }

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

  # if discrete, assign all observations in the same cluster the corresponding
  # graph
  if (discrete){
    str <- array(NA, dim = c(p, p, n))
    for (j in 1:length(z01_est)){
      str[,, z01 == z01_est[j]] <- out$str[,, j]
    }
    out$str <- str
  }

  # get performance, convert graphs to a sparse array, and return
  perf <- eval_est(out$str, true)
  out[names(perf)] <- perf
  out$str <- sp.array(out$str, n)
  out
}

# function to fit and evaluate results for varbvs
varbvs.eval <- function(X, Z, true){

  start <- Sys.time()

  # get dimensions of the data and fit varbvs
  n <- nrow(X)
  p <- ncol(X)

  # estimate a graph for each level of the covariate
  out <- list(str = array(NA, c(p, p, n)))
  for (z in unique(Z)){

    # fix the datapoints corresponding to z
    X_z <- X[Z == z, ]
    graph_z <- matrix(0, p, p)

    # fix each of the variables as the response and estimate the PIP
    for (k in 1:p){
      regr <- varbvs(X = X_z[ , -k],
                     Z = NULL,
                     y = X_z[ , k],
                     verbose = F)
      graph_z[k , -k] <- regr$pip
    }

    # symmetrize and threshold to get the graph
    graph_z <- (graph_z + t(graph_z)) / 2
    graph_z <- (graph_z > 0.5) * 1

    # add the graph to the graphs array
    out$str[,, Z == z] <- graph_z
  }

  # record time
  out$time <- as.numeric(Sys.time() - start, units = "secs")

  # get performance, convert graphs to a sparse array, and return
  perf <- eval_est(out$str, true)
  out[names(perf)] <- perf
  out$str <- sp.array(out$str, n)
  message("\nvarbvs complete ", Sys.time(), "\n")
  out
}
