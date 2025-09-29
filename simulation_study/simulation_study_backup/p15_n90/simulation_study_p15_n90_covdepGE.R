rm(list = ls())
library(covdepGE)

(now <- format(Sys.time(), "%Y%m%d_%H%M%S"))

# initialize storage for results, time, and progress tracking
set.seed(1)
n_trials <- 100
results <- vector("list", n_trials)
names(results) <- c(paste0("trial", 1:n_trials))
pb <- txtProgressBar(0, n_trials, style = 3)

# define data dimensions
p <- 15
(n <- 2 * 3 * p)
(nj <- n %/% 3)

# p <- 5
# n <- 180
# (nj <- n %/% 3)

# generate the data
data_list <- replicate(n_trials, generateData(p, nj, nj, nj), F)

# get number of available workers
(num_workers <- parallel::detectCores() - 5)

eval_est <- function(est, true){

  # get n
  n <- dim(est)[3]

  # get true number of edges and non-edges
  num_edge <- sum(true, na.rm = T)
  num_non <- sum(true == 0, na.rm = T)

  # calculate sensitivity, specificity, etc.
  true_edge <- sum(est == 1 & true == 1, na.rm = T)
  false_edge <- sum(est == 1 & true == 0, na.rm = T)
  true_non <- sum(est == 0 & true == 0, na.rm = T)
  false_non <- sum(est == 0 & true == 1, na.rm = T)
  sens <- true_edge / num_edge
  spec <- true_non / num_non

  list(sens = sens, spec = spec, TP_n = true_edge / n, FP_n = false_edge / n,
       TN_n = true_non / n, FN_n = false_non / n)
}

# function to turn an array into a list of sparse matrices
sp.array <- function(arr, n){
  lapply(1:n, function(l) Matrix::Matrix(arr[ , , l], sparse = T))
}

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

  # covert the unique graphs to a sparse array
  out$unique_graphs <- out$graphs$unique_graphs
  for (k in 1:length(out$unique_graphs)){
    out$unique_graphs[[k]]$graph <- Matrix::Matrix(
      out$unique_graphs[[k]]$graph, sparse = T)
  }

  # remove large objects, put the unique graphs back in the graphs sublist
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

# perform trials
for (j in 1:n_trials){

  # record the time the trial started
  trial_start <- Sys.time()

  # get the data
  data <- data_list[[j]]

  # convert the true precision to an array and then to a graph; mask diagonal
  prec <- array(unlist(data$true_precision), c(p, p, n))
  graph <- (prec != 0) * 1 + replicate(n, diag(rep(NA, p)) * 1)

  # fit covdepGE
  out_covdepGE <- tryCatch(covdepGE.eval(X = data$X,
                                         Z = data$Z,
                                         true = graph,
                                         n_workers = num_workers),
                           error = function(e) list(error = e))
  if (!is.null(out_covdepGE$error)) message(out_covdepGE$error)

  # save the trial and update the progress bar
  results[[j]] <- out_covdepGE
  setTxtProgressBar(pb, j)
  save(results, file = paste0("res_p", p, "_n", n, "_covdepGE_", now, ".Rda"))
}

save(results, file = paste0("res_p", p, "_n", n, "_covdepGE_", now, ".Rda"))
