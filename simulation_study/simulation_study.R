rm(list = ls())
library(covdepGE)
library(loggle)
library(HeteroGGM)
library(mclust)

(now <- format(Sys.time(), "%Y%m%d_%H%M%S"))

# initialize storage for results, time, and progress tracking
set.seed(1)
n_trials <- 100
results <- vector("list", n_trials + 1)
names(results) <- c(paste0("trial", 1:n_trials), "time")
results$time <- Sys.time()
pb <- txtProgressBar(0, n_trials, style = 3)

# define data dimensions
p <- 25
(n <- round(2 * 3 * p))
(nj <- n %/% 3)

# get number of available workers
(num_workers <- parallel::detectCores() - 5)

# function for surpressing cat
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

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

# perform trials
for (j in 1:n_trials){

  # record the time the trial started
  trial_start <- Sys.time()

  # generate the data and create storage for the models
  data <- generateData(p, nj, nj, nj)
  trial <- vector("list", 4)
  names(trial) <- c("data", "covdepGE", "loggle", "HeteroGGM")
  prec_arr_dim <- c(p, p, n)
  trial$data <- data

  # convert the true precision to an array and then to a graph; mask diagonal
  prec <- array(unlist(data$true_precision), prec_arr_dim)
  graph <- (prec != 0) * 1 + replicate(n, diag(rep(NA, p)) * 1)

  # fit each method, save details about results, time, etc.

  # covdepGE
  doParallel::registerDoParallel(num_workers)
  out_covdepGE <- tryCatch(
    covdepGE(data$X, data$Z, parallel = T, num_workers = num_workers),
    error = function(e) list(error = e))
  doParallel::stopImplicitCluster()
  if (names(out_covdepGE)[[1]] == "error"){
    out_covdepGE <- out_covdepGE$error
  }else{
    out_covdepGE$time <- as.numeric(out_covdepGE$model_details$elapsed, units = "secs")
    out_covdepGE$str <- array(unlist(out_covdepGE$graphs$graphs), dim = prec_arr_dim)
    out_covdepGE[c("sens", "spec")] <- eval_est(out_covdepGE$str, graph)[c("sens", "spec")]
  }
  trial$covdepGE <- out_covdepGE
  rm(list = "out_covdepGE")
  gc()

  # loggle
  doParallel::registerDoParallel(num_workers)
  start <- Sys.time()
  out_loggle <- tryCatch(quiet(loggle.cv(t(data$X), num.thread = num_workers)),
                         error = function(e) list(error = e))
  if (names(out_loggle)[[1]] == "error"){
    out_loggle <- out_loggle$error
  }else{
    out_loggle <- out_loggle$cv.select.result
    out_loggle$time <- as.numeric(Sys.time() - start, units = "secs")
    out_loggle$str <- array(unlist(lapply(lapply(
      out_loggle$adj.mat.opt, `-`, diag(p)), as.matrix)), dim = prec_arr_dim)
    out_loggle[c("sens", "spec")] <- eval_est(out_loggle$str, graph)[c("sens", "spec")]
  }
  doParallel::stopImplicitCluster()
  trial$loggle <- out_loggle
  rm(list = "out_loggle")
  gc()

  # HeteroGGM
  clust <- Mclust(data$Z, verbose = F)
  lambda <- genelambda.obo(lambda2_min = 0.15)
  start <- Sys.time()
  out_hetGGM <- tryCatch(GGMPF(lambda, data$X + clust$classification * 10, clust$G),
                         error = function(e) list(error = e))
  if (names(out_hetGGM)[[1]] == "error"){
    out_hetGGM <- out_hetGGM$error
  }else{
    out_hetGGM$time <- as.numeric(Sys.time() - start, units = "secs")
    out_hetGGM$str <- (out_hetGGM$Theta_hat.list[[out_hetGGM$Opt_num]] != 0) * 1
    out_hetGGM$str <- out_hetGGM$str - replicate(dim(out_hetGGM$str)[3], diag(p))
    out_hetGGM$str <- out_hetGGM$str[ , , out_hetGGM$member.list[[out_hetGGM$Opt_num]]]
    out_hetGGM[c("sens", "spec")] <- eval_est(out_hetGGM$str, graph)[c("sens", "spec")]
  }
  trial$HeteroGGM <- out_hetGGM
  rm(list = "out_hetGGM")
  gc()

  # record the elapsed time for the trial
  trial$data$trial_time <- as.numeric(Sys.time() - trial_start, units = "secs")

  # save the trial and update the progress bar
  results[[j]] <- trial
  setTxtProgressBar(pb, j)
  save(results, file = paste0("results_", now, ".Rda"))
}
