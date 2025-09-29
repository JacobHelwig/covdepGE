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
p <- 5
(n <- 2 * 3 * 3 * p)
(nj <- n %/% 3)

# p <- 5
# n <- 180
# (nj <- n %/% 3)

# overload generateData
# tanhp <- function(x,p=50)(exp(p * x) - 1) / (2 * (exp(p * x) + 1))
# ggplot() + geom_function(fun = function(Z) tanhp(Z + 1) + 0.5, n = 1000) + xlim(-5, 5) + geom_vline(xintercept = 1) + geom_vline(xintercept = -1)
# ggplot() + geom_function(fun = function(Z) -tanhp(Z - 1) + 0.5, n = 1000) + xlim(-5, 5) + geom_vline(xintercept = 1) + geom_vline(xintercept = -1)
generateData <- function(p = 5, n1 = 60, n2 = 60, n3 = 60){

  # create covariate for observations in each of the three intervals

  # define number of samples
  n <- sum(n1, n2, n3)

  # define the intervals
  limits1 <- c(-3, -1)
  limits2 <- c(-1, 1)
  limits3 <- c(1, 3)

  # generate Z

  # define the interval labels
  interval <- c(rep(1, n1), rep(2, n2), rep(3, n3))

  # draw the covariate values within each interval
  z1 <- sort(stats::runif(n1, limits1[1], limits1[2]))
  z2 <- sort(stats::runif(n2, limits2[1], limits2[2]))
  z3 <- sort(stats::runif(n3, limits3[1], limits3[2]))
  Z <- matrix(c(z1, z2, z3), n, 1)

  # create the precision matrices

  # define precision matrix function
  tanhp <- function(x,p=50)(exp(p * x) - 1) / (2 * (exp(p * x) + 1))

  # the shared part of the structure for all three intervals is a 2 on the
  # diagonal and a 1 in the (2, 3) position
  common_str <- diag(p)
  common_str[2, 3] <- 1

  # all matrices have a function of Z in the (1, 2) position and (1, 3)
  # positions; define structures for each of these components
  str12 <- str13 <- matrix(0, p, p)
  str12[1, 2] <- str13[1, 3] <- 1

  # define the precision matrices for each of the observations
  prec_mats <- lapply(Z, function(z) common_str +
                        ((-tanhp(z - 1) + 0.5) * str12) +
                        ((tanhp(z + 1) + 0.5) * str13))

  # symmetrize the precision matrices
  true_precision <- lapply(prec_mats, function(mat) t(mat) + mat) # lapply(1:length(true_precision), function(j)matViz(true_precision[[j]], incl_val = T) + ggtitle(paste(j, Z[j])))


  # invert the precision matrices to get the covariance matrices
  cov_mats <- lapply(true_precision, solve)

  # generate the data using the covariance matrices
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p)))

  return(list(X = data_mat, Z = Z, true_precision = true_precision,
              interval = interval))
}

# generate the data
data_list <- replicate(n_trials, generateData(p, nj, nj, nj), F)

# get number of available workers and register parallel backend
(num_workers <- parallel::detectCores() - 5)
doParallel::registerDoParallel(num_workers)

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
  clust <- Mclust(Z, verbose = F)
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

  # re-scale Z to [0, 1]
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

# function to fit and evaluate results for covdepGE
covdepGE.eval <- function(X, Z, true, n_workers){

  start <- Sys.time()

  # get dimensions of the data and fit covdepGE
  n <- nrow(X)
  p <- ncol(X)
  out <- covdepGE(X = X,
                  Z = Z)#,
                  #parallel = T,
                  #num_workers = n_workers)

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

functions <- c("aic_JGL", "eval_est", "JGL.eval", "sp.array", "tvmgm.eval", "covdepGE.eval")
packages <- c("JGL", "mclust", "mgm", "covdepGE")

# perform trials
results <- foreach(j = 1:n_trials, .export = functions,
                   .packages = packages)%dopar%
  {
    # record the time the trial started
    trial_start <- Sys.time()

    # get the data and create storage for the models (j=1)
    data <- data_list[[j]]
    trial <- vector("list", 3)
    names(trial) <- c("mgm", "JGL", "covdepGE")

    # convert the true precision to an array and then to a graph; mask diagonal
    prec <- array(unlist(data$true_precision), c(p, p, n))
    graph <- (prec != 0) * 1 + replicate(n, diag(rep(NA, p)) * 1)

    # fit each method

    # mgm
    out_mgm <- tryCatch(tvmgm.eval(X = data$X,
                                   Z = data$Z,
                                   true = graph),
                        error = function(e) list(error = e))
    if (!is.null(out_mgm$error)){
      message("mgm ERROR:", out_mgm$error)
      next
    }
    trial$mgm <- out_mgm
    rm(list = "out_mgm")
    gc()

    # JGL
    out_JGL <- tryCatch(JGL.eval(X = data$X,
                                 Z = data$Z,
                                 true = graph),
                        error = function(e) list(error = e))
    if (!is.null(out_JGL$error)){
      message("JGL ERROR:", out_JGL$error)
      next
    }
    trial$JGL <- out_JGL
    rm(list = "out_JGL")
    gc()

    # fit covdepGE
    out_covdepGE <- tryCatch(covdepGE.eval(X = data$X,
                                           Z = data$Z,
                                           true = graph,
                                           n_workers = num_workers),
                             error = function(e) list(error = e))
    if (!is.null(out_covdepGE$error)){
      message("covdepGE ERROR:", out_covdepGE$error)
      next
    }
    trial$covdepGE <- out_covdepGE
    rm(list = "out_covdepGE")
    gc()

    # return the trial
    message("\nTrial ", j, " complete ", Sys.time(), "\n")
    trial
  }

doParallel::stopImplicitCluster()

library(kableExtra)

names(results) <- 1:n_trials

# aggregate results by model and put into a models list
covdepGE_mods <- lapply(results, `[[`, "covdepGE")
jgl_mods <- lapply(results, `[[`, "JGL")
mgm_mods <- lapply(results, `[[`, "mgm")
mods <- list(covdepGE = covdepGE_mods,
             JGL = jgl_mods,
             mgm = mgm_mods)

# extract times, sensitivity, specificity, ect.
times <- lapply(mods, sapply, `[[`, "time")
sens <- lapply(mods, sapply, `[[`, "sens")
spec <- lapply(mods, sapply, `[[`, "spec")
TP <- lapply(mods, sapply, `[[`, "TP_n")
# TN <- lapply(mods, lapply, sapply, `[[`, "TN_n")
FP <- lapply(mods, sapply, `[[`, "FP_n")
# FN <- lapply(mods, lapply, sapply, `[[`, "FN_n")

# function to get the mean and standard deviation at a specified precision
mean_sd <- function(x, prec = 4, mean_format = "f", sd_format = "f") {
  paste0("$", formatC(mean(x), prec, format = mean_format), "(",
         formatC(sd(x), prec, format = sd_format), ")$")
}

# get summary stats for each
times_sum <- lapply(times, mean_sd)
sens_sum <- lapply(sens, mean_sd)
spec_sum <- lapply(spec, mean_sd)
TP_sum <- lapply(TP, mean_sd)
FP_sum <- lapply(FP, mean_sd)

# create a matrix for each experiment
exp_sum <- list("Sensitivity$(\\uparrow)$" = unlist(sens_sum),
                "Specificity$(\\uparrow)$" = unlist(spec_sum),
                # "TP/graph$(\\uparrow)$" = TP_sum,
                # "FP/graph$(\\downarrow)$" = FP_sum,
                "Time(s)$(\\downarrow)$" = unlist(times_sum))
res_mat <- as.matrix(data.frame(exp_sum))

# combine all of the matrices
pkg_names <- paste0("\\texttt{", row.names(res_mat), "}")
colnames(res_mat) <- c(names(exp_sum))

kbl(res_mat, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle")
