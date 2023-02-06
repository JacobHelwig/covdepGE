# ------------------------------------------------------------------------------
#
# Testing script for performing multiple trials with each of the hyperparameter
# specification strategies in the q=1 setting
#
# Arguments:
#
# n: sample size; should be divisible by 3. 225 by default
# p: number of variables. 10 by default
# n_trials: number of trials. 5 by default
#
# Examples:
#
# For command line usage:
#
# R.exe CMD BATCH --no-save --no-restore "--args n=225 p=10 n_trials=5" test.R test_n225_p10_ntrials5.Rout
#
# For running in an IDE, modify these arguments under the `IDE args` section
# below
#
# ------------------------------------------------------------------------------

library(covdepGE)
rm(list = ls())

cat("OS:", Sys.info()['sysname'], Sys.info()['release'], "\n")

# function for evaluating estimated graphs given ground truth
eval_est <- function(est, true){

  # get n and p
  n <- length(true)
  p <- nrow(true[[1]])

  # convert estimated graphs to an array
  est <- array(unlist(est), dim = c(p, p, n))

  # convert the true precision to an array and then to a graph; mask diagonal
  prec <- array(unlist(true), c(p, p, n))
  true <- (prec != 0) * 1 + replicate(n, diag(rep(NA, p)) * 1)

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

  list(sens = sens, spec = spec)
}

# parse command line arguments
args <- (commandArgs(TRUE))

# IDE args
if (interactive()){
  args <- c("n=225", "p=5", "n_trials=2")
}

if (length(args) > 0){
  cat("Args: ", args)
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
}

# check passed arguments
if (!("n" %in% ls())){
  n <- 225
}
if (!("p" %in% ls())){
  p <- 10
}
if (!("n_trials" %in% ls())){
  n_trials <- 5
}

if (n %% 3 != 0) stop("n must be divisible by 3")

message(paste0("Performing ", n_trials, " trials, n=", n, ", p=", p ))

n <- n / 3

# generate the data
set.seed(1)
data_list <- replicate(n_trials,
                       generateData(p = p, n1 = n, n2 = n, n3 = n),
                       simplify = F)

# create storage for results
metrics <- list(sens = NA, spec = NA, time = NA)
trial_list <- list(grid_search = metrics,
                   hybrid = metrics,
                   model_average = metrics)
results <- replicate(n_trials, trial_list, simplify = F)

# register parallel backend
num_workers <- min(p, parallel::detectCores() %/% 2)
message(paste("Registering parallel backend with", num_workers, "workers"))
doParallel::registerDoParallel(num_workers)

# run each trial
for (trial in 1:n_trials){

  # get data for this trial
  data <- data_list[[trial]]
  X <- data$X
  Z <- data$Z
  Omega <- data$true_precision

  # run each HP specification strategy
  for (strategy in names(results[[trial]])){

    # fit the model
    out <- covdepGE(X = X,
                    Z = Z,
                    hp_method = strategy,
                    parallel = TRUE,
                    num_workers = num_workers)

    # evaluate and save results
    perf <- eval_est(est = out$graphs$graphs, true = Omega)
    perf$time <- as.numeric(out$model_details$elapsed, units = "secs")
    results[[trial]][[strategy]][names(metrics)] <- perf[names(metrics)]

    message(paste(strategy, "complete", Sys.time()))

  }

  message(paste0("\nTrial ", trial, "/", n_trials, " complete ", Sys.time()))

}

# shut down the cluster
doParallel::stopImplicitCluster()

# aggregate and display results
prec <- 2

# process sensitivity results
sens <- sapply(results, sapply, `[[`, "sens")
sens_mean <- rowMeans(sens)
sens_mean <- sprintf(paste0("%.", prec, "f"), sens_mean * 100)
sens_sd  <- sprintf(paste0("%.", prec, "f"), apply(sens * 100, 1, sd))
sens_str <- paste0(sens_mean, " (", sens_sd, ")")

# process specificity results
spec <- sapply(results, sapply, `[[`, "spec")
spec_mean <- rowMeans(spec)
spec_mean <- sprintf(paste0("%.", prec, "f"), spec_mean * 100)
spec_sd  <- sprintf(paste0("%.", prec, "f"), apply(spec * 100, 1, sd))
spec_str <- paste0(spec_mean, " (", spec_sd, ")")

# process time results
time <- sapply(results, sapply, `[[`, "time")
time_mean <- rowMeans(time)
time_mean <- sprintf(paste0("%.", prec, "f"), time_mean)
time_sd  <- sprintf(paste0("%.", prec, "f"), apply(time, 1, sd))
time_str <- paste0(time_mean, " (", time_sd, ")")

# combine summary strings
summaries <- cbind.data.frame(sens_str, spec_str, time_str)
df <- data.frame(Strategy = row.names(spec))
df[ , c("Sensitivity(%)", "Specificity(%)", "Time(seconds)")] <- summaries

# display final result
df
