# function to turn an array into a list of sparse matrices
sp.array <- function(arr, n){
  lapply(1:n, function(l) Matrix::Matrix(arr[ , , l], sparse = T))
}

# function for surpressing cat
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

eval_est <- function(est, true){

  # get n and p
  n <- length(true)
  p <- nrow(true[[1]])

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

  list(sens = sens, spec = spec, TP_n = true_edge / n, FP_n = false_edge / n,
       TN_n = true_non / n, FN_n = false_non / n)
}

# function to perform trials
trials <- function(data_list, results, filename){

  # save sample data to results
  results$sample_data <- data_list[[1]]

  # get number of available workers and trials
  num_workers <- min(10, parallel::detectCores() - 5)
  n_trials <- length(data_list)

  # trials for loggle
  for (j in 1:n_trials){

    # record the time the trial started
    trial_start <- Sys.time()

    # get the data
    data <- data_list[[j]]

    # loggle
    out_loggle <- tryCatch(loggle.eval(X = data$X,
                                       Z = data$Z,
                                       true = data$true_precision,
                                       n_workers = num_workers),
                           error = function(e) list(error = e))
    if (!is.null(out_loggle$error)){
      message("loggle ERROR:", out_loggle$error)
    }

    results[[j]]$loggle <- out_loggle
    rm(list = "out_loggle")
    gc()

    # save the trial
    time_delta <- round(as.numeric(Sys.time() - trial_start, units = "mins"))
    message("\nloggle trial ", j, " complete; ", time_delta, " minutes elapsed\n")
    save(results, file = filename)
  }

  # save the results
  save(results, file = filename)

  # trials for mgm
  functions <- c("eval_est", "sp.array", "tvmgm.eval")
  packages <- "mgm"
  doParallel::registerDoParallel(num_workers)
  results_mgm <- foreach(j = 1:n_trials, .export = functions,
                         .packages = packages)%dopar%
    {

      # ensure that loggle produced results for the current trial
      if (!is.null(results[[j]]$loggle$error)){
        return(NA)
      }

      # record the time the trial started
      trial_start <- Sys.time()

      # get the data
      data <- data_list[[j]]

      # mgm
      out_mgm <- tryCatch(tvmgm.eval(X = data$X,
                                     Z = data$Z,
                                     true = data$true_precision),
                          error = function(e) list(error = e))
      if (!is.null(out_mgm$error)){
        message("mgm ERROR:", out_mgm$error)
      }

      time_delta <- round(as.numeric(Sys.time() - trial_start, units = "mins"))
      message("\nmgm trial ", j, " complete; ", time_delta, " minutes elapsed\n")
      out_mgm
    }

  # add the mgm results to overall results
  for (j in 1:n_trials){
    results[[j]]$mgm <- results_mgm[[j]]
  }

  # save the results
  save(results, file = filename)

  # check if varbvs trials should be performed
  if ("varbvs" %in% names(results$trial1)){

    # trials for varbvs
    functions <- c("eval_est", "sp.array", "varbvs.eval")
    packages <- "varbvs"
    doParallel::registerDoParallel(num_workers)
    results_varbvs <- foreach(j = 1:n_trials, .export = functions,
                              .packages = packages)%dopar%
      {

        # ensure that both loggle and mgm produced results for the current trial
        if (!is.null(results[[j]]$mgm$error) | !is.null(results[[j]]$loggle$error)){
          return(NA)
        }

        # record the time the trial started
        trial_start <- Sys.time()

        # get the data
        data <- data_list[[j]]

        # varbvs
        out_varbvs <- tryCatch(varbvs.eval(X = data$X,
                                           Z = data$Z,
                                           true = data$true_precision),
                               error = function(e) list(error = e))
        if (!is.null(out_varbvs$error)){
          message("varbvs ERROR:", out_varbvs$error)
        }

        time_delta <- round(as.numeric(Sys.time() - trial_start, units = "mins"))
        message("\nvarbvs trial ", j, " complete; ", time_delta, " minutes elapsed\n")
        out_varbvs
      }

    # add the varbvs results to overall results
    for (j in 1:n_trials){
      results[[j]]$varbvs <- results_varbvs[[j]]
    }

    # save the results
    save(results, file = filename)
  }

  # trials for covdepGE

  # get number of available workers
  (num_workers <- parallel::detectCores() - 5)
  for (j in 1:n_trials){

    # ensure that both loggle, mgm, and varbvs produced results for the current trial
    if (!is.null(results[[j]]$mgm$error) | !is.null(results[[j]]$loggle$error) | !is.null(results[[j]]$varbvs$error)){
      return(NA)
    }

    # record the time the trial started
    trial_start <- Sys.time()

    # get the data
    data <- data_list[[j]]

    # covdepGE
    out_covdepGE <- tryCatch(covdepGE.eval(X = data$X,
                                           Z = data$Z,
                                           true = data$true_precision,
                                           n_workers = num_workers),
                             error = function(e) list(error = e))
    if (!is.null(out_covdepGE$error)){
      message("covdepGE ERROR:", out_covdepGE$error)
    }

    results[[j]]$covdepGE <- out_covdepGE
    rm(list = "out_covdepGE")
    gc()

    time_delta <- round(as.numeric(Sys.time() - trial_start, units = "mins"))
    message("\ncovdepGE trial ", j, " complete; ", time_delta, " minutes elapsed\n")
    save(results, file = filename)
  }

  # save the results and stop the cluster
  save(results, file = filename)
  doParallel::stopImplicitCluster()
}
