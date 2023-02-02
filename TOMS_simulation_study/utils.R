library(doParallel)

# function for visualizing a list of graphs
graph_viz <- function(graphs){

  # find the unique graphs
  unique_graphs <- unique(graphs)

  # create a list where the j-th element is the j-th unique graph and the
  # indices of the observations corresponding to this graph
  unique_sum <- vector("list", length(unique_graphs))
  names(unique_sum) <- paste0("graph", 1:length(unique_graphs))

  # iterate over each of the unique graphs
  for (j in 1:length(unique_graphs)){

    # fix the unique graph
    graph <- unique_graphs[[j]]

    # find indices of the observations corresponding to this graph
    graph_inds <- which(sapply(graphs, identical, graph))

    # split up the contiguous subsequences of these indices
    cont_inds <- split(sort(graph_inds), cumsum(c(1, diff(sort(graph_inds))
                                                  != 1)))

    # create a character summary for each of the contiguous sequences
    inds_sum <- sapply(cont_inds, function(idx_seq) ifelse(length(
      idx_seq) > 3, paste0(min(idx_seq), ",...,", max(idx_seq)),
      paste0(idx_seq, collapse = ",")))

    # combine the summary
    inds_sum <- paste0(inds_sum, collapse = ",")

    # add the graph, indices, and summary to the unique graphs summary list
    unique_sum[[j]] <- list(graph = graph, indices = graph_inds,
                            ind_sum = inds_sum)
  }

  # create the titles for the plots
  titles <- paste("Graph", 1:length(unique_graphs))
  obs_sum <- sapply(unique_sum, `[[`, "ind_sum")
  titles <- paste0(titles, ", observations ", obs_sum)

  # create a visualization for each graph and store it in a list
  graph_viz <- lapply(1:length(unique_graphs), function(gr_idx) matViz(
    unique_graphs[[gr_idx]], color2 = "#500000") + ggplot2::ggtitle(
      titles[gr_idx]))

  graph_viz
}

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

# function for sorting based on bi-variate Z
sort_Z <- function(Z){

  sort_inds <- 1

  for (j in 2:nrow(Z)){

    # get index of the last observation and corresponding Z
    curr_ind <- sort_inds[length(sort_inds)]
    curr_Z <- matrix(Z[curr_ind , ], nrow(Z), 2, T)

    # get norm between current Z and all others; set the norm for those that have
    # already been sorted to Inf so that they are not sorted again
    norms <- rowSums((Z - curr_Z)^2)
    norms[sort_inds] <- Inf

    # choose the minimum norm as the next sorted observation
    sort_inds <- c(sort_inds, which.min(norms))
  }

  # plot(Z[sort_inds, ])
  # lines(Z[sort_inds, ])
  return(sort_inds)
}

# function for evaluating estimated graphs given ground truth
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
trials <- function(data_list, results, filename, skips, hp_method, max_iter_grid){

  # save sample data to results
  results$sample_data <- dim(data_list[[1]]$X)

  # get number of available workers and trials
  num_workers <- min(10, parallel::detectCores() - 5)
  n_trials <- length(data_list)

  # check if mgm trials should be performed
  if ("JGL" %in% names(results$trial1) & !("JGL" %in% skips)){

    # trials for JGL
    functions <- c("aic_JGL", "eval_est", "JGL.eval", "sp.array")
    packages <- c("JGL", "mclust")
    num_workers <- min(25, parallel::detectCores())
    doParallel::registerDoParallel(num_workers)
    results_JGL <- foreach(j = 1:n_trials, .export = functions,
                           .packages = packages)%dopar%
      {

        # record the time the trial started
        trial_start <- Sys.time()

        # get the data
        data <- data_list[[j]]

        # JGL
        out_JGL <- tryCatch(JGL.eval(X = data$X,
                                     Z = data$Z,
                                     true = data$true_precision),
                            error = function(e) list(error = e))
        if (!is.null(out_JGL$error)){
          message("JGL ERROR:", out_JGL$error)
        }

        time_delta <- round(as.numeric(Sys.time() - trial_start, units = "mins"))
        message("\nJGL trial ", j, " complete; ", time_delta, " minutes elapsed\n")
        out_JGL
      }

    # add JGL results to overall results
    for (j in 1:n_trials){
      results[[j]]$JGL <- results_JGL[[j]]
    }

    # save the results
    save(results, file = filename)

    message(paste("JGL finished", Sys.time()))
  }

  # check if mgm trials should be performed
  if ("mgm" %in% names(results$trial1) & !("mgm" %in% skips)){

    # trials for mgm
    functions <- c("eval_est", "sort_Z", "sp.array", "tvmgm.eval")
    packages <- "mgm"
    num_workers <- min(25, parallel::detectCores())
    doParallel::registerDoParallel(num_workers)
    results_mgm <- foreach(j = 1:n_trials, .export = functions,
                           .packages = packages)%dopar%
      {

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

    message(paste("mgm finished", Sys.time()))
  }

  doParallel::stopImplicitCluster()

  # check if covdepGE trials should be performed
  if (!("covdepGE" %in% skips)){

    # trials for covdepGE

    # get number of available workers
    (num_workers <- parallel::detectCores() - 5)
    doParallel::registerDoParallel(num_workers)
    for (j in 1:n_trials){

      # record the time the trial started
      trial_start <- Sys.time()

      # get the data
      data <- data_list[[j]]

      # covdepGE
      out_covdepGE <- tryCatch(covdepGE.eval(X = data$X,
                                             Z = data$Z,
                                             hp_method = hp_method,
                                             true = data$true_precision,
                                             n_workers = num_workers,
                                             max_iter_grid = max_iter_grid),
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

    message(paste("covdepGE finished", Sys.time()))
  }

  # check if covdepGE_sortZ trials should be performed
  if ("covdepGE_sortZ" %in% names(results$trial1) & !("covdepGE_sortZ" %in% skips)){
    for (j in 1:n_trials){

      # record the time the trial started
      trial_start <- Sys.time()

      # get the data
      data <- data_list[[j]]

      # covdepGE

      # sort observations in X and ground truth
      sort_inds <- sort_Z(data$Z)
      data$X <- data$X[sort_inds, ]
      data$true_precision <- data$true_precision[sort_inds]
      data$Z <- 1:nrow(data$X)
      out_covdepGE <- tryCatch(covdepGE.eval(X = data$X,
                                             Z = data$Z,
                                             hp_method = hp_method,
                                             true = data$true_precision,
                                             n_workers = num_workers),
                               error = function(e) list(error = e))
      if (!is.null(out_covdepGE$error)){
        message("covdepGE ERROR:", out_covdepGE$error)
      }

      results[[j]]$covdepGE_sortZ <- out_covdepGE
      rm(list = "out_covdepGE")
      gc()

      time_delta <- round(as.numeric(Sys.time() - trial_start, units = "mins"))
      message("\ncovdepGE_sortZ trial ", j, " complete; ", time_delta, " minutes elapsed\n")
      save(results, file = filename)
    }
  }

  # save the results and stop the cluster
  save(results, file = filename)
  doParallel::stopImplicitCluster()
}
