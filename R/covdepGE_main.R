## -----------------------------------------------------------------------------
#' @title covdepGE: Covariate Dependent Graph Estimation
#' @aliases covdepGE-method
#' @export
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' @description Model the conditional dependence structure of data as a function
#' of extraneous covariates as described in (1).
## -----------------------------ARGUMENTS---------------------------------------
#' @param data_mat \eqn{n} x \eqn{(p + 1)} matrix; data
#'
#' @param Z \eqn{n} x \eqn{p'} matrix; extraneous covariates
#'
#' @param tau positive numeric OR numeric vector of length \eqn{n} with positive
#' entries; bandwidth parameter. Greater values allow for more information to be
#' shared between individuals. Allows for global or individual-specific
#' specification. If `kde = T`, this argument is ignored. `0.1` by default
#'
#' @param kde logical; if `T`, use 2-step KDE methodology as described in
#' (2) to calculate individual-specific bandwidths in place of global bandwidth
#' parameter `tau`. `T` by default
#'
#' @param alpha numeric in \eqn{(0, 1)}; global initialization value for the
#' variational parameters `alpha_matrices` (approximates probabilities of
#' inclusion). `0.2` by default
#'
#' @param mu numeric; global initialization value for the variational parameters
#' `mu_matrices` (approximates posterior mean of regression coefficients). `0`
#' by default
#'
#' @param sigmasq_vec numeric vector of length `n_param` with positive
#' entries; candidate values of `sigmasq`, the error term variance. `0.5` by
#' default
#'
#' @param sigmabetasq_vec numeric vector of length `n_param` with positive
#' entries; candidate values of `sigmabeta_sq`, the slab variance. `NULL` by
#' default
#'
#' @param n_param positive integer; if `sigmabetasq_vec` is `NULL`, `n_param` is
#' the number of candidate `sigmabeta_sq` that will be auto-generated as:
#'
#' `sigmabetasq_vec <- exp(seq(log(var_max), log(var_min), length = n_param))`
#'
#' `8` by default
#'
#' @param pi_vec numeric vector of length `n_param` with entries in \eqn{(0, 1)};
#' candidate values of `pi`. `0.1` by default
#'
#' @param norm numeric in \eqn{[1, Inf]}; norm to use when calculating weights.
#' `Inf` results in infinity norm. `2` by default
#'
#' @param scale logical; if `T`, center and scale extraneous covariates to mean
#' 0, standard deviation 1 prior to calculating the weights. `T` by default
#'
#' @param tolerance positive numeric; end CAVI when the Frobenius norm of the
#' iteration-to-iteration change in the alpha matrix are within tolerance.
#' `1e-12` by default
#'
#' @param max_iter positive integer; if the tolerance criteria has not been met
#' by `max_iter_grid` iterations, end CAVI. `100` by default
#'
#' @param edge_threshold numeric in \eqn{(0, 1)}; when post-processing the
#' inclusion probabilities, an edge will be added to the graph if the
#' \eqn{(i, j)} edge has probability of inclusion greater than `edge_threshold`.
#' `0.5` by default
#'
#' @param sym_method character in \{`"mean"`, `"max"`, `"min"`\}; to symmetrize
#' the alpha matrices, the \eqn{i,j = j,i} entry is
#' `sym_method`\eqn{((i,j entry), (j,i entry))}. `"mean"` by default
#'
#' @param parallel logical; if `T`, grid search and CAVI for each variable will
#' be performed in parallel using `foreach`. Parallel backend may be registered
#' prior to making a call to `covdepGE`. If no active parallel backend can be
#' detected, then parallel backend will be automatically registered using
#' `doParallel::registerDoParallel(num_workers)`
#'
#' @param num_workers integer in \eqn{{1, 2,...,`parallel::detectCores()`}};
#' argument to `doParallel::registerDoParallel` if `parallel = T` and no
#' parallel backend is detected. `NULL` by default, which results in
#' `num_workers = floor(parallel::detectCores() / 2)`
#'
#' @param stop_cluster logical; if `T`, run `doParallel::stopImplicitCluster()`
#' after parallel exectution of CAVI for all variables has completed. This will
#' stop the cluster created implicitly by
#' `doParallel::registerDoParallel(num_workers)` and will shut down the unused
#' workers. Setting `F` is useful when making multiple calls to `covdepGE` with
#' `parallel = T`, as it avoids the overhead of creating a new cluster. `T` by
#' default
#'
#' @param warnings logical; if `T`, convergence and grid warnings will be
#' displayed. Convergence warnings occur when the tolerance exit condition has
#' not been met by `max_iter_grid` or `max_iter_final` iterations. Grid warnings
#' occur when, for either `sigmabetasq_vec` or `pi_vec`, the grid is longer than
#' 2 candidates, and the final CAVI selects a candidate value on the grid
#' boundary. `T` by default
#'
#' @param CS logical; if `T`, `pi_vec` and `sigma_sq` will be selected according
#' to Carbonetto-Stephens. `F` by default
## -----------------------------RETURNS-----------------------------------------
#' @return Returns `list` with the following values:
#'
#' 1. `graphs`: `list` of \eqn{n} \eqn{(p + 1)} x \eqn{(p + 1)} matrices; the
#' \eqn{l}-th matrix is the adjacency matrix for the \eqn{l}-th individual
#' (obtained from `inclusion_probs` according to `edge_threshold`)
#'
#' 2. `inclusion_probs`: list of \eqn{n} \eqn{(p + 1)} x \eqn{(p + 1)} matrices;
#' the \eqn{l}-th matrix is a symmetric matrix of inclusion probabilities for
#' the \eqn{l}-th individual (obtained by symmetrizing the `alpha_matrices`
#' according to `sym_method`)
#'
#' 3. `alpha_matrices`: list of \eqn{n} \eqn{(p + 1)} x \eqn{(p + 1)} matrices;
#' the \eqn{l}-th matrix is an asymmetric matrix of inclusion probabilities for
#' the \eqn{l}-th individual
#'
#' 4. `unique_graphs`: list of \eqn{g} lists; \eqn{g} is the number of
#' unique graphs. The \eqn{v}-th list has 3 values:
#' - `graph`: \eqn{(p + 1)} x \eqn{(p + 1)} matrix; the adjacency matrix for the
#'  \eqn{v}-th graph
#' - `individuals`: vector with entries in \eqn{{1,...,n}}; the individual
#' indices corresponding to the \eqn{v}-th graph
#' - `individuals_summary`: character; summarizes the individual indices in
#' `individuals`
#'
#' 5. `CAVI_details`: list of \eqn{(p + 1)} lists; the \eqn{j}-th list
#' corresponds to the \eqn{j}-th variable and contains the following values:
#'  - `sigmasq`, `sigmabeta_sq`, `pi`: numerics; the values of the
#'  hyperparameters that maximized the ELBO for the j-th variable
#'  - `ELBO`: numeric; the maximum value of ELBO for the final CAVI
#'  - `converged_iter`: numeric; the number of iterations to attain convergence
#'  for the final CAVI
#'  - `hyperparameters`: `n_param` x 4 matrix; each of the hyperparameter grid
##  points with the resulting ELBO
#'
#' 6. `model_details`: list with the following values:
#'  - `elapsed`: timediff; the amount of time to fit the model
#'  - `n`: integer; sample size
#'  - `p`: integer; the number of variables in the data minus one
#'  - `ELBO`: numeric; the total ELBO summed across the \eqn{p + 1} final models
#'  - `num_unique`: integer; the number of unique conditional dependence
#'  structures identified
#'  - `final_DNC`: integer; the number of variables for which the final CAVI did
#'  not attain convergence within `max_iter_final` iterations
#'  - `num_candidates`: integer; the number of grid points in the grid search
#'
#' 7. `weights`: \eqn{n} x \eqn{n} matrix; the \eqn{i, j} entry is the weighting
#' of the \eqn{i}-th individual with respect to the \eqn{j}-th individual using
#' the \eqn{j}-th individual's bandwidth
#'
#' 8. `bandwidths`: vector of length \eqn{n}; individual-specific bandwidths
#'
#' 9. `arguments`: vector; argument values passed to the current call to
#' `covdepGE`
## -----------------------------EXAMPLES----------------------------------------
#' @examples
#' set.seed(1)
#' n <- 100
#' p <- 4
#'
#' # generate the extraneous covariate
#' Z_neg <- sort(runif(n / 2) * -1)
#' Z_pos <- sort(runif(n / 2))
#' Z <- c(Z_neg, Z_pos)
#' summary(Z)
#'
#' # create true covariance structure for 2 groups: positive Z and negative Z
#' true_graph_pos <- true_graph_neg <- matrix(0, p + 1, p + 1)
#' true_graph_pos[1, 2] <- true_graph_pos[2, 1] <- 1
#' true_graph_neg[1, 3] <- true_graph_neg[3, 1] <- 1
#'
#' # visualize the true covariance structures
#' (gg_adjMat(true_graph_neg) +
#'     ggplot2::ggtitle("True graph for individuals with negative Z"))
#' (gg_adjMat(true_graph_pos, color1 = "steelblue") +
#'     ggplot2::ggtitle("True graph for individuals with positive Z"))
#'
#' # generate the covariance matrices as a function of Z
#' sigma_mats_neg <- lapply(Z_neg, function(z) z * true_graph_neg + diag(p + 1))
#' sigma_mats_pos <- lapply(Z_pos, function(z) z * true_graph_pos + diag(p + 1))
#' sigma_mats <- c(sigma_mats_neg, sigma_mats_pos)
#'
#' # generate the data using the covariance matrices
#' data_mat <- t(sapply(sigma_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))
#'
#' # visualize the sample correlation
#' gg_adjMat(abs(cor(data_mat[1:(n / 2), ])) - diag(p + 1))
#' gg_adjMat(abs(cor(data_mat[(n / 2 + 1):n, ])) - diag(p + 1),
#'           color1 = "dodgerblue")
#'
#' # estimate the covariance structure
#' out <- covdepGE(data_mat, Z)
#'
#' # analyze results
#' gg_adjMat(out, 1)
#' gg_adjMat(out, 50, color1 = "tomato")
#' gg_adjMat(out, 54, color1 = "steelblue")
#' gg_adjMat(out, 100, color1 = "dodgerblue")
#'
#' gg_inclusionCurve(out, 1, 2)
#' gg_inclusionCurve(out, 1, 3, point_color = "dodgerblue")
## -----------------------------REFERENCES--------------------------------------
#' @references
#' 1. Dasgupta S., Ghosh P., Pati D., Mallick B., *An approximate Bayesian
#' approach to covariate dependent graphical modeling*, 2021
#'
#' 2. Dasgupta S., Pati D., Srivastava A., *A Two-Step Geometric Framework For
#' Density Modeling*, Statistica Sinica, 2020
## -----------------------------------------------------------------------------
covdepGE <- function(data, Z, alpha = 0.2, mu = 0, ssq = NULL, ssq_p = NULL,
                     ssq_q = NULL, sbsq = NULL, sbsq_p = NULL, sbsq_q = NULL,
                     pip = NULL, pip_p = NULL, pip_q = NULL, nssq = 5, nsbsq = 5,
                     npip = 5, ssq_upper_mult = 4, var_lower = 1e-3,
                     grid_search = T, tau = 0.1, kde = T, norm = 2, scale = T,
                     elbo_tol = 1e-1, alpha_tol = 1e-5, max_iter = 100,
                     edge_threshold = 0.5, sym_method = "mean", parallel = F,
                     num_workers = NULL, stop_cluster = T, prog_bar = T,
                     warnings = T, CS = F, R = F){

  start_time <- Sys.time()

  # run compatibility checks
  # covdepGE_checks(data, Z, tau, kde, alpha, mu, ssq, sbsq, pi_vec, norm, scale, tolerance, max_iter, edge_threshold,
  #                 sym_method, parallel, num_workers, stop_cluster, warnings)

  # ensure that data_mat and Z are matrices
  data <- as.matrix(data)
  Z <- as.matrix(Z)

  # get sample size and number of parameters
  n <- nrow(data); p <- ncol(data) - 1

  # if the hyperparameters, importance densities, or prior densities are
  # m-vectors, change to m x p matrix
  if (is.vector(ssq)) ssq <- matrix(ssq, length(ssq), p + 1)
  if (is.vector(sbsq)) sbsq <- matrix(sbsq, length(sbsq), p + 1)
  if (is.vector(pip)) pip <- matrix(pip, length(pip), p + 1)
  if (is.vector(ssq_p)) ssq_p <- matrix(ssq_p, length(ssq_p), p + 1)
  if (is.vector(sbsq_p)) sbsq_p <- matrix(sbsq_p, length(sbsq_p), p + 1)
  if (is.vector(pip_p)) pip_p <- matrix(pip_p, length(pip_p), p + 1)
  if (is.vector(ssq_q)) ssq_q <- matrix(ssq_q, length(ssq_q), p + 1)
  if (is.vector(sbsq_q)) sbsq_q <- matrix(sbsq_q, length(sbsq_q), p + 1)
  if (is.vector(pip_q)) pip_q <- matrix(pip_q, length(pip_q), p + 1)

  # if the covariates should be centered and scaled, do so ([ , ] for attributes)
  if (scale) Z <- matrix(scale(Z)[ , ], n)

  # get weights
  D <- get_weights(Z, norm, kde, tau)
  bandwidths <- D$bandwidths
  D <- D$D

  # main loop over the predictors

  # check if CAVIs are to be parallelized
  if (parallel){

    # check to see if parallel backend has been registered
    registered <- tryCatch(
      {
        # return true if parallel backend is registered with more than 1 worker
        foreach::`%dopar%`(foreach::foreach(NULL), NULL)
        foreach::getDoParRegistered()
      },

      # return false if error
      error = function(msg) F,

      # return false if warning
      warning = function(msg) F)

    if (registered){
      if (warnings) message(paste("Detected", foreach::getDoParWorkers(),
                                  "workers"))
    }else{

      # otherwise, register parallel backend

      # if num_workers has not been provided, get the number of workers
      if (is.null(num_workers)){
        num_workers <- floor(parallel::detectCores() / 2)
      }

      # perform registration
      if (warnings) warning(paste(
        "No registered workers detected; registering doParallel with",
        num_workers, "workers"))
      doParallel::registerDoParallel(cores = num_workers)
    }

    # attempt to execute the variational update in parallel
    res <- tryCatch(
      {
        foreach::`%dopar%`(
          foreach::foreach(resp_index = 1:(p + 1), .packages = "covdepGE"),
          {

            # Set variable number `resp_index` as the response
            y <- data[, resp_index]

            # Set the remaining p variables as predictors
            X <- data[, -resp_index]

            # perform the grid search and final CAVI; save the results to res
            cavi_search(X, Z, D, y, alpha, mu, ssq, ssq_p, ssq_q, sbsq, sbsq_p,
                        sbsq_q, pip, pip_p, pip_q, nssq, nsbsq, npip,
                        ssq_upper_mult, var_lower, grid_search, elbo_tol,
                        alpha_tol, max_iter, warnings, resp_index, CS, R)
            }
          )
      },

      # if parallel execution did not finish successfully, display an error
      error = function(msg) stop(paste(
        "Parallel execution failed; error message: ", msg))
    )

    # shut down the workers if desired
    if(stop_cluster) doParallel::stopImplicitCluster()

  }else{

    # otherwise, CAVI will be executed sequentially

    # instantiate the progress bar
    if (prog_bar) pb <- utils::txtProgressBar(0, p + 1, style = 3)

    # list to store each of the results from cavi_search
    res <- vector("list", p + 1)

    for (resp_index in 1:(p + 1)) {

      # Set variable number `resp_index` as the response
      y <- data[, resp_index]

      # Set the remaining p variables as predictors
      X <- data[, -resp_index]

      # perform the grid search and final CAVI; save the results to res
      res[[resp_index]] <- cavi_search(X, Z, D, y, alpha, mu, ssq, ssq_p, ssq_q,
                                       sbsq, sbsq_p, sbsq_q, pip, pip_p, pip_q,
                                       nssq, nsbsq, npip, ssq_upper_mult,
                                       var_lower, grid_search, elbo_tol,
                                       alpha_tol, max_iter, warnings, resp_index,
                                       CS, R)

      # update the progress bar
      if (prog_bar) utils::setTxtProgressBar(pb, resp_index)
    }

    # close the progress bar
    if (prog_bar) close(pb)
  }

  # gather the cavi_details lists, hyperparameter_details lists, and
  # alpha_matrix matrices into their own lists/ vectors
  cavi_details <- lapply(res, `[[`, "cavi_details")
  hp <- lapply(res, `[[`, "hyperparameters")
  grid_sz <- hp[[1]]$grid_sz
  names(cavi_details) <- names(hp) <- paste0("Variable ", 1:length(cavi_details))
  alpha_matrices <- lapply(res, `[[`, "alpha_matrix")
  mu_matrices <- lapply(res, `[[`, "mu_matrix")

  # gather progress of ELBO and alpha
  progress <- lapply(res, `[[`, "progress")
  names(progress) <- names(hp)

  # organize the hyperparameters into matrices
  ssq <- sapply(hp, `[[`, "ssq")
  ssq_cands <- sapply(hp, `[[`, "ssq_cands")
  sbsq <- sapply(hp, `[[`, "sbsq")
  sbsq_cands <- sapply(hp, `[[`, "sbsq_cands")
  pip <- sapply(hp, `[[`, "pip")
  pip_cands <- sapply(hp, `[[`, "pip_cands")
  hyperparameters <- list(ssq = ssq, sbsq = sbsq, pip = pip, ssq_cands = ssq_cands,
                          sbsq_cands = sbsq_cands, pip_cands = pip_cands)

  # find the number of final CAVIs that converged
  converged <- sapply(cavi_details, `[[`, "converged")

  # if any of the CAVI for a variable did not converge, display a warning
  if (warnings & sum(!converged) > 0){

    sapply(which(!converged), function(var_index) warning(paste0(
      "Variable ", var_index, ": final CAVI did not converge in ",
      max_iter, " iterations")))
  }

  # Creating the graphs:
  # transform p + 1 n by n matrices to n p + 1 by p + 1 matrices using
  # alpha_matrices
  # the j, k entry in the l-th matrix is the probability of inclusion of an edge
  # between the j, k variables for the l-th individual
  incl_probs <- replicate(n, matrix(0, p + 1, p + 1), simplify = F)

  # iterate over the p matrices
  for (j in 1:(p + 1)){

    # fix the j-th alpha matrix
    alpha_mat_j <- alpha_matrices[[j]]

    # iterate over the rows of alpha_mat_j
    for (l in 1:n){

      # the j-th row of the l-th individual's graph is the l-th row of
      # alpha_mat_j with a 0 in the j-th position
      incl_probs[[l]][j, -j] <- alpha_mat_j[l,]
    }
  }

  # save the asymmetric matrices
  incl_probs_asym <- incl_probs

  # symmetrize the inclusion matrices according to the symmetrization method
  if (sym_method == "mean"){

    # take the mean of (i,j), (j,i) entries to symmetrize
    incl_probs <- lapply(incl_probs, function(mat) (mat + t(mat)) / 2)
  }else if (sym_method == "min"){

    # take the min of (i,j), (j,i) entries to symmetrize
    incl_probs <- lapply(incl_probs, function(mat) pmin(mat, t(mat)))
  }else if (sym_method == "max"){

    # take the max of (i,j), (j,i) entries to symmetrize
    incl_probs <- lapply(incl_probs, function(mat) pmax(mat, t(mat)))
  }

  # using the symmetrized graphs, if the probability of an edge is greater than
  # edge_threshold, denote an edge by 1; otherwise, 0
  graphs <- lapply(incl_probs, function(mat) (mat > edge_threshold) * 1)

  # find the unique graphs
  unique_graphs <- unique(graphs)

  # find the individuals corresponding to each of the unique graphs
  indv_graphs <- lapply(unique_graphs, function(unique_graph)
    which(sapply(graphs, function(indv_graph) identical(indv_graph, unique_graph))))

  # for each unique graph, create a summary of the individuals corresponding to
  # that graph
  indv_graphs_sum <- lapply(indv_graphs, function(indv_graph) paste0(sapply(
    split(sort(indv_graph), cumsum(c(1, diff(sort(indv_graph)) != 1))), function(
      idx_seq) ifelse(length(idx_seq) > 2, paste0(min(idx_seq), ",...,", max(
        idx_seq)), paste0(idx_seq, collapse = ","))), collapse = ","))

  # create a nested list where the j-th inner list has three values; the j-th
  # unique graph, the individuals corresponding to that graph, and a summary of
  # the individuals corresponding to that graph
  unique_graphs <- lapply(1:length(unique_graphs), function(gr_idx)
    list(graph = unique_graphs[[gr_idx]], individuals = indv_graphs[[gr_idx]],
         individuals_summary = indv_graphs_sum[[gr_idx]]))
  names(unique_graphs) <- paste0("graph", 1:length(unique_graphs))

  # calculate the total ELBO
  total_elbo <- sum(sapply(cavi_details, `[[`, "ELBO"))

  # create a list for model details
  model_details <- list(elapsed = NA, n = n, p = p, ELBO = total_elbo,
                        num_unique = length(unique_graphs),
                        converged = sum(converged), grid_size = grid_sz)

  # create a named vector to return the function arguments
  args <- c(kde = kde, norm = norm, scale = scale, elbo_tol = elbo_tol,
            alpha_tol = alpha_tol, edge_threshold = edge_threshold,
            sym_method = sym_method, parallel = parallel)

  # record the elapsed time and add it to the model details
  model_details[["elapsed"]] <- Sys.time() - start_time

  # define the list of return values
  ret <- list(graphs = graphs, inclusion_probs = incl_probs,
              alpha_matrices = incl_probs_asym, unique_graphs = unique_graphs,
              mu_matrices = mu_matrices, hyperparameters = hyperparameters,
              CAVI_details = cavi_details, model_details = model_details,
              weights = D, bandwidths = bandwidths, arguments = args,
              progress = progress)

  # define the class of the return values
  class(ret) <- c("covdepGE", "list")

  return(ret)
}

## -----------------------------------------------------------------------------
#' @title print.covdepGE
#' @export
#' @rdname covdepGE
## -----------------------------DESCRIPTION-------------------------------------
## S3 method for printing an object of class `covdepGE`
## -----------------------------ARGUMENTS---------------------------------------
#' @param x object of class covdepGE; the return of the covdepGE function
#' @param ... additional arguments will be ignored
## -----------------------------------------------------------------------------
print.covdepGE <- function(x, ...){

  cat("                      Covariate Dependent Graphical Model\n\n")

  spc <- 80

  with(x$model_details,
       {
         # print ELBO and number of unique graphs
         elbo_str <- paste0("Model ELBO: ", round(ELBO, 2))
         cat(sprintf(paste0("%-s%", spc - nchar(elbo_str), "s"), elbo_str,
                     paste0("Unique conditional dependence structures: ",
                            num_unique, "\n")))

         # print data dimensions and the grid size
         data_dim <- paste0("n: ", n, ", variables: ", p + 1)
         cat(sprintf(paste0("%-s%", spc - nchar(data_dim), "s"), data_dim,
                     paste0("Hyperparameter grid size: ", grid_size,
                            " points\n")))

         # print the number of converged final CAVIs
         cat("CAVI converged for ", converged, "/", p + 1,
             " variables\n\n", sep = "")

         # print time to fit
         time_units <- attr(elapsed, "units")
         duration <- as.numeric(elapsed)
         cat("Model fit completed in", round(duration, 3), time_units, "\n")
       }
  )
}
