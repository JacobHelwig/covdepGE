## -----------------------------------------------------------------------------
#' @title covdepGE: Covariate Dependent Graph Estimation
#' @export
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' @description Model the conditional dependence structure of data as a function
#' of extraneous covariates as described in (1). The final model uses the
#' candidate pair of `pi` and `sigmabeta_sq` that maximize the ELBO for each
#' variable
## -----------------------------ARGUMENTS---------------------------------------
#' @param data_mat \eqn{n} x \eqn{(p + 1)} matrix; data
#'
#' @param Z \eqn{n} x \eqn{p'} matrix; extraneous covariates
#'
#' @param tau scalar in \eqn{(0, Inf)} OR \eqn{n} x \eqn{1} vector, entries in
#' \eqn{(0, Inf)}; bandwidth parameter. Greater values allow for more
#' information to be shared between individuals. Allows for global or
#' individual-specific specification. If `kde = T`, this argument is ignored.
#' `0.1` by default
#'
#' @param kde logical scalar; if `T`, use 2-step KDE methodology as described in
#' (2) to calculate individual-specific bandwidths in place of global bandwidth
#' parameter `tau.` `T` by default
#'
#' @param alpha scalar in \eqn{[0, 1]}; global initialization value for the
#' variational parameters `alpha_matrices` (approximates probabilities of
#' inclusion). `0.2` by default
#'
#' @param mu scalar; global initialization value for the variational parameters
#' `mu_matrices` (approximates regression coefficients). 0 by default
#'
#' @param sigmasq scalar in \eqn{(0, Inf)}; Error term variance for
#' spike-and-slab. Algorithm scales this value by individual-specific weights.
#' `0.5` by default
#'
#' @param sigmabetasq_vec `n_sigma` x \eqn{1} vector, entries in\eqn{ (0, Inf)};
#' candidate values of `sigmabeta_sq`, the slab variance. `NULL` by default
#'
#' @param var_min scalar in \eqn{(0, Inf)}; if `sigmabetasq_vec` is `NULL`,
#'  `var_min` is the lower bound of the auto-generated `sigmabetasq_vec`. `0.01`
#'  by default
#'
#' @param var_max scalar in \eqn{(0, Inf)}; if `sigmabetasq_vec` is `NULL`,
#' `var_max` is the upper bound of the auto-generated `sigmabetasq_vec`. `10` by
#' default
#'
#' @param n_sigma scalar in \eqn{{1, 2,...}}; if `sigmabetasq_vec` is `NULL`,
#'  `n_sigma` is the number of candidate `sigmabetasq` that will be
#'  auto-generated as:
#'
#' `sigmabetasq_vec <- exp(seq(log(var_max), log(var_min), length = n_sigma))`
#'
#' `8` by default
#'
#' @param pi_vec `n_pi` x \eqn{1} vector, entries in \eqn{[0, 1]}; candidate
#' values of `pi`. `0.1` by default
#'
#' @param norm scalar in \eqn{[1, Inf]}; norm to use when calculating weights.
#' `Inf` results in infinity norm. `2` by default
#'
#' @param scale logical scalar; if `T`, center and scale extraneous covariates
#' to mean 0, standard deviation 1 prior to calculating the weights. `T` by
#' default
#'
#' @param tolerance scalar in \eqn{(0, Inf)}; end the variational update loop
#'  when the square root of the sum of squared changes to the elements of the
#'  alpha matrix are within tolerance. `1e-12` by default
#'
#' @param max_iter scalar in \eqn{{1, 2,...}}; if the tolerance criteria has not
#' been met by `max_iter` iterations, end the variational update loop. `1e4` by
#' default
#'
#' @param edge_threshold scalar in \eqn{(0, 1)}; when post-processing the
#'  inclusion probabilities, an edge will be added to the graph if the
#'  \eqn{(i, j)} edge has probability of inclusion greater than `edge_threshold`.
#'  `0.5` by default
#'
#' @param sym_method character scalar in \{`"mean"`, `"max"`, `"min"`\}; to
#' symmetrize the alpha matrices, the \eqn{i,j = j,i} entry is
#' `sym_method`\eqn{((i,j entry), (j,i entry))}. `"mean"` by default
#'
#' @param parallel logical scalar; if `T`, the variational updates for each
#' response will be performed in parallel using `foreach`. Parallel backend may
#' be registered prior to making a call to `covdepGE`. If no active parallel
#' backend can be detected, then parallel backend will be automatically
#' registered using `doParallel::registerDoParallel(num_workers)`
#'
#' @param num_workers scalar in {1, 2,...,parallel::detectCores()}; number of
#' workers to pass to `doParallel::registerDoParallel` if `parallel = T` and no
#' parallel backend is detected. `NULL` by default, which results in
#' `num_workers = floor(parallel::detectCores() / 2)`
#'
#' @param stop_cluster: logical scalar; if `T`, run
#' `doParallel::stopImplicitCluster()` after parallel exectution of variational
#' updates has completed. This will stop the cluster created implicitly by
#' `doParallel::registerDoParallel(num_workers)` and will shut down the unused
#' workers. Setting `F` is useful when making multiple calls to `covdepGE` with
#' `parallel = T`, as it avoids the overhead of creating a new cluster. `T` by
#' default
#'
#' @param monitor_final_elbo logical scalar; if `T`, the ELBO history for the
#' final model will be returned. `F` by default
#'
#' @param monitor_cand_elbo logical scalar; if `T`, the ELBO history for each of
#' the candidate models that do not attain convergence within `max_iter`
#' iterations will be returned
#'
#' @param monitor_period scalar in \eqn{{1, 2,..., max_iter}}; the periodicity
#' with which the ELBO is recorded if `monitor_final_elbo` or `monitor_cand_elbo`
#' is `T`. `1` by default
#'
#' @param print_time logical scalar; if `T`, function run time is printed. `F`
#' by default
#'
#' @param warnings logical scalar; if `T`, convergence and grid warnings will be
#'  displayed. Convergence warnings occur when the tolerance exit condition has
#'  not been met by `max_iter` iterations. Grid warnings occur when, for either
#'  `sigmabetasq_vec` or `pi_vec`, the grid is longer than 2 candidates, and the
#'  final model selects a candidate value on the grid boundary. `T` by default
#'
#' @param CS logical scalar; if `T`, `pi_vec` and `sigma_sq` will be scalars
#' selected according to Carbonetto-Stephens. `F` by default
## -----------------------------RETURNS-----------------------------------------
#' @return Returns `list` with the following values:
#'
#' 1. `graphs`: `list` of \eqn{n} \eqn{(p + 1)} x \eqn{(p + 1)} matrices; the
#' \eqn{l}-th matrix is the adjacency matrix for the \eqn{l}-th individual
#' (obtained from `inclusion_probs` according to `edge_threshold`)
#'
#' 2. `inclusion_probs`: list of \eqn{n} \eqn{(p + 1)} x \eqn{(p + 1)} matrices;
#'  the \eqn{l}-th matrix is a symmetric matrix of inclusion probabilities for
#'  the \eqn{l}-th individual (obtained by symmetrizing the `alpha_matrices`
#'  according to `sym_method`)
#'
#' 3. `alpha_matrices`: list of \eqn{n} \eqn{(p + 1)} x \eqn{(p + 1)} matrices;
#' the \eqn{l}-th matrix is an asymmetric matrix of inclusion probabilities for
#'  the \eqn{l}-th individual
#'
#' 4. `unique_graphs`: list of \eqn{g} lists; \eqn{g} is the number of
#' unique graphs. The \eqn{v}-th list has 3 values:
#' - `graph`: \eqn{(p + 1)} x \eqn{(p + 1)} matrix; the adjacency matrix for the
#'  \eqn{v}-th graph
#' - `individuals`: \eqn{u} x \eqn{1} vector; \eqn{u} is the number of
#'  individuals corresponding to the \eqn{v}-th graph. The elements of this
#'  vector are the individual indices corresponding to the \eqn{v}-th graph
#' - `individuals_summary`: character scalar; summarizes the individual indices
#' in `individuals`
#'
#' 5. `VB_details`: list of \eqn{(p + 1)} `lists`; the \eqn{j}-th list corresponds to
#'  the \eqn{j}-th predictor and contains 6 values:
#'  - `sigmabeta_sq`, `pi`: scalars; the final values of `pi`
#'  and `sigmabeta_sq` that maximized ELBO over all individuals with the
#'  \eqn{j}-th predictor fixed as the response
#'  - `ELBO`: scalar; the maximum value of ELBO for the final model
#'  - `converged_iter`: scalar; the number of iterations to attain convergence
#'  for the final model
#'  - `ELBO_history`: vector; ELBO history by iteration for the final model. If
#'  `monitor_final_elbo` is `F`, then this value will be `NULL`
#'  - `non_converged`: matrix; each row corresponds to the ELBO history for each
#'  of the candidate models that did not converge. If `monitor_cand_elbo` is
#'  `F`, then the ELBO history is omitted, and only the non-convergent
#'  `sigmabeta_sq` and `pi` values are provided. If all pairs resulted in
#'  convergence, then this value is `NULL`
#'
#' 6. `weights`: \eqn{n} x \eqn{n} matrix; the \eqn{j, i} entry is the weighting
#'  of the \eqn{j}-th individual with respect to the \eqn{i}-th individual using
#'  the \eqn{i}-th individual's bandwidth
#'
#' 7. `bandwidths`: \eqn{n} x \eqn{1} vector; individual-specific bandwidths
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
covdepGE <- function(data_mat, Z, tau = 0.1, kde = T, alpha = 0.2, mu = 0,
                     sigmasq = 0.5, sigmabetasq_vec = NULL, var_min = 0.01,
                     var_max = 10, n_sigma = 8, pi_vec = 0.1, norm = 2,
                     scale = T, tolerance = 1e-12, max_iter = 1e4,
                     edge_threshold = 0.5, sym_method = "mean", parallel = F,
                     num_workers = NULL, stop_cluster = T, monitor_final_elbo = F,
                     monitor_cand_elbo = F, monitor_period = 1, print_time = F,
                     warnings = T, CS = F){

  start_time <- Sys.time()

  # run compatibility checks
  covdepGE_checks(data_mat, Z, tau, kde, alpha, mu, sigmasq, sigmabetasq_vec,
                  var_min, var_max, n_sigma, pi_vec, norm, scale, tolerance,
                  max_iter, edge_threshold, sym_method, parallel, num_workers,
                  stop_cluster, monitor_final_elbo, monitor_cand_elbo,
                  monitor_period, print_time, warnings)

  # ensure that data_mat and Z are matrices
  data_mat <- as.matrix(data_mat)
  Z <- as.matrix(Z)

  # get sample size and number of parameters
  n <- nrow(data_mat); p <- ncol(data_mat) - 1

  # if the covariates should be centered and scaled, do so
  if (scale) Z <- matrix(scale(Z)[ , ], n)

  # get weights
  D <- get_weights(Z, norm, kde, tau)
  bandwidths <- D$bandwidths
  D <- D$D

  # if sigmabetasq_vec is NULL, instantiate the grid
  if(is.null(sigmabetasq_vec)){
    sigmabetasq_vec <- exp(seq(log(var_max), log(var_min), length = n_sigma))
  }

  # main loop over the predictors

  # check if the variational updates are to be parallelized
  if (parallel){

    # check to see if parallel backend has been registered and is active
    if (foreach::getDoParRegistered() & nrow(showConnections()) != 0){
      if (warnings) message(paste("Detected", foreach::getDoParWorkers(),
                                  "registered workers on an active cluster"))
    }else{

      # otherwise, register parallel backend

      # if num_workers has not been provided, get the number of workers
      if (is.null(num_workers)){
        num_workers <- floor(parallel::detectCores() / 2)
      }

      # registration
      if (warnings) warning(paste("Active parallel backend not detected; registering doParallel with",
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
            y <- data_mat[, resp_index]

            # Set the remaining p variables as predictors
            X_mat <- data_mat[, -resp_index]

            E <- stats::rnorm(n, 0, 1) # removing this causes discrepency in discrete case

            # If CS, choose pi and sigmasq according to the Carbonetto-Stephens model
            if (CS){
              idmod <- varbvs::varbvs(X_mat, y, Z = Z[ , 1], verbose = FALSE)
              sigmasq <- mean(idmod$sigma)
              pi_vec <- mean(1 / (1 + exp(-idmod$logodds))) # need to convert to log base 10
            }

            # perform the variational updates and save the results to res
            var_updates(X_mat, Z, D, y, alpha, mu, sigmasq, sigmabetasq_vec,
                        pi_vec, tolerance, max_iter, monitor_final_elbo,
                        monitor_cand_elbo, monitor_period, warnings, resp_index)
            }
          )
      },
      error = function(msg) stop(paste(
        "Parallel execution failed; error message: ", msg))
    )

    # shut down the workers if desired
    if(stop_cluster) doParallel::stopImplicitCluster()

  }else{

    # otherwise, the variational update will be executed sequentially

    # list to store each of the results from var_updates
    res <- vector("list", p + 1)

    for (resp_index in 1:(p + 1)) {

      # Set variable number `resp_index` as the response
      y <- data_mat[, resp_index]

      # Set the remaining p variables as predictors
      X_mat <- data_mat[, -resp_index]

      E <- stats::rnorm(n, 0, 1) # removing this causes discrepency in discrete case

      # If CS, choose pi and sigmasq according to the Carbonetto-Stephens model
      if (CS){
        idmod <- varbvs::varbvs(X_mat, y, Z = Z[ , 1], verbose = FALSE)
        sigmasq <- mean(idmod$sigma)
        pi_vec <- mean(1 / (1 + exp(-idmod$logodds))) # need to convert to log base 10
      }

      # perform the variational updates and save the results to res
      res[[resp_index]] <- var_updates(X_mat, Z, D, y, alpha, mu, sigmasq,
                                       sigmabetasq_vec, pi_vec, tolerance,
                                       max_iter, monitor_final_elbo,
                                       monitor_cand_elbo, monitor_period,
                                       warnings, resp_index)
    }
  }

  # gather the VB_details lists, alpha_matrix matrices, and warnings_vec vectors
  # into their own lists/ vectors
  VB_details <- lapply(res, `[[`, "VB_details")
  alpha_matrices <- lapply(res, `[[`, "alpha_matrix")
  warnings_vec <- unlist(lapply(res, `[[`, "warnings"))

  if (warnings){

    # if there are any warnings to be displayed in warnings_vec, display
    # them
    for (warning in warnings_vec){
      warning(warning)
    }

    # grid warnings - for each of the candidate grids (pi and sigmabeta_sq), if
    # points along the grid boundaries were selected and the grid had more than
    # 2 candidates, display a warning

    # sigmabetasq_vec
    if (length(sigmabetasq_vec) > 2){

      # get the selected values of sigmabeta_sq for each response
      final_sigmabeta_sq <- unlist(lapply(VB_details, `[[`, "sigmabeta_sq"))

      # count the number of final_sigmabeta_sq that were on the boundary of the
      # grid
      grid_boundary <- sigmabetasq_vec[c(1, length(sigmabetasq_vec))]
      on_boundary <- sum(final_sigmabeta_sq %in% grid_boundary)

      # if any of the final sigma were on the boundary, display a warning
      if (on_boundary > 0){
        warning(paste0("For ", on_boundary, "/", p + 1,
                       " responses, the selected value of sigmabeta_sq was on the grid boundary. See return value VB_details"))
      }
    }

    # pi_vec
    if (length(pi_vec) > 2){

      # get the selected values of pi for each response
      final_pi <- unlist(lapply(VB_details, `[[`, "pi"))

      # count the number of final_pi that were on the boundary of the grid
      grid_boundary <- pi_vec[c(1, length(pi_vec))]
      on_boundary <- sum(final_pi %in% grid_boundary)

      # if any of the final pi were on the boundary, display a warning
      if (on_boundary > 0){
        warning(paste0("For ", on_boundary, "/", p + 1,
                       " responses, the selected value of pi was on the grid boundary. See return value VB_details"))
      }
    }
  }

  # Creating the graphs:
  # transform p + 1 n by n matrices to n p + 1 by p + 1 matrices using alpha_matrices
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

  # stop timer and see how much time has elapsed
  if (print_time) print(Sys.time() - start_time)

  return(list(graphs = graphs, inclusion_probs = incl_probs,
              alpha_matrices = incl_probs_asym, unique_graphs = unique_graphs,
              VB_details = VB_details, weights = D, bandwidths = bandwidths))
}
