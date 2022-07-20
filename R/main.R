## -----------------------------------------------------------------------------
#' @title covdepGE: Covariate Dependent Graph Estimation
#' @aliases covdepGE-method
#' @export
## -----------------------------DESCRIPTION-------------------------------------
#' @description Model the conditional dependence structure of data as a function
#' of extraneous covariates as described in (1).
## -----------------------------ARGUMENTS---------------------------------------
#' @param data `n x p numeric matrix`; data
#'
#' @param Z `n x q numeric matrix`; extraneous covariates
#'
#' @param hp_method `character` in
#' `c("grid_search", "model_average", "hybrid")`; method for setting
#' hyperparameter values based on the hyperparameter grid. The grid will be
#' generated as the Cartesian product of `ssq`, `sbsq`, and `pip`.
#'
#' If `"grid_search"`, then the point in the hyperparameter grid that maximizes
#' the ELBO across all individuals for a given variable `y` fixed as the
#' response will be selected.
#'
#' If `"model_average`, then all posterior quantities for each variable `y`
#' fixed as the response will be a convex combination of the models resulting
#' from each point in the hyperparameter grid, where the weightings are both
#' individual and variable specific. Unnormalized weights are calculated using
#' the exponentiated ELBO
#'
#' If `"hybrid"`, then `pip` will be averaged over as with `"model_average"`,
#' while a single point in the grid defined by the Cartesian product of `ssq`
#' and `sbsq` will be selected for each variable `y` fixed as the response and
#' point in `pip`. `"hybrid"` by default
#'
#' @param ssq `NULL` OR `numeric vector` with positive entries; candidate values
#' of the hyperparameter `sigma^2` (prior residual variance). If `NULL`, `ssq`
#' will be generated for each variable `y` fixed as the response as:
#'
#' `ssq <- seq(ssq_lower, ssq_upper, length.out = nssq)`
#'
#' `NULL` by default
#'
#' @param sbsq `NULL` OR `numeric vector` with positive entries; candidate values
#' of the hyperparameter `sigma^2_beta` (prior slab variance). If `NULL`, `sbsq`
#' will be generated for each variable `y` fixed as the response as:
#'
#' `sbsq <- seq(sbsq_lower, sbsq_upper, length.out = nsbsq)`
#'
#' `NULL` by default
#'
#' @param pip `NULL` OR `numeric vector` with entries in `(0, 1)`; candidate
#' values of the hyperparameter `pi` (prior inclusion probability). If `NULL`,
#' `pip` will be generated for each variable `y` fixed as the response as:
#'
#' `pip <- seq(pip_lower, pi_upper, length.out = npip)`
#'
#' `pi_upper`
#'
#' `NULL` by default
#'
#' @param nssq  positive integer; number of points in `ssq` if `ssq` is `NULL`.
#' `5` by default
#'
#' @param nsbsq positive integer; number of points in `sbsq` if `sbsq` is
#' `NULL`. `5` by default
#'
#' @param npip positive integer; number of points in `pip` if `pip` is `NULL`.
#' `5` by default
#'
#' @param ssq_mult positive `numeric`; if `ssq` is `NULL`, then for each variable
#' `y` fixed as the response:
#'
#' `ssq_upper <- ssq_upper_mult * var(y)`
#'
#' Then, `ssq_upper` will be the greatest value in `ssq` for variable `y`. `1.5`
#' by default
#'
#' @param ssq_lower positive `numeric`; if `ssq` is `NULL`, then `ssq_lower` will
#' be the least value in `ssq`. `1e-5` by default
#'
#' @param snr_upper positive `numeric`; if `sbsq` is `NULL`, then for each
#' variable `y` fixed as the response:
#'
#' `s2_sum <- sum(apply(X, 2, var))`
#'
#' `sbsq_upper <- snr_upper / (pip_upper * s2_sum)`
#'
#' Then, `sbsq_upper` will be the greatest value in `sbsq` for variable `y`.
#' `25` by default
#'
#' @param sbsq_lower positive `numeric`; if `sbsq` is `NULL`, then `sbsq_lower`
#' will be the least value in `sbsq`. `1e-5` by default
#'
#' @param pip_lower `numeric` in `(0, 1)`; if `pip` is `NULL`, then
#' `pip_lower` will be the least value in `pip`. `1e-5` by default
#'
#' @param pip_upper `NULL` OR  `numeric` in`(0, 1)`; if `pip` is `NULL`, then
#' `pip_upper` will be the greatest value in `pip`. If `sbsq` is `NULL`,
#' `pip_upper` will be used to find the greatest value in `sbsq`. If `NULL`,
#' `pip_upper` will be generated for each variable `y` fixed as the response as:
#'
#' `lasso <- glmnet::cv.glmnet(X, y)`
#'
#' `non0 <- sum(coef(lasso, s = "lambda.1se")[-1] != 0)`
#'
#' `non0 <- min(max(non0, 1), p - 1)`
#'
#' `pip_upper <- non0 / p`
#'
#' `NULL` by default
#'
#' @param tau `NULL` OR positive `numeric` OR `numeric vector` of length `n`
#' with positive entries; bandwidth parameter. Greater values allow for more
#' information to be shared between individuals. Allows for global or
#' individual-specific specification. If `NULL`, use 2-step KDE methodology as
#' described in (2) to calculate individual-specific bandwidths. `NULL` by
#' default
#'
#' @param norm `numeric` in `[1, Inf]`; norm to use when calculating weights.
#' `Inf` results in infinity norm. `2` by default
#'
#' @param center_data `logical`; if `T`, center `data` column-wise to mean `0`.
#' `T` by default
#'
#' @param scale_Z `logical`; if `T`, center and scale `Z` column-wise to mean `0`,
#' standard deviation 1 prior to calculating the weights. `T` by default
#'
#' @param alpha_tol positive `numeric`; end CAVI when the Frobenius norm of the
#' change in the alpha `matrix` is within `alpha_tol`. `1e-5` by default
#'
#' @param max_iter positive integer; if a tolerance criteria has not been met
#' by `max_iter_grid` iterations, end CAVI. `100` by default
#'
#' @param edge_threshold `numeric` in `(0, 1)`; a graph for each individual
#' will be constructed by including an edge between variable `i` and
#' variable `j` if, and only if, the `(i,j)` entry of the symmetrized
#' posterior inclusion probability `matrix` corresponding to the individual is
#' greater than `edge_threshold`. `0.5` by default
#'
#' @param sym_method `character` in \{`"mean"`, `"max"`, `"min"`\`; to symmetrize
#' the posterior inclusion probability `matrix` for each individual, the
#' `(i,j)` and `(j,i)` entries will be post-processed as
#' `sym_method``((i,j entry), (j,i entry))`. `"mean"` by default
#'
#' @param parallel `logical`; if `T`, hyperparameter selection and CAVI for each
#' of the `p` variables will be performed in parallel using `foreach`.
#' Parallel backend may be registered prior to making a call to `covdepGE`. If
#' no active parallel backend can be detected, then parallel backend will be
#' automatically registered using `doParallel::registerDoParallel(num_workers)`
#'
#' @param num_workers `NULL` OR positive integer less than or equal to
#' `parallel::detectCores()`; argument to `doParallel::registerDoParallel` if
#' `parallel = T` and no parallel backend is detected. If `NULL`, then:
#'
#' `num_workers <- floor(parallel::detectCores() / 2)`
#'
#' `NULL` by default
#'
#' @param prog_bar `logical`; if `T`, then a progress bar will be displayed
#' denoting the number of remaining variables to fix as the response and perform
#' CAVI. If `parallel`, no progress bar will be displayed. `T` by default
#'
## -----------------------------RETURNS-----------------------------------------
#' @return Returns `list` with the following values:
#'=
#'
#' 1. `graphs`: `list` with the following values:
#'
#' - `graphs`: `list` of `n p x p numeric` matrices; the
#' `l`-th `matrix` is the adjacency `matrix` for the `l`-th individual
#'
#' - `unique_graphs`: `list`; the `l`-th element is a `list` containing the
#' `l`-th unique graph and the individual(s) corresponding to this graph
#'
#' - `inclusion_probs_sym`: `list` of `n p x p numeric`
#' matrices; the `l`-th `matrix` is the symmetrized posterior inclusion
#' probability `matrix` for the `l`-th individual
#'
#' - `inclusion_probs_asym`: `list` of `n p x p numeric`
#' matrices; the `l`-th `matrix` is the posterior inclusion probability `matrix`
#' for the `l`-th individual prior to symmetrization
#'
#'
#' 2. `variational_params`: `list` with the following values:
#'
#' - `alpha`: `list` of `p n x (p - 1) numeric` matrices; the
#' `(i, j)` entry of the `k`-th `matrix` is the variational approximation
#' to the posterior inclusion probability of the `j`-th variable in a
#' weighted regression with variable `k` fixed as the response, where the
#' weights are taken with respect to individual `i`
#'
#' - `mu`: `list` of `p n x (p - 1) numeric` matrices; the
#' `(i, j)` entry of the `k`-th `matrix` is the variational approximation
#' to the posterior slab mean for the `j`-th variable in a weighted
#' regression with variable `k` fixed as the response, where the weights are
#' taken with respect to individual `i`
#'
#' - `ssq_var`: `list` of `p n x (p - 1) numeric` matrices; the
#' `(i, j)` entry of the `k`-th `matrix` is the variational approximation
#' to the posterior slab variance for the `j`-th variable in a weighted
#' regression with variable `k` fixed as the response, where the weights are
#' taken with respect to individual `i`
#'
#'
#' 3. `hyperparameters`: `list` of `p` lists; the `j`-th `list` has the
#' following values for variable `j` fixed as the response:
#'
#' - `grid`: `matrix` of candidate hyperparameter values, corresponding ELBO,
#' and iterations to converge
#'
#' - `final`: the final hyperparameters chosen by grid search
#'
#' 4. `model_details`: `list` with the following values:
#'
#' - `elapsed`: amount of time to fit the model
#'
#' - `n`: number of individuals
#'
#' - `p`: number of variables
#'
#' - `ELBO`: ELBO summed across all individuals and variables. If
#' `hp_method` is `"model_average"` or `"hybrid"`, this ELBO is averaged across
#' the hyperparameter grid using the model averaging weights
#'
#' - `num_unique`: number of unique graphs
#'
#' - `grid_size`: number of points in the hyperparameter grid
#'
#' - `args`: `list` containing all passed arguments of `length` 1
#'
#'
#' 5. `weights`: `list` with the following values:
#'
#' - `weights`: `n x n numeric matrix`. The `(i, j)` entry is the weight of the
#' `i`-th individual with respect to the `j`-th individual using the `j`-th
#' individual's bandwidth
#'
#' - `bandwidths`: `numeric vector` of length `n`. The `i`-th entry is the
#' bandwidth for the `i`-th individual
## -----------------------------EXAMPLES----------------------------------------
#' @examples
## -----------------------------REFERENCES--------------------------------------
#' @references
#' 1. Dasgupta S., Ghosh P., Pati D., Mallick B., *An approximate Bayesian
#' approach to covariate dependent graphical modeling*, 2021
#'
#' 2. Dasgupta S., Pati D., Srivastava A., *A Two-Step Geometric Framework For
#' Density Modeling*, Statistica Sinica, 2020
## -----------------------------------------------------------------------------
covdepGE <- function(data, Z, hp_method = "hybrid", ssq = NULL, sbsq = NULL,
                     pip = NULL, nssq = 3, nsbsq = 3, npip = 3, ssq_mult = 1.5,
                     ssq_lower = 1e-5, snr_upper = 25, sbsq_lower = 1e-5,
                     pip_lower = 1e-5, pip_upper = NULL, tau = NULL, norm = 2,
                     center_data = T, scale_Z = T, alpha_tol = 1e-5,
                     max_iter = 100, edge_threshold = 0.5, sym_method = "mean",
                     parallel = F, num_workers = NULL, prog_bar = T){

  start_time <- Sys.time()

  # ensure that data_mat and Z are matrices
  data <- as.matrix(data)
  Z <- as.matrix(Z)

  # get sample size and number of parameters
  n <- nrow(data)
  p <- ncol(data)

  # verify that data and Z have the same number of observations
  if (nrow(Z) != n){
    stop(paste0("Number of observations in data (", n,
                ") and number of observations in Z (", nrow(Z),
                ") do not match"))
  }

  # if the covariates should be centered and scaled, do so ([ , ] for attributes)
  if (scale_Z) Z <- matrix(scale(Z)[ , ], n)

  # if the data should be centered, do so
  if (center_data) data <- matrix(scale(data, T, F)[ , ], n)

  # get weights
  D <- get_weights(Z, norm, tau)
  bandwidths <- D$bandwidths
  D <- D$D

  # `list` for the weights and bandwidths
  weights = list(weights = D, bandwidths = bandwidths)

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
      message(paste("Detected", foreach::getDoParWorkers(), "workers"))
    }else{

      # otherwise, register parallel backend
      warning(paste(
        "No registered workers detected; registering doParallel with",
        num_workers, "workers"))

      # if num_workers has not been provided, get the number of workers
      if (is.null(num_workers)){
        num_workers <- floor(parallel::detectCores() / 2)
      }

      # perform registration
      doParallel::registerDoParallel(cores = num_workers)
    }

    # attempt to execute the variational update in parallel
    res <- tryCatch(
      {
        foreach::`%dopar%`(
          foreach::foreach(resp_index = 1:(p), .packages = "covdepGE"),
          {

            # Set variable number `resp_index` as the response
            y <- data[, resp_index]

            # Set the remaining p variables as predictors
            X <- data[, -resp_index, drop = F]

            # perform the grid search and final CAVI; save the results to res
            cavi(X, Z, D, y, hp_method, ssq, sbsq, pip, nssq, nsbsq, npip,
                 ssq_mult, ssq_lower, snr_upper, sbsq_lower, pip_lower,
                 pip_upper, alpha_tol, max_iter)
            }
          )
      },

      # if parallel execution did not finish successfully, display an error
      error = function(msg) stop(paste(
        "Parallel execution failed; error message: ", msg))
    )

    # shut down the cluster
    doParallel::stopImplicitCluster()

  }else{

    # otherwise, CAVI will be executed sequentially

    # instantiate the progress bar
    if (prog_bar) pb <- utils::txtProgressBar(0, p, style = 3)

    # `list` to store each of the results from cavi_search
    res <- vector("list", p)

    for (resp_index in 1:p) {

      # Set variable number `resp_index` as the response
      y <- data[, resp_index]

      # Set the remaining p variables as predictors
      X <- data[, -resp_index, drop = F]

      # perform the grid search and final CAVI; save the results to res
      res[[resp_index]] <- cavi(X, Z, D, y, hp_method, ssq, sbsq, pip, nssq,
                                nsbsq, npip, ssq_mult, ssq_lower, snr_upper,
                                sbsq_lower, pip_lower, pip_upper, alpha_tol,
                                max_iter)

      # update the progress bar
      if (prog_bar) utils::setTxtProgressBar(pb, resp_index)
    }

    # close the progress bar
    if (prog_bar) close(pb)
  }

  # gather the variational parameter matrices, hyperparameter details, and
  # final elbo for each variable into lists/ vectors
  alpha_matrices <- lapply(res, `[[`, "alpha_matrix")
  mu_matrices <- lapply(res, `[[`, "mu_matrix")
  ssqv_matrices <- lapply(res, `[[`, "ssqv_matrix")
  hp <- lapply(res, `[[`, "hyperparameters")
  elbo <- sapply(res, `[[`, "elbo")

  # name elements of these lists by variable
  names(hp) <- names(alpha_matrices) <- names(mu_matrices) <- names(
    ssqv_matrices) <- paste0("variable", 1:(p))

  # `list` for the variational parameters
  var_mats <- list(alpha = alpha_matrices, mu = mu_matrices,
                   ssq_var = ssqv_matrices)

  # calculate the total ELBO
  total_elbo <- sum(elbo)

  # get the grid size; if hp_method is "hybrid", include the number of pip
  grid_sz <- nrow(hp[[1]]$grid) * ifelse(hp_method == "hybrid",
                                         nrow(hp[[1]]$final), 1)

  # Graph post-processing

  # transform p n by n matrices to n p by p matrices using alpha_matrices
  # the j, k entry in the l-th matrix is the probability of inclusion of an edge
  # between the j, k variables for the l-th individual
  incl_probs <- replicate(n, matrix(0, p, p), simplify = F)

  # iterate over the p matrices
  for (j in 1:p){

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

  }else{

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

  # create a nested `list` where the j-th inner `list` has three values; the j-th
  # unique graph, the individuals corresponding to that graph, and a summary of
  # the individuals corresponding to that graph
  unique_graphs <- lapply(1:length(unique_graphs), function(gr_idx)
    list(graph = unique_graphs[[gr_idx]], individuals = indv_graphs[[gr_idx]],
         individuals_summary = indv_graphs_sum[[gr_idx]]))
  names(unique_graphs) <- paste0("graph", 1:length(unique_graphs))

  # create a `list` for graphs and inclusion probabilties
  graphs = list(graphs = graphs, unique_graphs = unique_graphs,
                inclusion_probs_sym = incl_probs,
                inclusion_probs_asym = incl_probs_asym)

  # create a `list` to return the scalar function arguments
  args <- list(hp_method = hp_method, nssq = nssq, nsbsq = nsbsq, npip = npip,
               ssq_mult = ssq_mult, ssq_lower = ssq_lower,
               snr_upper = snr_upper, sbsq_lower = sbsq_lower,
               pip_lower = pip_lower, pip_upper = pip_upper, norm = norm,
               center_data = center_data, scale_Z = scale_Z,
               alpha_tol = alpha_tol, max_iter = max_iter,
               edge_threshold = edge_threshold, sym_method = sym_method,
               parallel = parallel, num_workers = num_workers,
               prog_bar = prog_bar)

  # create a `list` for model details
  model_details <- list(elapsed = NA, n = n, p = p, ELBO = total_elbo,
                        num_unique = length(unique_graphs), grid_size = grid_sz,
                        args = args)

  # record the elapsed time and add it to the model details
  model_details[["elapsed"]] <- Sys.time() - start_time

  # define the `list` of return values
  ret <- list(graphs = graphs, variational_params = var_mats,
              hyperparameters = hp, model_details = model_details,
              weights = weights)

  # define the class of the return values
  class(ret) <- c("covdepGE", "list")

  return(ret)
}

## -----------------------------------------------------------------------------
#' @title print.covdepGE
#' @export
#' @rdname covdepGE
## -----------------------------DESCRIPTION-------------------------------------
#' S3 method for printing an object of class `covdepGE`
## -----------------------------ARGUMENTS---------------------------------------
#' @param x object of class `covdepGE`; the return of the `covdepGE` function
#'
#' @param ... additional arguments will be ignored
## -----------------------------------------------------------------------------
print.covdepGE <- function(x, ...){

  cat("                      Covariate Dependent Graphical Model\n\n")

  spc <- 80

  with(x$model_details,
       {
         # print ELBO and number of unique graphs
         elbo_str <- paste0("ELBO: ", round(ELBO, 2))
         cat(sprintf(paste0("%-s%", spc - nchar(elbo_str), "s"), elbo_str,
                     paste0("# Unique Graphs: ", num_unique, "\n")))

         # print data dimensions and the grid size
         data_dim <- paste0("n: ", n, ", variables: ", p)
         cat(sprintf(paste0("%-s%", spc - nchar(data_dim), "s"), data_dim,
                     paste0("Hyperparameter grid size: ", grid_size,
                            " points\n")))

         # print time to fit
         time_units <- attr(elapsed, "units")
         duration <- as.numeric(elapsed)
         cat("Model fit completed in", round(duration, 3), time_units, "\n")
       }
  )
}
