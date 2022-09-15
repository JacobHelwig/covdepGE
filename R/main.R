## -----------------------------------------------------------------------------
## Distributed under GPL (≥ 3) license
#'
#' @title Covariate Dependent Graph Estimation
#' @aliases covdepGE-method
#' @export
## -----------------------------DESCRIPTION-------------------------------------
#' @description Model the conditional dependence structure of `X` as a function
#' of `Z` as described in (1)
## -----------------------------ARGUMENTS---------------------------------------
#' @param X \eqn{n \times p}{n x p} numeric matrix; data matrix. For best
#' results, \eqn{n} should be greater than \eqn{p}
#'
#' @param Z `NULL` OR \eqn{n \times q}{n x q} numeric matrix; extraneous
#' covariates. If `NULL`, `Z` will be treated as constant for all observations,
#' i.e.:
#'
#' ```
#' Z <- rep(0, nrow(X))
#' ```
#'
#' If `Z` is constant, the estimated graph will be homogeneous throughout the
#' data. `NULL` by default
#'
#' @param hp_method `character` in `c("grid_search","model_average","hybrid")`;
#' method for selecting hyperparameters from the the hyperparameter grid. The
#' grid will be generated as the Cartesian product of `ssq`, `sbsq`, and `pip`.
#' Fix \eqn{X_j}{Xj}, the \eqn{j}-th column of `X`, as the response; then, the
#' hyperparameters will be selected as follows:
#'
#'  \itemize{
#'    \item If `"grid_search"`, the point in the hyperparameter grid that
#'    maximizes the total ELBO summed across all \eqn{n} regressions will be
#'    selected
#'    \item If `"model_average"`, then all posterior quantities will be an
#'    average of the variational estimates resulting from the model fit for each
#'    point in the hyperparameter grid. The unnormalized averaging weights for
#'    each of the \eqn{n} regressions are the exponentiated ELBO
#'    \item If `"hybrid"`, then models will be averaged over `pip` as in
#'    `"model_average"`, with \eqn{\sigma^2}{sigma^2} and
#'    \eqn{\sigma_\beta^2}{sigma_beta^2} chosen for each \eqn{\pi}{pi} in `pip`
#'    by maximizing the total ELBO over the grid defined by the Cartesian
#'    product of `ssq` and `sbsq` as in `"grid_search"`
#' }
#'
#' `"hybrid"` by default
#'
#' @param ssq `NULL` OR numeric vector with positive entries; candidate values
#' of the hyperparameter \eqn{\sigma^2}{sigma^2} (prior residual variance). If
#' `NULL`, `ssq` will be generated for each variable \eqn{X_j}{Xj} fixed as the
#' response as:
#'
#' ```
#' ssq <- seq(ssq_lower, ssq_upper, length.out = nssq)
#' ```
#'
#' `NULL` by default
#'
#' @param sbsq `NULL` OR numeric vector with positive entries; candidate values
#' of the hyperparameter \eqn{\sigma_\beta^2}{sigma_beta^2} (prior slab
#' variance). If `NULL`, `sbsq` will be generated for each variable
#' \eqn{X_j}{Xj} fixed as the response as:
#'
#' ```
#' sbsq <- seq(sbsq_lower, sbsq_upper, length.out = nsbsq)
#' ```
#'
#' `NULL` by default
#'
#' @param pip `NULL` OR numeric vector with entries in \eqn{(0, 1)}; candidate
#' values of the hyperparameter \eqn{\pi}{pi} (prior inclusion probability). If
#' `NULL`, `pip` will be generated for each variable \eqn{X_j}{Xj} fixed as the
#' response as:
#'
#' ```
#' pip <- seq(pip_lower, pi_upper, length.out = npip)
#' ```
#' `NULL` by default
#'
#' @param nssq  positive integer; number of points to generate for `ssq` if
#' `ssq` is `NULL`. `5` by default
#'
#' @param nsbsq positive integer; number of points to generate for `sbsq` if
#' `sbsq` is `NULL`. `5` by default
#'
#' @param npip positive integer; number of points to generate for `pip` if `pip`
#' is `NULL`. `5` by default
#'
#' @param ssq_mult positive numeric; if `ssq` is `NULL`, then for each variable
#' \eqn{X_j}{Xj} fixed as the response:
#'
#' ```
#' ssq_upper <- ssq_mult * stats::var(X_j)
#' ```
#'
#' Then, `ssq_upper` will be the greatest value in `ssq` for variable
#' \eqn{X_j}{Xj}. `1.5` by default
#'
#' @param ssq_lower positive numeric; if `ssq` is `NULL`, then `ssq_lower` will
#' be the least value in `ssq`. `1e-5` by default
#'
#' @param snr_upper positive numeric; upper bound on the signal-to-noise ratio.
#' If `sbsq` is `NULL`, then for each variable \eqn{X_j}{Xj} fixed as the
#' response:
#'
#' ```
#' s2_sum <- sum(apply(X, 2, stats::var))
#' sbsq_upper <- snr_upper / (pip_upper * s2_sum)
#' ```
#'
#' Then, `sbsq_upper` will be the greatest value in `sbsq`. `25` by default
#'
#' @param sbsq_lower positive numeric; if `sbsq` is `NULL`, then `sbsq_lower`
#' will be the least value in `sbsq`. `1e-5` by default
#'
#' @param pip_lower numeric in \eqn{(0, 1)}; if `pip` is `NULL`, then
#' `pip_lower` will be the least value in `pip`. `1e-5` by default
#'
#' @param pip_upper `NULL` OR  numeric in \eqn{(0, 1)}; if `pip` is `NULL`, then
#' `pip_upper` will be the greatest value in `pip`. If `sbsq` is `NULL`,
#' `pip_upper` will be used to calculate `sbsq_upper`. If `NULL`, `pip_upper`
#' will be calculated for each variable \eqn{X_j}{Xj} fixed as the response as:
#'
#' ```
#' lasso <- glmnet::cv.glmnet(X, X_j)
#' non0 <- sum(glmnet::coef.glmnet(lasso, s = "lambda.1se")[-1] != 0)
#' non0 <- min(max(non0, 1), p - 1)
#' pip_upper <- non0 / p
#' ```
#' `NULL` by default
#'
#' @param tau `NULL` OR positive numeric OR numeric vector of length \eqn{n}
#' with positive entries; bandwidth parameter. Greater values allow for more
#' information to be shared between observations. Allows for global or
#' observation-specific specification. If `NULL`, use 2-step KDE methodology as
#' described in (2) to calculate observation-specific bandwidths. `NULL` by
#' default
#'
#' @param norm numeric in \eqn{[1, \infty]}{[1, Inf]}; norm to use when
#' calculating weights. `Inf` results in infinity norm. `2` by default
#'
#' @param center_X logical; if `TRUE`, center `X` column-wise to mean \eqn{0}.
#' `TRUE` by default
#'
#' @param scale_Z logical; if `TRUE`, center and scale `Z` column-wise to mean
#' \eqn{0}, standard deviation \eqn{1} prior to calculating the weights. `TRUE`
#' by default
#'
#' @param alpha_tol positive numeric; end CAVI when the Frobenius norm of the
#' change in the alpha matrix is within `alpha_tol`. `1e-5` by default
#'
#' @param max_iter positive integer; if tolerance criteria has not been met by
#' `max_iter` iterations, end CAVI. `100` by default
#'
#' @param max_iter_grid positive integer; if tolerance criteria has not been
#' met by `max_iter_grid` iterations during grid search, end CAVI. After grid
#' search has completed, CAVI is performed with the final hyperparameters
#' selected by grid search for at most `max_iter` iterations. Does not apply to
#' `hp_method = "model_average"`. `10` by default
#'
#' @param edge_threshold numeric in \eqn{(0, 1)}; a graph for each observation
#' will be constructed by including an edge between variable \eqn{i} and
#' variable \eqn{j} if, and only if, the \eqn{(i, j)} entry of the symmetrized
#' posterior inclusion probability matrix corresponding to the observation is
#' greater than `edge_threshold`. `0.5` by default
#'
#' @param sym_method `character` in `c("mean","max","min")`; to symmetrize
#' the posterior inclusion probability matrix for each observation, the
#' \eqn{(i, j)} and \eqn{(j, i)} entries will be post-processed as `sym_method`
#' applied to the \eqn{(i, j)} and \eqn{(j, i)} entries. `"mean"` by default
#'
#' @param parallel logical; if `TRUE`, hyperparameter selection and CAVI for
#' each of the \eqn{p} variables will be performed in parallel using `foreach`.
#' Parallel backend may be registered prior to making a call to `covdepGE`. If
#' no active parallel backend can be detected, then parallel backend will be
#' automatically registered using:
#'
#' ```
#' doParallel::registerDoParallel(num_workers)
#' ```
#'
#' `FALSE` by default
#'
#' @param num_workers `NULL` OR positive integer less than or equal to
#' `parallel::detectCores()`; argument to `doParallel::registerDoParallel` if
#' `parallel = TRUE` and no parallel backend is detected. If `NULL`, then:
#'
#' ```
#' num_workers <- floor(parallel::detectCores() / 2)
#' ```
#'
#' `NULL` by default
#'
#' @param prog_bar logical; if `TRUE`, then a progress bar will be displayed
#' denoting the number of remaining variables to fix as the response and perform
#' CAVI. If `parallel`, no progress bar will be displayed. `TRUE` by default
## -----------------------------RETURNS-----------------------------------------
#' @return Returns object of class `covdepGE` with the following values:
#'
#'  \item{graphs}{list with the following values:
#'
#'    \itemize{
#'      \item `graphs`: list of \eqn{n} numeric matrices of dimension
#'      \eqn{p \times p}{p x p}; the \eqn{l}-th matrix is the adjacency matrix
#'      for the \eqn{l}-th observation
#'      \item `unique_graphs`: list; the \eqn{l}-th element is a list containing
#'      the \eqn{l}-th unique graph and the indices of the observation(s)
#'      corresponding to this graph
#'      \item `inclusion_probs_sym`: list of \eqn{n} numeric matrices of
#'      dimension \eqn{p \times p}{p x p}; the \eqn{l}-th matrix is the
#'      symmetrized posterior inclusion probability matrix for the \eqn{l}-th
#'      observation
#'      \item `inclusion_probs_asym`: list of \eqn{n} numeric matrices of
#'      dimension \eqn{p \times p}{p x p}; the \eqn{l}-th matrix is the
#'      posterior inclusion probability matrix for the \eqn{l}-th observation
#'      prior to symmetrization
#'    }
#'  }
#'
#'  \item{variational_params}{list with the following values:
#'
#'    \itemize{
#'      \item `alpha`: list of \eqn{p} numeric matrices of dimension
#'      \eqn{n \times (p - 1)}{n x (p - 1)}; the \eqn{(i, j)} entry of the
#'      \eqn{k}-th matrix is the variational approximation to the posterior
#'      inclusion probability of the \eqn{j}-th variable in a weighted
#'      regression with variable \eqn{k} fixed as the response, where the
#'      weights are taken with respect to observation \eqn{i}
#'      \item `mu`: list of \eqn{p} numeric matrices of dimension
#'      \eqn{n \times (p - 1)}{n x (p - 1)}; the \eqn{(i, j)} entry of the
#'      \eqn{k}-th matrix is the variational approximation to the posterior slab
#'      mean for the \eqn{j}-th variable in a weighted regression with variable
#'      \eqn{k} fixed as the response, where the weights are taken with respect
#'      to observation \eqn{i}
#'      \item `ssq_var`: list of \eqn{p} numeric
#'      matrices of dimension \eqn{n \times (p - 1)}{n x (p - 1)}; the
#'      \eqn{(i, j)} entry of the \eqn{k}-th matrix is the variational
#'      approximation to the posterior slab variance for the \eqn{j}-th variable
#'      in a weighted regression with variable \eqn{k} fixed as the response,
#'      where the weights are taken with respect to observation \eqn{i}
#'    }
#'  }
#'
#'  \item{hyperparameters}{list of \eqn{p} lists; the \eqn{j}-th list has the
#'  following values for variable \eqn{j} fixed as the response:
#'
#'    \itemize{
#'      \item `grid`: matrix of candidate hyperparameter values, corresponding
#'      ELBO, and iterations to converge
#'      \item `final`: the final hyperparameters chosen by grid search and the
#'      ELBO and iterations to converge for these hyperparameters
#'    }
#'  }
#'
#'  \item{model_details}{list with the following values:
#'
#'    \itemize{
#'      \item `elapsed`: amount of time to fit the model
#'      \item `n`: number of observations
#'      \item `p`: number of variables
#'      \item `ELBO`: ELBO summed across all observations and variables. If
#'      `hp_method` is `"model_average"` or `"hybrid"`, this ELBO is averaged
#'      across the hyperparameter grid using the model averaging weights for
#'      each variable
#'      \item `num_unique`: number of unique graphs
#'      \item `grid_size`: number of points in the hyperparameter grid
#'      \item `args`: list containing all passed arguments of length \eqn{1}
#'    }
#'  }
#'
#'  \item{weights}{list with the following values:
#'
#'    \itemize{
#'      \item `weights`: \eqn{n\times n}{n x n} numeric matrix. The \eqn{(i, j)}
#'      entry is the similarity weight of the \eqn{i}-th observation with
#'      respect to the \eqn{j}-th observation using the \eqn{j}-th observation's
#'      bandwidth
#'      \item `bandwidths`: numeric vector of length \eqn{n}. The \eqn{i}-th
#'      entry is the bandwidth for the \eqn{i}-th observation
#'    }
#'  }
## -----------------------------EXAMPLES----------------------------------------
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # get the data
#' set.seed(12)
#' data <- generateData()
#' X <- data$X
#' Z <- data$Z
#' interval <- data$interval
#' prec <- data$true_precision
#'
#' # get overall and within interval sample sizes
#' n <- nrow(X)
#' n1 <- sum(interval == 1)
#' n2 <- sum(interval == 2)
#' n3 <- sum(interval == 3)
#'
#' # visualize the distribution of the extraneous covariate
#' ggplot(data.frame(Z = Z, interval = as.factor(interval))) +
#'   geom_histogram(aes(Z, fill = interval), color = "black", bins = n %/% 5)
#'
#' # visualize the true precision matrices in each of the intervals
#'
#' # interval 1
#' matViz(prec[[1]], incl_val = TRUE) +
#'   ggtitle(paste0("True precision matrix, interval 1, observations 1,...,", n1))
#'
#' # interval 2 (varies continuously with Z)
#' cat("\nInterval 2, observations ", n1 + 1, ",...,", n1 + n2, sep = "")
#' int2_mats <- prec[interval == 2]
#' int2_inds <- c(5, n2 %/% 2, n2 - 5)
#' lapply(int2_inds, function(j) matViz(int2_mats[[j]], incl_val = TRUE) +
#'          ggtitle(paste("True precision matrix, interval 2, observation", j + n1)))
#'
#' # interval 3
#' matViz(prec[[length(prec)]], incl_val = TRUE) +
#'   ggtitle(paste0("True precision matrix, interval 3, observations ",
#'                  n1 + n2 + 1, ",...,", n1 + n2 + n3))
#'
#' # fit the model and visualize the estimated graphs
#' (out <- covdepGE(X, Z))
#' plot(out)
#'
#' # visualize the posterior inclusion probabilities for variables (1, 3) and (1, 2)
#' inclusionCurve(out, 1, 2)
#' inclusionCurve(out, 1, 3)
#' }
## -----------------------------REFERENCES--------------------------------------
#' @references
#' (1) Sutanoy Dasgupta, Peng Zhao, Prasenjit Ghosh, Debdeep Pati, and Bani
#' Mallick. An approximate Bayesian approach to covariate-dependent graphical
#' modeling. pages 1–59, 2022.
#'
#' (2) Sutanoy Dasgupta, Debdeep Pati, and Anuj Srivastava. A Two-Step Geometric
#' Framework For Density Modeling. *Statistica Sinica*, 30(4):2155–2177, 2020.
## -----------------------------------------------------------------------------
covdepGE <- function(X, Z = NULL, hp_method = "hybrid", ssq = NULL, sbsq = NULL,
                     pip = NULL, nssq = 5, nsbsq = 5, npip = 5, ssq_mult = 1.5,
                     ssq_lower = 1e-5, snr_upper = 25, sbsq_lower = 1e-5,
                     pip_lower = 1e-5, pip_upper = NULL, tau = NULL, norm = 2,
                     center_X = TRUE, scale_Z = TRUE, alpha_tol = 1e-5,
                     max_iter_grid = 10, max_iter = 100, edge_threshold = 0.5,
                     sym_method = "mean", parallel = FALSE, num_workers = NULL,
                     prog_bar = TRUE){

  start_time <- Sys.time()

  # if Z is NULL, make constant
  if (is.null(Z)){
    Z <- rep(0, nrow(X))
  }

  # ensure that X and Z are matrices
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  # get sample size
  n <- nrow(X)
  p <- ncol(X)

  # verify that X and Z have the same number of observations
  if (nrow(Z) != n){
    stop(paste0("Number of observations in X (", n,
                ") and number of observations in Z (", nrow(Z),
                ") do not match"))
  }

  # if the covariates should be centered and scaled, do so ([ , ] for attributes)
  if (scale_Z){

    # if there is only one unique value in Z, warn and do not scale
    if (any(apply(Z, 2, stats::sd) == 0)){
      warning("Cannot scale constant Z")
      Z <- matrix(scale(Z, scale = FALSE)[ , ], n)
    } else{

      # otherwise, center and scale Z
      Z <- matrix(scale(Z)[ , ], n)
    }
  }

  # if X should be centered, do so
  if (center_X) X <- matrix(scale(X, scale = FALSE)[ , ], n)

  # get weights
  D <- get_weights(Z, norm, tau)
  bandwidths <- D$bandwidths
  D <- D$D

  # create list for the weights and bandwidths
  weights <- list(weights = D, bandwidths = bandwidths)

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
      error = function(msg) FALSE,

      # return false if warning
      warning = function(msg) FALSE)

    # display a message that registered workers have been detected
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
          foreach::foreach(j = 1:p, .packages = "covdepGE"),
          {

            # Set variable number \eqn{j} as the response
            y <- X[, j]

            # Set the remaining p variables as predictors
            X_j <- X[, -j, drop = FALSE]

            # perform CAVI and save results to res
            cavi(X_j, Z, D, y, hp_method, ssq, sbsq, pip, nssq, nsbsq, npip,
                 ssq_mult, ssq_lower, snr_upper, sbsq_lower, pip_lower,
                 pip_upper, alpha_tol, max_iter, max_iter_grid)
            }
          )
      },

      # if parallel execution did not finish successfully, display an error
      error = function(msg) stop(paste(
        "Parallel execution failed; error message:", msg))
    )

    # shut down the cluster
    doParallel::stopImplicitCluster()

  }else{

    # otherwise, CAVI will be executed sequentially

    # instantiate the progress bar
    if (prog_bar) pb <- utils::txtProgressBar(0, p, style = 3)

    # list to store each of the results from cavi_search
    res <- vector("list", p)

    for (j in 1:p) {

      # Set variable number \eqn{j} as the response
      y <- X[, j]

      # Set the remaining p variables as predictors
      X_j <- X[, -j, drop = FALSE]

      # perform CAVI and save results to res
      res[[j]] <- cavi(X_j, Z, D, y, hp_method, ssq, sbsq, pip, nssq, nsbsq,
                       npip, ssq_mult, ssq_lower, snr_upper, sbsq_lower,
                       pip_lower, pip_upper, alpha_tol, max_iter, max_iter_grid)

      # update the progress bar
      if (prog_bar) utils::setTxtProgressBar(pb, j)
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

  # create list for the variational parameters
  var_mats <- list(alpha = alpha_matrices, mu = mu_matrices,
                   ssq_var = ssqv_matrices)

  # calculate the total ELBO
  total_elbo <- sum(elbo)

  # get the grid size
  grid_sz <- nrow(hp[[1]]$grid)

  # Graph post-processing

  # transform p n by n matrices to n p by p matrices using alpha_matrices
  # the j, k entry in the l-th matrix is the probability of inclusion of an edge
  # between the j, k variables for the l-th observation
  incl_probs <- replicate(n, matrix(0, p, p), simplify = FALSE)

  # iterate over the p matrices
  for (j in 1:p){

    # fix the j-th alpha matrix
    alpha_mat_j <- alpha_matrices[[j]]

    # iterate over the rows of alpha_mat_j
    for (l in 1:n){

      # the j-th row of the l-th observation's graph is the l-th row of
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

  # create a list for graphs and inclusion probabilties
  graphs <- list(graphs = graphs, unique_graphs = unique_sum,
                 inclusion_probs_sym = incl_probs,
                 inclusion_probs_asym = incl_probs_asym)

  # create a list to return the scalar function arguments
  args <- list(hp_method = hp_method, nssq = nssq, nsbsq = nsbsq, npip = npip,
               ssq_mult = ssq_mult, ssq_lower = ssq_lower,
               snr_upper = snr_upper, sbsq_lower = sbsq_lower,
               pip_lower = pip_lower, pip_upper = pip_upper, norm = norm,
               center_X = center_X, scale_Z = scale_Z, alpha_tol = alpha_tol,
               max_iter = max_iter, max_iter_grid = max_iter_grid,
               edge_threshold = edge_threshold, sym_method = sym_method,
               parallel = parallel, num_workers = num_workers,
               prog_bar = prog_bar)

  # create a list for model details
  model_details <- list(elapsed = NA, n = n, p = p, ELBO = total_elbo,
                        num_unique = length(unique_graphs), grid_size = grid_sz,
                        args = args)

  # record the elapsed time and add it to the model details
  model_details[["elapsed"]] <- Sys.time() - start_time

  # define the list of return values
  ret <- list(graphs = graphs, variational_params = var_mats,
              hyperparameters = hp, model_details = model_details,
              weights = weights)

  # define the class of the return values
  class(ret) <- c("covdepGE", "list")

  return(ret)
}

## -----------------------------------------------------------------------------
#' Distributed under GPL (≥ 3) license
#'
#' @title Print the Return of `covdepGE`
#' @export
#' @rdname covdepGE
## -----------------------------DESCRIPTION-------------------------------------
#' S3 method for printing an object of class `covdepGE`
## -----------------------------ARGUMENTS---------------------------------------
#' @param x object of `class` `covdepGE`; the return of the `covdepGE` function
#'
#' @param ... additional arguments will be ignored
## -----------------------------------------------------------------------------
print.covdepGE <- function(x, ...){
  summary(x)
}

## -----------------------------------------------------------------------------
#' @title Summarize the Return of `covdepGE`
#' @export
#' @rdname covdepGE
## -----------------------------DESCRIPTION-------------------------------------
#' S3 method for producing a summary of an object of class `covdepGE`
## -----------------------------ARGUMENTS---------------------------------------
#' @param object object of `class` `covdepGE`; the return of the `covdepGE`
#' function
#'
#' @param ... additional arguments will be ignored
## -----------------------------------------------------------------------------
summary.covdepGE <- function(object, ...){

  cat("                      Covariate Dependent Graphical Model\n\n")
  spc <- 80
  with(object$model_details,
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
