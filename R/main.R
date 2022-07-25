## -----------------------------------------------------------------------------
#' @title Covariate Dependent Graph Estimation
#' @aliases covdepGE-method
#' @export
## -----------------------------DESCRIPTION-------------------------------------
#' @description Model the conditional dependence structure of `X` as a function
#' of `Z` as described in (1).
## -----------------------------ARGUMENTS---------------------------------------
#' @param X `n` `x` `p` `numeric` `matrix`; data `matrix`
#'
#' @param Z `n` `x` `q` `numeric` `matrix`; extraneous covariates
#'
#' @param hp_method `character` in
#' `c("grid_search", "model_average", "hybrid")`; method for selecting
#' hyperparameters from the the hyperparameter grid. The grid will be generated
#' as the Cartesian product of `ssq`, `sbsq`, and `pip`. Fix `X_j`, the `j-th`
#' column of `X`, as the response; then, the hyperparameters will be selected as
#' follows:
#'
#'  \itemize{
#'    \item If `"grid_search"`, the point in the hyperparameter grid that
#'    maximizes the total ELBO across summed across all `n` regressions will be
#'    selected
#'    \item If `"model_average`, then all posterior quantities will be a convex
#'    combination of the variational estimates resulting from the model fit for
#'    each point in the hyperparameter grid. Note that the weightings are
#'    observation-specific, as unnormalized weights are calculated using the
#'    exponentiated ELBO
#'    \item If `"hybrid"`, then `pip` will be averaged over as with
#'    `"model_average"`, while a single point in the grid defined by the
#'    Cartesian product of `ssq` and `sbsq` will be selected via grid search for
#'    each point in `pip`
#' }
#'
#' `"hybrid"` by default
#'
#' @param ssq `NULL` OR `numeric vector` with positive entries; candidate values
#' of the hyperparameter `sigma^2` (prior residual variance). If `NULL`, `ssq`
#' will be generated for each variable `X_j` fixed as the response as:
#'
#' ```
#' ssq <- seq(ssq_lower, ssq_upper, length.out = nssq)
#' ```
#'
#' `NULL` by default
#'
#' @param sbsq `NULL` OR `numeric vector` with positive entries; candidate
#' values of the hyperparameter `sigma^2_beta` (prior slab variance). If `NULL`,
#' `sbsq` will be generated for each variable `X_j` fixed as the response as:
#'
#' ```
#' sbsq <- seq(sbsq_lower, sbsq_upper, length.out = nsbsq)
#' ```
#'
#' `NULL` by default
#'
#' @param pip `NULL` OR `numeric vector` with entries in `(0,1)`; candidate
#' values of the hyperparameter `pi` (prior inclusion probability). If `NULL`,
#' `pip` will be generated for each variable `X_j` fixed as the response as:
#'
#' ```
#' pip <- seq(pip_lower, pi_upper, length.out = npip)
#' ```
#' `NULL` by default
#'
#' @param nssq  positive `integer`; number of points to generate for `ssq` if
#' `ssq` is `NULL`. `5` by default
#'
#' @param nsbsq positive `integer`; number of points to generate for `sbsq` if
#' `sbsq` is `NULL`. `5` by default
#'
#' @param npip positive `integer`; number of points to generate for `pip` if
#' `pip` is `NULL`. `5` by default
#'
#' @param ssq_mult positive `numeric`; if `ssq` is `NULL`, then for each variable
#' `X_j` fixed as the response:
#'
#' ```
#' ssq_upper <- ssq_mult * stats::var(y)
#' ```
#'
#' Then, `ssq_upper` will be the greatest value in `ssq` for variable `X_j`.
#' `1.5` by default
#'
#' @param ssq_lower positive `numeric`; if `ssq` is `NULL`, then `ssq_lower`
#' will be the least value in `ssq`. `1e-5` by default
#'
#' @param snr_upper positive `numeric`; upper bound on the signal to noise
#' ratio. If `sbsq` is `NULL`, then for each variable `X_j` fixed as the
#' response:
#'
#' ```
#' s2_sum <- sum(apply(X, 2, stats::var))
#' sbsq_upper <- snr_upper / (pip_upper * s2_sum)
#' ```
#'
#' Then, `sbsq_upper` will be the greatest value in `sbsq` for variable `X_j`.
#' `25` by default
#'
#' @param sbsq_lower positive `numeric`; if `sbsq` is `NULL`, then `sbsq_lower`
#' will be the least value in `sbsq`. `1e-5` by default
#'
#' @param pip_lower `numeric` in `(0,1)`; if `pip` is `NULL`, then
#' `pip_lower` will be the least value in `pip`. `1e-5` by default
#'
#' @param pip_upper `NULL` OR  `numeric` in`(0,1)`; if `pip` is `NULL`, then
#' `pip_upper` will be the greatest value in `pip`. If `sbsq` is `NULL`,
#' `pip_upper` will be used to calculate `sbsq_upper`. If `NULL`, `pip_upper`
#' will be calculated for each variable `X_j` fixed as the response as:
#'
#' ```
#' lasso <- glmnet::cv.glmnet(X, y)
#' non0 <- sum(glmnet::coef.glmnet(lasso, s = "lambda.1se")[-1] != 0)
#' non0 <- min(max(non0, 1), p - 1)
#' pip_upper <- non0 / p
#' ```
#' `NULL` by default
#'
#' @param tau `NULL` OR positive `numeric` OR `numeric vector` of length `n`
#' with positive entries; bandwidth parameter. Greater values allow for more
#' information to be shared between observations. Allows for global or
#' observation-specific specification. If `NULL`, use 2-step KDE methodology as
#' described in (2) to calculate observation-specific bandwidths. `NULL` by
#' default
#'
#' @param norm `numeric` in `[1,Inf]`; norm to use when calculating weights.
#' `Inf` results in infinity norm. `2` by default
#'
#' @param center_X `logical`; if `T`, center `X` column-wise to mean `0`.
#' `T` by default
#'
#' @param scale_Z `logical`; if `T`, center and scale `Z` column-wise to mean
#' `0`, standard deviation `1` prior to calculating the weights. `T` by default
#'
#' @param alpha_tol positive `numeric`; end CAVI when the Frobenius norm of the
#' change in the alpha `matrix` is within `alpha_tol`. `1e-5` by default
#'
#' @param max_iter positive `integer`; if tolerance criteria has not been met by
#' `max_iter` iterations, end CAVI. `100` by default
#'
#' @param max_iter_grid positive `integer`; if tolerance criteria has not been
#' met by `max_iter_grid` iterations during grid search, end CAVI. After grid
#' search has completed, CAVI is performed with the final hyperparameters
#' selected by grid search for at most `max_iter` iterations. Does not apply to
#' `hp_method = "model_average"`. `10` by default
#'
#' @param edge_threshold `numeric` in `(0,1)`; a graph for each observation
#' will be constructed by including an edge between variable `i` and
#' variable `j` if, and only if, the `(i,j)` entry of the symmetrized
#' posterior inclusion probability `matrix` corresponding to the observation is
#' greater than `edge_threshold`. `0.5` by default
#'
#' @param sym_method `character` in `c("mean"`, `"max"`, `"min")`; to symmetrize
#' the posterior inclusion probability `matrix` for each observation, the
#' `(i,j)` and `(j,i)` entries will be post-processed as `sym_method` applied to
#' the `(i,j)` and `(j,i)` entries. `"mean"` by default
#'
#' @param parallel `logical`; if `T`, hyperparameter selection and CAVI for each
#' of the `p` variables will be performed in parallel using `foreach`.
#' Parallel backend may be registered prior to making a call to `covdepGE`. If
#' no active parallel backend can be detected, then parallel backend will be
#' automatically registered using:
#'
#' ```
#' doParallel::registerDoParallel(num_workers)
#' ```
#'
#' @param num_workers `NULL` OR positive `integer` less than or equal to
#' `parallel::detectCores()`; argument to `doParallel::registerDoParallel` if
#' `parallel = T` and no parallel backend is detected. If `NULL`, then:
#'
#' ```
#' num_workers <- floor(parallel::detectCores() / 2)
#' ```
#'
#' `NULL` by default
#'
#' @param prog_bar `logical`; if `T`, then a progress bar will be displayed
#' denoting the number of remaining variables to fix as the response and perform
#' CAVI. If `parallel`, no progress bar will be displayed. `T` by default
#'
## -----------------------------RETURNS-----------------------------------------
#' @return Returns `list` with the following values:
#'
#' \enumerate{
#'
#'  \item `graphs`: `list` with the following values:
#'
#'    \itemize{
#'      \item `graphs`: `list` of `n` `p` `x` `p` `numeric` matrices; the `l`-th
#'      `matrix` is the adjacency `matrix` for the `l`-th observation
#'      \item `unique_graphs`: `list`; the `l`-th element is a `list` containing
#'      the `l`-th unique graph and the indices of the observation(s)
#'      corresponding to this graph
#'      \item `inclusion_probs_sym`: `list` of `n` `p` `x` `p` `numeric`
#'      matrices; the `l`-th `matrix` is the symmetrized posterior inclusion
#'      probability `matrix` for the `l`-th observation
#'      \item `inclusion_probs_asym`: `list` of `n` `p` `x` `p` `numeric`
#'      matrices; the `l`-th `matrix` is the posterior inclusion probability
#'      `matrix` for the `l`-th observation prior to symmetrization
#'    }
#'
#'  \item `variational_params`: `list` with the following values:
#'
#'    \itemize{
#'      \item `alpha`: `list` of `p` `n` `x` `(p-1)` `numeric` matrices; the
#'      `(i,j)` entry of the `k`-th `matrix` is the variational approximation to
#'      the posterior inclusion probability of the `j`-th variable in a weighted
#'      regression with variable `k` fixed as the response, where the weights
#'      are taken with respect to observation `i`
#'      \item `mu`: `list` of `p` `n` `x` `(p-1)` `numeric` matrices; the
#'      `(i,j)` entry of the `k`-th `matrix` is the variational approximation to
#'      the posterior slab mean for the `j`-th variable in a weighted regression
#'      with variable `k` fixed as the response, where the weights are taken
#'      with respect to observation `i`
#'      \item `ssq_var`: `list` of `p` `n` `x` `(p-1)` `numeric` matrices; the
#'      `(i,j)` entry of the `k`-th `matrix` is the variational approximation
#'      to the posterior slab variance for the `j`-th variable in a weighted
#'      regression with variable `k` fixed as the response, where the weights
#'      are taken with respect to observation `i`
#'    }
#'
#'  \item `hyperparameters`: `list` of `p` lists; the `j`-th `list` has the
#'  following values for variable `j` fixed as the response:
#'
#'    \itemize{
#'      \item `grid`: `matrix` of candidate hyperparameter values, corresponding
#'      ELBO, and iterations to converge
#'      \item `final`: the final hyperparameters chosen by grid search, and the
#'      ELBO and iterations to converge for these hyperparameters
#'    }
#'
#'  \item `model_details`: `list` with the following values:
#'
#'    \itemize{
#'      \item `elapsed`: amount of time to fit the model
#'      \item `n`: number of observations
#'      \item `p`: number of variables
#'      \item `ELBO`: ELBO summed across all observations and variables. If
#'      `hp_method` is `"model_average"` or `"hybrid"`, this ELBO is averaged
#'      across the hyperparameter grid using the model averaging weights
#'      \item `num_unique`: number of unique graphs
#'      \item `grid_size`: number of points in the hyperparameter grid
#'      \item `args`: `list` containing all passed arguments of `length 1`
#'    }
#'
#'  \item `weights`: `list` with the following values:
#'
#'  \itemize{
#'    \item `weights`: `n` `x` `n` `numeric` `matrix`. The `(i,j)` entry is the
#'    weight of the `i`-th observation with respect to the `j`-th observation
#'    using the `j`-th observation's bandwidth
#'    \item `bandwidths`: `numeric vector` of length `n`. The `i`-th entry is
#'    the bandwidth for the `i`-th observation
#'  }
#'  }
## -----------------------------EXAMPLES----------------------------------------
#' @examples
#'
#' library(ggplot2)
#'
#' # get the data
#' set.seed(1)
#' data <- generateData()
#' X <- data$data
#' Z <- data$covts
#' interval <- data$interval
#' prec <- data$true_precision
#'
#' # get overall and within interval sample sizes
#' n <- nrow(X)
#' n1 <- sum(interval == 1)
#' n2 <- sum(interval == 2)
#'
#' # visualize the distribution of the extraneous covariate
#' ggplot(data.frame(Z = Z, interval = as.factor(interval))) +
#' geom_histogram(aes(Z, fill = interval), color = "black", bins = n %/% 5)
#'
#' # visualize the true precision matrices in each of the intervals
#'
#' # interval 1
#' matViz(prec[[1]], incl_val = TRUE) +
#' ggtitle("True precision matrix, interval 1")
#'
#' # interval 2 (varies continuously with Z)
#' int2_mats <- prec[interval == 2]
#' int2_inds <- c(5, n2 %/% 2, n2 - 5)
#' lapply(int2_inds, function(j) matViz(int2_mats[[j]], incl_val = TRUE) +
#' ggtitle(paste("True precision matrix, interval 2, observation", j)))
#'
#' # interval 3
#' matViz(prec[[length(prec)]], incl_val = TRUE) +
#' ggtitle("True precision matrix, interval 3")
#'
#' # fit the model and visualize the estimated precision matrices
#' (out <- covdepGE(X, Z))
#' plot(out)
#'
#' # visualize the inclusion probabilities for variables (1, 3) and
#' inclusionCurve(out, 1, 2)
#' inclusionCurve(out, 1, 3)
#'
## -----------------------------DETAILS-----------------------------------------
#' @details
#' # Overview
#'
#' Suppose that `X` is a `p`-dimensional data `matrix` with `n` observations and
#' that `Z` is a `q`-dimensional extraneous covariate, also with `n`
#' observations, where the `l`-th observation in `Z` is associated with the
#' `l`-th observation in `X`. Further suppose that the `l`-th row of `X` follows
#' a `p`-dimensional Gaussian distribution with mean `0` and precision matrix
#' `Omega(z_l)`, where `z_l` is the `l`-th entry of `Z` and `Omega` is a
#' continuous function mapping from the space of extraneous covariates to the
#' space of `p` `x` `p` non-singular matrices. Then, for the `l`-th observation,
#' the `(j,k)` entry of `Omega(z_l)` is non-zero if, and only if, variable `j`
#' and variable `k` are dependent given the remaining variables in `X`.
#'
#' Given data satisfying these assumptions, the `covdepGE` function employs the
#' algorithm described in (1) to estimate a graphical representation of the
#' structure of `Omega` for each of the observations in `X` as a continuous
#' function of `Z`. This graph contains an undirected edge between two variables
#' `X_j` and `X_k` if, and only if, `X_j` and `X_k` are conditionally dependent
#' given the remaining variables.
#'
#' # Graph Estimation
#'
#' Graphs are constructed by fixing each of the columns `X_j` of `X` as the
#' response and performing a spike-and-slab regression using the remaining
#' variables `X_k` in `X` as predictors. To determine if an edge should be added
#' between `X_j` and `X_k`, the posterior inclusion probability of `X_k` in a
#' regression with `X_j` fixed as the response (`PIP(X_k)`) and vice versa
#' (`PIP(X_j)`) are symmetrized according to `sym_method` (e.g., by taking the
#' mean of `PIP(X_j)` and `PIP(X_k)`). If the symmetrized PIP is greater than
#' `edge_threshold`, an edge will be included between `X_j` and `X_k`.
#'
#' To model `Omega` as a function of `Z`, `n` weighted spike-and-slab
#' regressions are performed for each variable `X_j` fixed as the response. The
#' similarity weights for the `l`-th regression are taken with respect to
#' observation `l` such that observations having similar values of `Z` will have
#' larger weights.
#'
#' # Variational Inference
#'
#' Spike-and-slab posterior quantities are estimated using a variational
#' approximation. Coordinate Ascent Variational Inference (CAVI) is performed
#' for each of the weighted regressions to select the variational parameters
#' that maximize the ELBO. The parameters for each of the regression
#' coefficients are the mean and variance of the slab (`mu` and `ssq_var`,
#' respectively) and the probability that the coefficient is non-zero (`alpha`).
#'
#' CAVI for the `n` regressions is performed simultaneously for variable `X_j`
#' fixed as the response. With each of the `n` sets of `alpha` as the rows of an
#' `n` `x` `(p-1)` `matrix`, the CAVI for variable `X_j` is ended for all `n`
#' regressions when the Frobenius norm of the change in the `alpha` `matrix` is
#' less than `alpha_tol` or after `max_iter` iterations of CAVI have been
#' performed.
#'
#' Note that since the regressions performed for variable `X_j` and `X_k` fixed
#' as the response are independent of each other, they may be performed in
#' parallel by setting `parallel = T`. Registering parallel backend with greater
#' than `p` workers offers no benefit, since each worker takes on one variable
#' to fix as the response and perform the `n` regressions.
#'
#' # Hyperparameter specification
#'
#' Each regression requires the specification of 3 hyperparameters: `pi` (the
#' prior probability of inclusion), `sigma^2` (the prior residual variance), and
#' `sigma^2_beta` (the prior variance of the slab). `covdepGE` offers 3 methods
#' for hyperparameter specification via the `hp_method` argument: `grid_search`,
#' `model_average`, and `hybrid`. Empirically, `grid search` offers the best
#' sensitivity and `model_average` offers the best specificity, while `hybrid`
#' sits between the other two methods in both metrics.
#'
#' The hyperparameter candidate grid is generated by taking the Cartesian
#' product between `ssq`, `sbsq`, and `pip` (candidate values for `sigma^2`,
#' `sigma^2_beta`, and `pi`, respectively). Each of the methods gives an
#' approach for selecting points from this grid.
#'
#' In `grid_search`, the point from the grid that produces the model that has
#' the greatest total ELBO is selected, where the total ELBO is calculated by
#' summing the ELBO for each of the `n` regressions for a variable `X_j` fixed
#' as the response. Thus, all observations use the same set of hyperparameters
#' for the regression on `X_j`.
#'
#' Instead of selecting only one model as in `grid_search`, models are averaged
#' over in `model_average`. With `X_j` fixed as the response, the unnormalized
#' weights for each grid point used to perform this averaging is calculated by
#' exponentiating the ELBO for each of the `n` regressions. Note that since the
#' ELBO for a given grid point will vary across the `n` regressions due to
#' differing similarity weights, each of the `n` sets of averaging weights will
#' be unique.
#'
#' Finally, `hybrid` combines `grid_search` and `model_average`. Fixing `X_j` as
#' the response, for each `pi` candidate in `pip`, the point in the grid defined
#' by the Cartesian product of `ssq` and `sbsq` is selected by maximizing the
#' total ELBO summed across the `n` regressions. The resulting models for each
#' of the `pi` candidates are then averaged using the exponentiated ELBO for
#' each of the `n` regressions as the unnormalized averaging weights.
#'
#' Note that in the search step of `grid_search` and `hybrid`, CAVI for each of
#' the candidates is performed for at most `max_iter_grid` iterations. A second
#' CAVI is then performed for `max_iter` iterations using the `n` models that
#' maximized the total ELBO in the first step. Setting `max_iter_grid` to be
#' less than `max_iter` (as is the default) will result in a more efficient
#' search.
#'
#' ## Candidate grid generation
#'
#' The candidate grids (`ssq`, `sbsq`, and `pip`) may be passed as arguments,
#' however, by default, these grids are generated automatically. Each of the
#' grids are spaced uniformly between an upper end point and a lower end point.
#' The number of points in each grid is `5` by default. Grids include end
#' points, and the number of points in each grid is controlled by the arguments
#' `nssq`, `nsbsq`, and `npip`. The lower endpoints (`ssq_lower`, `sbsq_lower`,
#' and `pip_lower`) are all `1e-5` by default. The upper endpoints are
#' calculated dependent on the variable `X_j` fixed as the response.
#'
#' `ssq_upper` is simply the variance of `X_j` times `ssq_mult`. By default,
#' `ssq_mult` is `1.5`.
#'
#' `pip_upper` is calculated by regressing the remaining variables on `X_j`
#' using LASSO. The shrinkage hyperparameter for LASSO is chosen to be
#' `lambda.1se`. The number of non-zero coefficients estimated by LASSO is then
#' divided by `p-1` to calculate `pip_upper`. Note that if the LASSO estimate to
#' the number of non-zero coefficients is `0` or `p-1`, this estimate is changed
#' to `1` or `p-2` (respectively) to ensure that `pip_upper` is greater than `0`
#' and less than `1`.
#'
#' Finally, an upper bound is induced on `sigma^2_beta` by deriving a rough
#' upper bound for the signal-to-noise ratio that depends on `sigma^2_beta`. Let
#' `sum_S^2` be the sum of the sample variances of the columns of the predictors
#' `X’`. Under the simplifying assumptions that the expected values of `X’` and
#' the spike-and-slab regression coefficients `beta` are `0` and that `X’` and
#' `beta` are independent, the variance of the dot product of `X’` with `beta`
#' is `pi*sigma^2*sigma^2_beta*sum_S^2`. Thus, the signal-to-noise ratio under
#' these assumptions is given by `pi*sigma^2_beta*sum_S^2`, and so replacing
#' `pi` with `pip_upper` and `sigma^2_beta` with `sbsq_upper` gives an upper
#' bound on the signal to noise ratio. Setting this bound equal to `snr_upper`
#' gives an expression for `sbsq_upper`.
#'
#' # Weights
#'
#' The similarity weight for individual `k` with respect to individual `l` is
#' calculated as `dnorm(Norm(z_l, z_k), tau_l)`, where `Norm(z_l, z_k)` denotes
#' the norm specified by the `norm` argument of the values of `Z` for the `l`-th
#' and `k`-th observations, and `tau_l` is the bandwidth for the `l`-th
#' observation. `tau` may be passed as an argument, however, by default, it is
#' estimated using the methodology given in (2). (2) describes a two-step
#' approach for density estimation, where in the first step, an initial estimate
#' is calculated using Silverman’s rule of thumb for initializing bandwidth
#' values, and in the second step, the density is refined by updating the
#' bandwidth values. This methodology is used here to estimate the density of
#' `Z`, and the updated bandwidths from the second step are used for `tau`.
## -----------------------------REFERENCES--------------------------------------
#' @references
#' 1. Dasgupta S., Ghosh P., Pati D., Mallick B., *An approximate Bayesian
#' approach to covariate dependent graphical modeling*, 2021
#'
#' 2. Dasgupta S., Pati D., Srivastava A., *A Two-Step Geometric Framework For
#' Density Modeling*, Statistica Sinica, 2020
## -----------------------------------------------------------------------------
covdepGE <- function(X, Z, hp_method = "hybrid", ssq = NULL, sbsq = NULL,
                     pip = NULL, nssq = 5, nsbsq = 5, npip = 5, ssq_mult = 1.5,
                     ssq_lower = 1e-5, snr_upper = 25, sbsq_lower = 1e-5,
                     pip_lower = 1e-5, pip_upper = NULL, tau = NULL, norm = 2,
                     center_X = T, scale_Z = T, alpha_tol = 1e-5,
                     max_iter_grid = 10, max_iter = 100, edge_threshold = 0.5,
                     sym_method = "mean", parallel = F, num_workers = NULL,
                     prog_bar = T){

  start_time <- Sys.time()

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
  if (scale_Z) Z <- matrix(scale(Z)[ , ], n)

  # if X should be centered, do so
  if (center_X) X <- matrix(scale(X, T, F)[ , ], n)

  # get weights
  D <- get_weights(Z, norm, tau)
  bandwidths <- D$bandwidths
  D <- D$D

  # create list for the weights and bandwidths
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

            # Set variable number `j` as the response
            y <- X[, j]

            # Set the remaining p variables as predictors
            X_j <- X[, -j, drop = F]

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

      # Set variable number `j` as the response
      y <- X[, j]

      # Set the remaining p variables as predictors
      X_j <- X[, -j, drop = F]

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
  incl_probs <- replicate(n, matrix(0, p, p), simplify = F)

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

  # create a `list` for graphs and inclusion probabilties
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
