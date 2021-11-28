## _____________________________________________________________________________
## _____________________________covdepGE________________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## Function to model the conditional dependence structure of data as a function
## of extraneous covariates as described in (1)
## "An approximate Bayesian approach to covariate dependent graphical modeling".
## The final model uses the candidate pair of pi and sigmabeta_sq that maximize
## the ELBO for each variable
## -----------------------------ARGUMENTS---------------------------------------
## data_mat: n x (p + 1) matrix; data
##
## Z: n x p' matrix; extraneous covariates
##
## tau: scalar in (0, Inf) OR n x 1 vector, entries in (0, Inf); bandwidth
## parameter. Greater values allow for more information to be shared between
## individuals. Allows for global or individual-specific specification. If
## kde = T, this argument is ignored. 0.1 by default
##
## kde: logical scalar; if T, use 2-step KDE methodology as described in (2) to
## calculate individual-specific bandwidths in place of global bandwidth
## parameter tau. T by default
##
## alpha: scalar in [0, 1]; global initialization value for the variational
## parameters alpha (approximates probabilities of inclusion). 0.2 by default
##
## mu: scalar; global initialization value for the variational parameters mu
## (approximates regression coefficients). 0 by default
##
## sigmasq: scalar in (0, Inf); variance hyperparameter for spike-and-slab.
## 0.5 by default
##
## sigmabetasq_vec: n_sigma x 1 vector, entries in (0, Inf); candidate values of
## sigmabeta_sq, the slab variance. NULL by default
##
## var_min: scalar in (0, Inf); if sigmabetasq_vec is NULL, var_min is the lower
## bound of the auto-generated sigmabetasq_vec. 0.01 by default
##
## var_max: scalar in (0, Inf); if sigmabetasq_vec is NULL, var_max is the upper
## bound of the auto-generated sigmabetasq_vec. 10 by default
##
## n_sigma: scalar in {1, 2,...}; if sigmabetasq_vec is NULL, it will be
## auto-generated as:
## sigmabetasq_vec <- exp(seq(log(var_max), log(var_min), length = n_sigma))
## 8 by default.
##
## pi_vec: n_pi x 1 vector, entries in [0, 1]; candidate values of pi. 0.2 by
## default
##
## norm: scalar in [1, Inf]; norm to use when calculating weights. Inf results
## in infinity norm. 2 by default
##
## scale: logical scalar; if T, center and scale extraneous covariates to mean 0,
## standard deviation 1 prior to calculating the weights. T by default
##
## tolerance: scalar in (0, Inf); end the variational update loop when the
## square root of the sum of squared changes to the elements of the alpha matrix
## are within tolerance. 1e-9 by default
##
## max_iter: scalar in {1, 2,...} if the tolerance criteria has not been met by
## max_iter iterations, end the variational update loop. 100 by default
##
## edge_threshold: scalar in (0, 1); when processing the inclusion
## probabilities, an edge will be added to the graph if the (i, j) edge has
## probability of inclusion greater than edge_threshold. 0.5 by default
##
## sym_method: character scalar in {"mean", "max", "min"}; to symmetrize the
## alpha matrices, the i,j = j,i entry is sym_method((i,j entry), (j,i entry)).
## "mean" by default
##
## print_time: logical scalar; if T, function run time is printed. F by default
##
## warnings: logical scalar; if T, convergence and grid warnings will be displayed.
## Convergence warnings occur when the tolerance exit condition has not been
## met by max_iter iterations. Grid warnings occur when, for either
## sigmabetasq_vec or pi_vec, the grid is longer than 2 candidates, and the
## final model selects a candidate value on the grid boundary. T by default
##
## CS: logical scalar; if T, pi_vec and sigma_sq will be scalars selected according to
## Carbonetto-Stephens
##
## -----------------------------RETURNS-----------------------------------------
## graphs: list of n (p + 1) x (p + 1) matrices; the l-th matrix is the
## adjacency matrix for the l-th individual (obtained from inclusion_probs
## according to edge_threshold)
##
## inclusion_probs: list of n (p + 1) x (p + 1) matrices; the l-th element is a
## symmetric matrix of inclusion probabilities for the l-th individual
## (obtained by symmetrizing alpha_matrices according to sym_method)
##
## alpha_matrices: list of n (p + 1) x (p + 1) matrices; the l-th element is an
## asymmetric matrix of inclusion probabilities
##
## ELBO: list of (p + 1) lists; the j-th list corresponds to the j-th predictor
## and contains 3 elements - the final values of pi and sigmabeta_sq that
## maximized ELBO over all individuals with the j-th predictor fixed as the
## response and the maximum value of ELBO
##
## weights: n x n matrix; j, i entry is the weighting of the j-th individual
## with respect to the i-th individual using the i-th individual's bandwidth
##
## bandwidths: n x 1 vector; individual-specific bandwidths
##
#' Title
#'
#' @param data_mat
#' @param Z
#' @param tau
#' @param kde
#' @param alpha
#' @param mu
#' @param sigmasq
#' @param sigmabetasq_vec
#' @param var_min
#' @param var_max
#' @param n_sigma
#' @param pi_vec
#' @param norm
#' @param scale
#' @param tolerance
#' @param max_iter
#' @param edge_threshold
#' @param sym_method
#' @param print_time
#' @param warnings
#' @param CS
#'
#' @return
#' @export
#'
#' @examples
covdepGE <- function(data_mat, Z, tau = 0.1, kde = T, alpha = 0.2, mu = 0,
                     sigmasq = 0.5, sigmabetasq_vec = NULL, var_min = 0.01,
                     var_max = 10, n_sigma = 8, pi_vec = 0.2, norm = 2,
                     scale = T, tolerance = 1e-9, max_iter = 100,
                     edge_threshold = 0.5, sym_method = "mean", print_time = F,
                     warnings = T, CS = F){

  start_time <- Sys.time()

  # run compatibility checks
  covdepGE_checks(data_mat, Z, tau, kde, alpha, mu, sigmasq, sigmabetasq_vec,
                  var_min, var_max, n_sigma, pi_vec, norm, scale, tolerance,
                  max_iter, edge_threshold, sym_method, print_time, warnings)

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

  # List for the variable-specific inclusion probability matrix; the j-th element
  # in the list is a n by p matrix; in this matrix, the l-th row corresponds to
  # the probabilties of inclusion for the l-th individual, with the j-th
  # predictor fixed as the response
  alpha_matrices <- vector("list", p + 1)

  # List for saving the final ELBO for each of the p responses
  ELBO_p <- vector("list", p + 1)
  names(ELBO_p) <- paste("Response", 1:(p + 1))

  # if sigmabetasq_vec is NULL, instantiate the grid
  if(is.null(sigmabetasq_vec)){
    sigmabetasq_vec <- exp(seq(log(var_max), log(var_min), length = n_sigma))
  }

  # main loop over the predictors
  for (resp_index in 1:(p + 1)) {

    # Set variable number `resp_index` as the response
    y <- data_mat[, resp_index]

    # Set the remaining p variables as predictors
    X_mat <- data_mat[, -resp_index]

    # instantiate initial values of variational parameters; the l, j entry is
    # the variational approximation to the j-th parameter in a regression with
    # the resp_index predictor fixed as the response with weightings taken with
    # respect to the l-th individual
    alpha_mat <- matrix(alpha, n, p)
    mu_mat <- matrix(mu, n, p)

    E <- stats::rnorm(n, 0, 1) # removing this causes discrepency in discrete case

    # If CS, choose pi and sigmasq according to the Carbonetto-Stephens model
    if (CS){
      idmod <- varbvs::varbvs(X_mat, y, Z = Z[ , 1], verbose = FALSE)
      sigmasq <- mean(idmod$sigma)
      pi_vec <- mean(1 / (1 + exp(-idmod$logodds))) # need to convert to log base 10
    }

    # loop to optimize sigmabeta_sq; for each pair of candidate values of sigma in
    # sigmavec, pi in pi_vec, store the resulting ELBO
    sigma_loop_out <- sigma_loop_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq,
                                   sigmabetasq_vec, pi_vec, tolerance, max_iter)

    # total number of models fit by sigma_loop_c
    total_models <- length(pi_vec) * length(sigmabetasq_vec)

    # get the resulting ELBO and the number of converged models
    elbo_sigmaXpi <- sigma_loop_out[["elbo_grid"]]
    converged <- sigma_loop_out[["num_converged"]]

    # if any of the models did not converge, display a warning
    if (converged < total_models & warnings){
      warning(paste0("Response ", resp_index, ": ", (total_models - converged),
                     "/", total_models, " candidate models did not converge in ",
                     max_iter, " iterations"))
    }

    # Select the value of sigma_beta that maximizes the ELBO
    sigmabeta_sq <- sigmabetasq_vec[which(elbo_sigmaXpi
                                          == max(elbo_sigmaXpi), T)[,"row"]]

    # Select the value of pi that maximizes the ELBO
    pi_est <- pi_vec[which(elbo_sigmaXpi == max(elbo_sigmaXpi), T)[,"col"]]

    # fit another model using these values of sigma_beta and pi_est
    result <- cov_vsvb_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq,
                         pi_est, tolerance, max_iter)

    # if the final model did not converge, display a warning
    if (!result$converged & warnings){
      warning(paste0("Response ", resp_index, ": final model did not converge in ",
                     max_iter, " iterations"))
    }

    # save the final ELBO
    ELBO_p[[resp_index]] <- list("sigma^2_beta" = sigmabeta_sq, "pi" = pi_est,
                                 "ELBO" = result$var.elbo)

    # var.alpha is an n by p matrix; the i,j-th entry is the probability of
    # inclusion for the i-th individual for the j-th variable according to the
    # regression on y
    alpha_matrices[[resp_index]] <- result$var.alpha
  }

  if (warnings){

    # grid warnings - for each of the candidate grids (pi and sigmabeta_sq), if
    # points along the grid boundaries were selected and the grid had more than
    # 2 candidates, display a warning

    # sigmabetasq_vec
    if (length(sigmabetasq_vec) > 2){

      # get the selected values of sigmabeta_sq for each response
      final_sigmabeta_sq <- unlist(lapply(ELBO_p, `[[`, 1))

      # count the number of final_sigmabeta_sq that were on the boundary of the
      # grid
      grid_boundary <- sigmabetasq_vec[c(1, length(sigmabetasq_vec))]
      on_boundary <- sum(final_sigmabeta_sq %in% grid_boundary)

      # if any of the final sigma were on the boundary, display a warning
      if (on_boundary > 0){
        warning(paste0("For ", on_boundary, "/", p + 1,
                       " responses, the selected value of sigmabeta_sq was on the grid boundary. See return value ELBO for details"))
      }
    }

    # pi_vec
    if (length(pi_vec) > 2){

      # get the selected values of pi for each response
      final_pi <- unlist(lapply(ELBO_p, `[[`, 2))

      # count the number of final_pi that were on the boundary of the grid
      grid_boundary <- pi_vec[c(1, length(pi_vec))]
      on_boundary <- sum(final_pi %in% grid_boundary)

      # if any of the final pi were on the boundary, display a warning
      if (on_boundary > 0){
        warning(paste0("For ", on_boundary, "/", p + 1,
                       " responses, the selected value of pi was on the grid boundary. See return value ELBO for details"))
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

  # stop timer and see how much time has elapsed
  if (print_time) print(Sys.time() - start_time)

  return(list(graphs = graphs, inclusion_probs = incl_probs,
              alpha_matrices = incl_probs_asym, ELBO = ELBO_p, weights = D,
              bandwidths = bandwidths))
}
