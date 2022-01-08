## -----------------------------------------------------------------------------
#' cavi_search
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' Performs CAVI and grid search for a fixed data matrix and response for n
#' linear regressions, where the l-th regression is weighted with respect to the
#' l-th individual. Returns matrix of posterior inclusion probabilities, details
#' on the variational updates, a vector of warnings, and the number of final
#' CAVI that DNC
## -----------------------------ARGUMENTS---------------------------------------
#' X_mat: n x p matrix; predictors
#'
#' Z: n x p' matrix; extraneous covariates
#'
#' D: n x n matrix; weights (k,l entry is the weight of the k-th individual
#' with respect to the l-th individual using the l-th individual's bandwidth)
#'
#' y: vector of length n; response
#'
#' alpha: numeric in [0, 1]; global initialization value for the variational
#' parameters alpha_matrices (approximates probabilities of inclusion). 0.2 by
#' default
#'
#' mu: numeric; global initialization value for the variational parameters
#' mu_matrices (approximates regression coefficients). 0 by default
#'
#' sigmasq: positive numeric; Error term variance for spike-and-slab. Algorithm
#' scales this value by individual-specific weights. 0.5 by default
#'
#' sigmabetasq_vec: vector of length n_sigma with positive entries; candidate
#' values of sigmabeta_sq, the slab variance. NULL by default
#'
#' pi_vec: vector of length n_pi with entries in [0, 1]; candidate values of pi.
#' 0.1 by default
#'
#' tolerance: positive numeric; end CAVI when the Frobenius norm of the
#' iteration-to-iteration change in the alpha matrix are within tolerance.
#' `1e-12` by default
#'
#' max_iter_grid: positive integer; during the grid search, if the
#' tolerance criteria has not been met by `max_iter_grid` iterations, end the
#' CAVI. `1e4` by default
#'
#' max_iter_final: positive integer; for the final CAVI, if the tolerance
#' criteria has not been met by `max_iter_final` iterations, end the CAVI. `1e4`
#' by default
#'
#' monitor_final_elbo: logical; if T, the ELBO history for the final CAVI will
#' be returned. F by default
#'
#' monitor_grid_elbo: logical; if T, the ELBO history for each of the grid
#' points that do not attain convergence within max_iter_grid iterations
#' will be returned
#'
#' monitor_period: integer in
#' \eqn{{1, 2,..., `min(max_iter_grid, max_iter_final)`}}; the periodicity with
#' which the ELBO is recorded if monitor_final_elbo or monitor_grid_elbo is T.
#' 1 by default
#'
#' warnings: logical; if T, convergence and grid warnings will be
#' displayed. Convergence warnings occur when the tolerance exit condition has
#' not been met by max_iter_grid or max_iter_final iterations. Grid warnings
#' occur when, for either sigmabetasq_vec or pi_vec, the grid is longer than 2
#' candidates, and the final CAVI uses a candidate value on the grid boundary.
#' T by default
#'
#' resp_index: integer in {1,...,p + 1}; the index of the column that is y in
#' data_mat
#'
#' CS: logical; if T, pi_vec and sigma_sq will be selected
#' according to Carbonetto-Stephens. F by default
## -----------------------------RETURNS-----------------------------------------
#' Returns `list` with the following values:
#'
#' 1. alpha_matrix: n x p matrix; the l, j entry is the variational
#' approximation to the posterior inclusion probability of the j-th variable in
#' a regression with the y fixed as the response with weightings taken with
#  respect to the l-th individual
#'
#' 2. CAVI_details: list with 6 values:
#'  - sigmabeta_sq, pi: numerics; the grid point that maximized the ELBO for
#'  the j-th variable
#'  - ELBO: numeric; the maximum value of ELBO for the final CAVI
#'  - converged_iter: integer; the number of iterations to attain convergence
#'  for the final CAVI
#'  - ELBO_history: vector; ELBO history by iteration for the final CAVI. If
#'  monitor_final_elbo is F, then this value will be NULL
#'  - non_converged: matrix; each row corresponds to the ELBO history for each
#'  of the grid points that did not converge. If monitor_grid_elbo is F,
#'  then the ELBO history is omitted, and only the non-convergent sigmabeta_sq
#'  and pi values are provided. If all pairs resulted in convergence, then this
#'  value is NULL
#'
#' 3. warnings_vec: character vector; Vector of convergence warnings to be
#' displayed in covdepGE_main
#'
#' 4. final_DNC: integer; number of final CAVIs that did not converge
## -----------------------------------------------------------------------------
cavi_search <- function(X_mat, Z, D, y, alpha, mu, sigmasq, sigmabetasq_vec,
                        pi_vec, tolerance, max_iter_grid, max_iter_final,
                        monitor_final_elbo, monitor_grid_elbo, monitor_period,
                        warnings, resp_index, CS){

  # get the dimensions of the data
  n <- nrow(X_mat)
  p <- ncol(X_mat)

  # instantiate initial values of variational parameters; the l, j entry is
  # the variational approximation to the j-th parameter in a regression with
  # the resp_index predictor fixed as the response with weightings taken with
  # respect to the l-th individual
  alpha_mat <- matrix(alpha, n, p)
  mu_mat <- matrix(mu, n, p)

  # If CS, choose pi and sigmasq according to the Carbonetto-Stephens model
  if (CS){
    set.seed(resp_index)
    idmod <- varbvs::varbvs(X_mat, y, Z = Z[ , 1], verbose = FALSE)
    sigmasq <- mean(idmod$sigma)
    pi_vec <- mean(1 / (1 + exp(-idmod$logodds))) # need to convert to log base 10
  }

  # loop to optimize sigmabeta_sq; run CAVI for each grid points; store the
  # resulting ELBO
  grid_search_out <- grid_search_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq,
                                   sigmabetasq_vec, pi_vec, tolerance,
                                   max_iter_grid, monitor_grid_elbo,
                                   monitor_period)

  # total number of grid points
  grid_size <- length(pi_vec) * length(sigmabetasq_vec)

  # get the resulting ELBO and the number of converged candidates
  elbo_sigmaXpi <- grid_search_out[["elbo_grid"]]
  converged <- grid_search_out[["num_converged"]]

  # if any of the cavi did not converge, display a warning
  warning_vec <- c()
  if (converged < grid_size & warnings){
    warning_vec <- paste0("Variable ", resp_index,
                          ": CAVI did not converge in ", max_iter_grid,
                          " iterations for ", (grid_size - converged),
                          "/", grid_size, " grid search candidates")
  }

  # Select the value of sigma_beta that maximizes the ELBO
  sigmabeta_sq <- sigmabetasq_vec[which(elbo_sigmaXpi
                                        == max(elbo_sigmaXpi), T)[,"row"]][1]

  # Select the value of pi that maximizes the ELBO
  pi_est <- pi_vec[which(elbo_sigmaXpi == max(elbo_sigmaXpi), T)[,"col"]][1]

  # run CAVI using these values of sigmabeta_sq and pi_est
  result <- cov_vsvb_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq,
                       pi_est, tolerance, max_iter_final, monitor_final_elbo,
                       monitor_period)

  # count 1 if the final CAVI DNC
  final_dnc <- 0

  # if the final CAVI did not converge, display a warning
  if (result$converged_iter == max_iter_final & warnings){
    warning_vec <- c(warning_vec, (paste0("Variable ", resp_index,
                   ": final CAVI did not converge in ", max_iter_final,
                   " iterations")))
  }else if (result$converged_iter == max_iter_final){
    final_dnc <- 1
  }

  # if there are any non-convergent grid points, format these values
  if (nrow(grid_search_out$non_converged_pairs) > 0){

    # get the non-converged sigmabeta_sq and pi and format
    nc_ssb_sq <- grid_search_out$non_converged_pairs[ , 1]
    nc_pi <- grid_search_out$non_converged_pairs[ , 2]
    nc_ssb_sq <- ifelse(round(nc_ssb_sq, 3) == 0, formatC(
      nc_ssb_sq, format = "e", digits = 0), round(nc_ssb_sq, 3))
    nc_pi <- ifelse(round(nc_pi, 3) == 0, formatC(
      nc_pi, format = "e", digits = 0), round(nc_pi, 3))
    formatted_nc_vals <- paste0("slab var: ", nc_ssb_sq, ", pi: ", nc_pi)

    # if the elbo history for non-convergent CAVI has been recorded, use the
    # formatted values as row names for the non_converged_elbo matrix. Otherwise,
    # simply set the non_converged_elbo matrix equal to the formatted values
    if (monitor_grid_elbo){

      # use the first row (iteration index) as the col names; drop the first
      # row
      colnames(grid_search_out$non_converged_elbo) <-
        grid_search_out$non_converged_elbo[1, ]
      grid_search_out$non_converged_elbo <-
        grid_search_out$non_converged_elbo[-1, , drop = F]

      # set the row names
      row.names(grid_search_out$non_converged_elbo) <- formatted_nc_vals

    }else{
      grid_search_out$non_converged_elbo <- formatted_nc_vals

      # if there are no non_converged values, set to NULL
      if(nrow(grid_search_out$non_converged_pairs) == 0){
        grid_search_out$non_converged_elbo <- NULL
      }
    }
  }else{

    # otherwise, all pairs resulted in convergence; set this value to NULL
    grid_search_out$non_converged_elbo <- NULL
  }

  # if the elbo history for the final CAVI has been recorded, set the first
  # column (iteration index) to the column names and remove it; then, add
  # a rowname to denote ELBO
  if (monitor_final_elbo){
    colnames(result$elbo_history) <- result$elbo_history[1, ]
    result$elbo_history <- result$elbo_history[-1, , drop = F]
    rownames(result$elbo_history) <- "ELBO"
  }else{

    # otherwise, set this value to NULL
    result$elbo_history <- NULL
  }

  # save the variational bayes details
  cavi_details <- list("sigmabeta_sq" = sigmabeta_sq, "pi" = pi_est,
                       "ELBO" = result$var_elbo,
                       "converged_iter" = result$converged_iter,
                       "ELBO_history" = result$elbo_history,
                       "non_converged" = grid_search_out$non_converged_elbo)

  # var.alpha is an n by p matrix; the i,j-th entry is the probability of
  # inclusion for the i-th individual for the j-th variable according to the
  # regression on y
  alpha_matrix <- result$var_alpha

  return(list(alpha_matrix = alpha_matrix, cavi_details = cavi_details,
              warnings = warning_vec, final_dnc = final_dnc))
}
