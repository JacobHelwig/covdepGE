## -----------------------------------------------------------------------------
#' var_updates
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' Performs variational updates for a fixed data matrix and response for n
#' linear regressions, where the l-th regression is weighted with respect to the
#' l-th individual. Returns matrix of posterior inclusion probabilities, details
#' on the variational updates and a vector of warnings
## -----------------------------ARGUMENTS---------------------------------------
#' X_mat: n x p matrix; predictors
#'
#' Z: n x p' matrix; extraneous covariates
#'
#' D: n x n matrix; weights (k,l entry is the weight of the k-th individual
#' with respect to the l-th individual using the l-th individual's bandwidth)
#'
#' y: n x 1 vector; response
#'
#' alpha: scalar in [0, 1]; global initialization value for the variational
#' parameters alpha_matrices (approximates probabilities of inclusion). 0.2 by
#' default
#'
#' mu: scalar; global initialization value for the variational parameters
#' mu_matrices (approximates regression coefficients). 0 by default
#'
#' sigmasq: scalar in (0, Inf); Error term variance for spike-and-slab.
#' Algorithm scales this value by individual-specific weights. 0.5 by default
#'
#' sigmabetasq_vec: n_sigma x 1 vector, entries in (0, Inf); candidate values
#' of sigmabeta_sq, the slab variance. NULL by default
#'
#' pi_vec: n_pi x 1 vector, entries in [0, 1]; candidate values of pi. 0.1 by
#' default
#'
#' tolerance: scalar in (0, Inf); end the variational update loop when the
#' square root of the sum of squared changes to the elements of the alpha matrix
#' are within tolerance. 1e-12 by default
#'
#' max_iter: scalar in {1, 2,...; if the tolerance criteria has not been met by
#' max_iter iterations, end the variational update loop. 1e4 by default
#'
#' monitor_final_elbo: logical scalar; if T, the ELBO history for the final
#' model will be returned. F by default
#'
#' monitor_cand_elbo: logical scalar; if T, the ELBO history for each of the
#' candidate models that do not attain convergence within max_iter iterations
#' will be returned
#'
#' monitor_period: scalar in {1, 2,..., max_iter; the periodicity with which the
#' ELBO is recorded if monitor_final_elbo or monitor_cand_elbo is T. 1 by
#' default
#'
#' warnings: logical scalar; if T, convergence and grid warnings will be
#' displayed. Convergence warnings occur when the tolerance exit condition has
#' not been met by max_iter iterations. Grid warnings occur when, for either
#' sigmabetasq_vec or pi_vec, the grid is longer than 2 candidates, and the
#' final model selects a candidate value on the grid boundary. T by default
#'
#' resp_index: scalar; the index of the column that is y in data_mat
## -----------------------------RETURNS-----------------------------------------
#' @return Returns `list` with the following values:
#'
#' 1. alpha_matrix: n x p matrix; the l, j entry is the variational
#' approximation to the posterior inclusion probability of the j-th variable in
#' a regression with the y fixed as the response with weightings taken with
#  respect to the l-th individual
#'
#' 2. VB_details: list with 6 values:
#'  - sigmabeta_sq, pi: scalars; the final values of pi and sigmabeta_sq that
#'  maximized ELBO over all individuals with the j-th predictor fixed as the
#'  response
#'  - ELBO: scalar; the maximum value of ELBO for the final model
#'  - converged_iter: scalar; the number of iterations to attain convergence
#'  for the final model
#'  - ELBO_history: vector; ELBO history by iteration for the final model. If
#'  monitor_final_elbo is F, then this value will be NULL
#'  - non_converged: matrix; each row corresponds to the ELBO history for each
#'  of the candidate models that did not converge. If monitor_cand_elbo is F,
#'  then the ELBO history is omitted, and only the non-convergent sigmabeta_sq
#'  and pi values are provided. If all pairs resulted in convergence, then this
#'  value is NULL
#'
#' 3. warnings_vec: character vector; Vector of convergence warnings to be
#' displayed in covdepGE_main
## -----------------------------------------------------------------------------
var_updates <- function(X_mat, Z, D, y, alpha, mu, sigmasq, sigmabetasq_vec,
                        pi_vec, tolerance, max_iter, monitor_final_elbo,
                        monitor_cand_elbo, monitor_period, warnings, resp_index){

  # get the dimensions of the data
  n <- nrow(X_mat)
  p <- ncol(X_mat)

  # instantiate initial values of variational parameters; the l, j entry is
  # the variational approximation to the j-th parameter in a regression with
  # the resp_index predictor fixed as the response with weightings taken with
  # respect to the l-th individual
  alpha_mat <- matrix(alpha, n, p)
  mu_mat <- matrix(mu, n, p)

  # loop to optimize sigmabeta_sq; for each pair of candidate values of sigma in
  # sigmavec, pi in pi_vec, store the resulting ELBO
  sigma_loop_out <- sigma_loop_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq,
                                 sigmabetasq_vec, pi_vec, tolerance, max_iter,
                                 monitor_cand_elbo, monitor_period)

  # total number of models fit by sigma_loop_c
  total_models <- length(pi_vec) * length(sigmabetasq_vec)

  # get the resulting ELBO and the number of converged models
  elbo_sigmaXpi <- sigma_loop_out[["elbo_grid"]]
  converged <- sigma_loop_out[["num_converged"]]

  # if any of the models did not converge, display a warning
  warning_vec <- c()
  if (converged < total_models & warnings){
    warning_vec <- paste0("Response ", resp_index, ": ",
                          (total_models - converged), "/", total_models,
                          " candidate models did not converge in ", max_iter,
                          " iterations")
  }

  # Select the value of sigma_beta that maximizes the ELBO
  sigmabeta_sq <- sigmabetasq_vec[which(elbo_sigmaXpi
                                        == max(elbo_sigmaXpi), T)[,"row"]][1]

  # Select the value of pi that maximizes the ELBO
  pi_est <- pi_vec[which(elbo_sigmaXpi == max(elbo_sigmaXpi), T)[,"col"]][1]

  # fit another model using these values of sigma_beta and pi_est
  result <- cov_vsvb_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq,
                       pi_est, tolerance, max_iter, monitor_final_elbo,
                       monitor_period)

  # if the final model did not converge, display a warning
  if (result$converged_iter == max_iter & warnings){
    warning_vec <- c(warning_vec, (paste0("Response ", resp_index,
                   ": final model did not converge in ", max_iter,
                   " iterations")))
  }

  # if there are any non-convergent sigmabeta_sq and pi pairs, format these values
  if (nrow(sigma_loop_out$non_converged_pairs) > 0){

    # get the non-converged sigmabeta_sq and pi and format
    nc_ssb_sq <- sigma_loop_out$non_converged_pairs[ , 1]
    nc_pi <- sigma_loop_out$non_converged_pairs[ , 2]
    nc_ssb_sq <- ifelse(round(nc_ssb_sq, 3) == 0, formatC(
      nc_ssb_sq, format = "e", digits = 0), round(nc_ssb_sq, 3))
    nc_pi <- ifelse(round(nc_pi, 3) == 0, formatC(
      nc_pi, format = "e", digits = 0), round(nc_pi, 3))
    formatted_nc_vals <- paste0("slab var: ", nc_ssb_sq, ", pi: ", nc_pi)

    # if the elbo history for non-convergent models has been recorded, use the
    # formatted values as row names for the non_converged_elbo matrix. Otherwise,
    # simply set the non_converged_elbo matrix equal to the formatted values
    if (monitor_cand_elbo){

      # use the first row (iteration index) as the col names; drop the first
      # row
      colnames(sigma_loop_out$non_converged_elbo) <-
        sigma_loop_out$non_converged_elbo[1, ]
      sigma_loop_out$non_converged_elbo <-
        sigma_loop_out$non_converged_elbo[-1, , drop = F]

      # set the row names
      row.names(sigma_loop_out$non_converged_elbo) <- formatted_nc_vals

    }else{
      sigma_loop_out$non_converged_elbo <- formatted_nc_vals

      # if there are no non_converged values, set to NULL
      if(nrow(sigma_loop_out$non_converged_pairs) == 0){
        sigma_loop_out$non_converged_elbo <- NULL
      }
    }
  }else{

    # otherwise, all pairs resulted in convergence; set this value to NULL
    sigma_loop_out$non_converged_elbo <- NULL
  }

  # if the elbo history for the final model has been recorded, set the first
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
  VB_details <- list("sigmabeta_sq" = sigmabeta_sq, "pi" = pi_est,
                     "ELBO" = result$var_elbo,
                     "converged_iter" = result$converged_iter,
                     "ELBO_history" = result$elbo_history,
                     "non_converged" = sigma_loop_out$non_converged_elbo)

  # var.alpha is an n by p matrix; the i,j-th entry is the probability of
  # inclusion for the i-th individual for the j-th variable according to the
  # regression on y
  alpha_matrix <- result$var_alpha

  return(list(alpha_matrix = alpha_matrix, VB_details = VB_details,
              warnings = warning_vec))
}
