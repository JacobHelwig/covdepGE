## _____________________________________________________________________________
## _____________________________covdepGE________________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## Function to estimate the covariance structure for n individuals using
## variational Bayes
## Optimizes the spike and slab paramter sigma_beta_sq by choosing the value
## for each predictor that maximizes the ELBO
## -----------------------------ARGUMENTS---------------------------------------
## data_mat: n by p data matrix
## Z: n by p' matrix of extraneous covariates
## tau: bandwidth parameter; greater values allow for more information to be
## shared from other individuals when estimating the graph for a fixed
## individual. 0.1 by default.
## default
## sigmavec: candidate values of sigmabeta_sq
## tolerance: end the variational update loop when the square root of the sum of
## squared changes to the elements of the alpha matrix are within tolerance
## max_iter: if the tolerance criteria has not been met by max_iter iterations,
## end the variational update loop
## -----------------------------RETURNS-----------------------------------------
## TBD
## _____________________________________________________________________________
#' covdepGE
#'
#' @param data_mat
#' @param Z
#' @param tau
#' @param sigmavec
#' @param tolerance
#' @param max_iter
#'
#' @return
#' @export
#'
#' @examples
covdepGE <- function(data_mat, Z, tau = 0.1,
                     sigmavec = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10),
                     tolerance = 1e-9, max_iter = 100, print_time = F){

  start_time <- Sys.time()

  # get sample size and number of parameters
  n <- nrow(data_mat); p <- ncol(data_mat) - 1

  # D is a symmetric n by n matrix of weights; the i, j entry is the similarity
  # between individuals i and j
  D <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      # this is the spectral norm
      D[j, i] <- stats::dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, tau)
      D[i, j] <- D[j, i]
    }
  }

  # Scale weights to sum to n
  D <- n * (D) * matrix(1 / colSums(D), n, n, T)

  # List for the variable-specific inclusion probability matrix; the j-th element
  # in the list is a n by p matrix corresponding to the k-th predictor;
  # in this matrix, the l-th row corresponds to the dependence structure for the
  # l-th individual, with the j-th predictor fixed as the response
  graph_list <- vector("list", p + 1)

  # main loop over the predictors
  for (resp_index in 1:(p + 1)) {

    # Set variable number `resp_index` as the response
    y <- data_mat[, resp_index]

    # Set the remaining p variables as predictors
    X_mat <- data_mat[, -resp_index]

    # instantiate initial values of variational parameters and hyperparmeters
    alpha_mat <- matrix(0.2, n, p) # argument?
    sigmabeta_sq <- 3 # argument?
    sigmasq <- 1 # argument?
    E <- rnorm(n, 0, sigmasq) # removing this causes discrepency in discrete case
    # dont delete S_sq <- matrix(sigmasq * (DXtX + 1 / sigmabeta_sq)^(-1), n, p) # should be byrow = T?
    S_sq <- t(sigmasq * (t(X_mat^2) + 1 / sigmabeta_sq)^(-1))
    mu_mat <- matrix(0, n, p, byrow = TRUE) # argument?

    # Setting hyperparameter values for sigmasq and the probability of inclusion
    # according to the Carbonetto Stephens model
    idmod <- varbvs::varbvs(X_mat, y, Z = Z[, 1], verbose = FALSE)
    sigmasq <- mean(idmod$sigma)
    pi_est <- mean(1 / (1 + exp(-idmod$logodds)))

    # loop to optimize sigma; for each value of sigma in sigmavec, store the
    # resulting ELBO
    elbo_sigma <- sigma_loop_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmavec,
                               pi_est, tolerance, max_iter)

    # Select the value of sigma_beta that maximizes the ELBO
    sigmabeta_sq <- sigmavec[which.max(elbo_sigma)]

    # fit another model using this value of sigma_beta
    result <- cov_vsvb_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq,
                         pi_est, tolerance, max_iter)

    # n by p matrix; the i,j-th entry is the probability of inclusion for the
    # i-th individual for the j-th variable according to the regression on y
    graph_list[[resp_index]] <- result$var.alpha
  }

  # stop timer and see how much time has elapsed
  end_time <- Sys.time()
  if (print_time) print(end_time - start_time)

  return(graph_list)
}
