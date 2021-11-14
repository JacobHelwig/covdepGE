## _____________________________________________________________________________
## _____________________________covdepGE________________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## Function to estimate the covariance structure for n individuals using
## variational Bayes
## Optimizes the spike and slab parameter sigma_beta_sq by choosing the value
## for each predictor that maximizes the ELBO; can also optimize pi if candidate
## values are specified
## -----------------------------ARGUMENTS---------------------------------------
## data_mat: n by p data matrix
## Z: n by p' matrix of extraneous covariates
## tau: bandwidth parameter; greater values allow for more information to be
## shared from other individuals when estimating the graph for a fixed
## individual. 0.1 by default.
## alpha: global initialization value for the variational parameter alpha
## (approximates probability of inclusion)
## mu: global initialization value for the variational parameter mu
## (approximates regression coefficients)
## sigmavec: candidate values of sigmabeta_sq
## pivec: candidate values of pi. If NULL, uses varsbvs to generate scaler pi
## scale: boolean that dictates whether the extraneous covariates should be
## centered and scaled prior to calculating the weights
## tolerance: end the variational update loop when the square root of the sum of
## squared changes to the elements of the alpha matrix are within tolerance
## max_iter: if the tolerance criteria has not been met by max_iter iterations,
## end the variational update loop
## -----------------------------RETURNS-----------------------------------------
## TBD
## -----------------------------TODO--------------------------------------------
## 1. symmetrization method - mean, min, max
## 2. norm - l2, l1, linf?
## 3. Change alpha matrix return to return asymmetric inclusion probabilties
## _____________________________________________________________________________
#' Title
#'
#' @param data_mat
#' @param Z
#' @param tau
#' @param alpha
#' @param mu
#' @param sigmavec
#' @param pi_vec
#' @param scale
#' @param tolerance
#' @param max_iter
#' @param edge_threshold
#' @param print_time
#'
#' @return
#' @export
#'
#' @examples
covdepGE <- function(data_mat, Z, tau = 0.1, alpha = 0.2, mu = 0,
                      sigmavec = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10),
                      pi_vec = seq(0.1, 0.9, 0.1), scale = T, tolerance = 1e-9,
                      max_iter = 100, edge_threshold = 0.5, print_time = F){

  start_time <- Sys.time()

  # get sample size and number of parameters
  n <- nrow(data_mat); p <- ncol(data_mat) - 1

  # if the covariates should be centered and scaled, do so
  if (scale) Z <- matrix(scale(Z)[ , ], n)

  # D is a symmetric n by n matrix of weights; the i, j entry is the similarity
  # between individuals i and j
  D <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in i:n) {

      # take the l2 norm
      diff_vec <- Z[i, ] - Z[j, ]
      diff_norm <- sqrt(as.numeric(crossprod(diff_vec)))

      # given the l2 norm, find the weight using the normal density
      D[j, i] <- stats::dnorm(diff_norm, 0, tau)
      D[i, j] <- D[j, i]
    }
  }

  # Scale weights to sum to n
  D <- n * (D) * matrix(1 / colSums(D), n, n, T)

  # List for the variable-specific inclusion probability matrix; the j-th element
  # in the list is a n by p matrix; in this matrix, the l-th row corresponds to
  # the probabilties of inclusion for the l-th individual, with the j-th
  # predictor fixed as the response
  alpha_matrices <- vector("list", p + 1)

  # List for saving the final ELBO for each of the p responses
  ELBO_p <- vector("list", p + 1)
  names(ELBO_p) <- paste("Response", 1:(p + 1))

  # main loop over the predictors
  for (resp_index in 1:(p + 1)) {

    # Set variable number `resp_index` as the response
    y <- data_mat[, resp_index]

    # Set the remaining p variables as predictors
    X_mat <- data_mat[, -resp_index]

    # instantiate initial values of variational parameters
    alpha_mat <- matrix(alpha, n, p)
    mu_mat <- matrix(mu, n, p)

    E <- rnorm(n, 0, 1) # removing this causes discrepency in discrete case

    # Setting hyperparameter values for sigmasq and the probability of inclusion
    # according to the Carbonetto-Stephens model
    idmod <- varbvs::varbvs(X_mat, y, Z = Z[ , 1], verbose = FALSE)
    sigmasq <- mean(idmod$sigma)
    if (is.null(pi_vec)){
      pi_vec <- mean(1 / (1 + exp(-idmod$logodds))) # need to convert to log base 10
    }

    # loop to optimize sigma; for each pair of candidate values of sigma in
    # sigmavec, pi in pi_vec, store the resulting ELBO
    elbo_sigmaXpi <- sigma_loop_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq,
                                  sigmavec, pi_vec, tolerance, max_iter)

    # Select the value of sigma_beta that maximizes the ELBO
    sigmabeta_sq <- sigmavec[which(elbo_sigmaXpi == max(elbo_sigmaXpi), T)[,"row"]]

    # Select the value of pi that maximizes the ELBO
    pi_est <- pi_vec[which(elbo_sigmaXpi == max(elbo_sigmaXpi), T)[,"col"]]

    # fit another model using these values of sigma_beta and pi_est
    result <- cov_vsvb_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq,
                         pi_est, tolerance, max_iter)

    # save the final ELBO
    ELBO_p[[resp_index]] <- list("sigma^2_beta" = sigmabeta_sq, "pi" = pi_est,
                                 "ELBO" = result$var.elbo)

    # var.alpha is an n by p matrix; the i,j-th entry is the probability of
    # inclusion for the i-th individual for the j-th variable according to the
    # regression on y
    alpha_matrices[[resp_index]] <- result$var.alpha
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

  # symmetrize the inclusion matrices
  incl_probs_asym <- incl_probs # return this instead of the alpha matrices
  incl_probs <- lapply(incl_probs, function(mat) (mat + t(mat)) / 2)

  # if the probability of an edge is greater than edge_threshold, denote an
  # edge by a 1; otherwise, 0
  graphs <- lapply(incl_probs, function(mat) (mat > edge_threshold) * 1)

  # stop timer and see how much time has elapsed
  if (print_time) print(Sys.time() - start_time)

  return(list(graphs = graphs, inclusion_probs = incl_probs,
              alpha_matrices = alpha_matrices, ELBO = ELBO_p))
}
