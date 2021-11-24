setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev")
rm(list = ls())
source("generate_data.R")

## _____________________________________________________________________________
## _____________________________covdepGE________________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## Function to model the conditional dependence structure of data as a function
## of individual-specific extranesous covariates as described in (1)
## "An approximate Bayesian approach to covariate dependent graphical modeling"
## -----------------------------ARGUMENTS---------------------------------------
## data_mat: n by (p + 1) matrix; data
## Z: n by p' matrix; extraneous covariates
## tau: scalar in (0, Inf); global bandwidth parameter. greater values allow for
## more information to be shared between individuals. 0.1 by default.
## kde: boolean; if T, use 2-step KDE methodology described in (2) to calculate
## individual-specific bandwidths in place of global bandwidth parameter tau.
## T by default
## alpha: scalar in [0, 1]; global initialization value for the variational
## parameter alpha (approximates probability of inclusion). 0.2 by default
## mu: scalar; global initialization value for the variational parameter mu
## (approximates regression coefficients). 0 by default
## sigmasq: scalar in (0, Inf); variance hyperparameter for spike-and-slab.
## 0.5 by default
## sigmabetasq_vec: n_sigma x 1 vector, entries in (0, Inf); candidate values of
## sigmabeta_sq. NULL by default
## var_min: scalar in (0, Inf); if sigmabetasq_vec is NULL, var_min is the lower
## bound of the auto-generated sigmabetasq_vec. 0.01 by default
## var_max: scalar in (0, Inf); if sigmabetasq_vec is NULL, var_max is the upper
## bound of the auto-generated sigmabetasq_vec. 10 by default
## n_sigma: scalar in {1, 2,...}; if sigmabetasq_vec is NULL, auto-generate it
## as sigmabetasq_vec <- exp(seq(log(var_max), log(var_min), length = n_sigma))
## 8 by default. Ex: for default values, sigmabeta_sq is:
## > round(exp(seq(log(var_max), log(var_min), length = n_sigma)), 3)
## [1] 10.000  3.728  1.389  0.518  0.193  0.072  0.027  0.010
## pi_vec: n_pi x 1 vector; candidate values of pi. 0.2 by default
## norm: scalar in [1, Inf]; norm to use when calculating weights. Inf results
## in infinity norm. 2 by default
## scale: boolean; if T, center and scale extraneous covariates to 0 mean,
## standard deviation 1 prior to calculating the weights. T by default
## tolerance: scalar in (0, Inf); end the variational update loop when the
## square root of the sum of squared changes to the elements of the alpha matrix
## are within tolerance. 1e-9 by default
## max_iter: scalar in {1, 2,...} if the tolerance criteria has not been met by
## max_iter iterations, end the variational update loop. 100 by default
## edge_threshold: scalar in [0, 1]; when processing the inclusion
## probabilities, add an edge to the graph if the (i, j) edge has probability
## of inclusion greater than edge_threshold. 0.5 by default
## sym_method: string in {"mean", "max", "min"}; to symmetrize the alpha
## matrices, the i,j = j,i entry is sym_method((i,j entry), (j,i) entry). "mean"
## by default
## -----------------------------RETURNS-----------------------------------------
## graphs: list of n (p + 1) x (p + 1) matrices; the l-th element is the graph
## for the l-th individual (obtained from inclusion_probs according to
## edge_threshold)
## inclusion_probs: list of n (p + 1) x (p + 1) matrices; the l-th element is a
## symmetric matrix of inclusion probabilities for the l-th individual
## (obtained by symmetrizing alpha_matrices according to sym_method)
## alpha_matrices: list of n (p + 1) x (p + 1) matrices; the l-th element is an
## asymmetric matrix of inclusion probabilities
## ELBO: list of (p + 1) lists; the j-th list corresponds to the j-th predictor
## and contains 3 elements - the final values of pi and sigmabeta_sq that
## maximized ELBO over all individuals with the j-th predictor fixed as the
## response and the maximum value of ELBO
covdepGE1 <- function(data_mat, Z, tau = 0.1, alpha = 0.2, mu = 0, sigmasq = 0.5,
                      sigmabetasq_vec = NULL, var_min = 0.01, var_max = 10,
                      n_sigma = 8, pi_vec = 0.2, norm = 2, scale = T,
                      tolerance = 1e-9, max_iter = 100, edge_threshold = 0.5,
                      sym_method = "mean", print_time = F, CS = F){

  start_time <- Sys.time()

  # get sample size and number of parameters
  n <- nrow(data_mat); p <- ncol(data_mat) - 1

  # if the covariates should be centered and scaled, do so
  if (scale) Z <- matrix(scale(Z)[ , ], n)

  # D is a n by n matrix of weights; the l, k entry is the weighting
  # between individuals of individual k with respect to individual l.
  # If KDE = T, then the bandwidth parameter (tau) used is the bandwidth
  # corresponding to individual l
  D <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in i:n) {

      # take the p-norm
      diff_vec <- Z[i, ] - Z[j, ]

      if (norm == 2){

        # take the 2-norm, use crossprod
        diff_norm <- sqrt(as.numeric(crossprod(diff_vec)))
      }else if (is.infinite(norm)){

        # take the infinity norm
        diff_norm <- max(abs(diff_vec))
      }else{

        # take the p-norm
        diff_norm <- (sum(abs(diff_vec)^norm))^(1 / norm)
      }


      # given the norm, find the weight using the normal density
      D[j, i] <- stats::dnorm(diff_norm, 0, tau)
      D[i, j] <- D[j, i]
    }
  }

  # Scale weights to sum to n down the columns
  D <- n * (D) * matrix(1 / colSums(D), n, n, T)

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

    E <- rnorm(n, 0, 1) # removing this causes discrepency in discrete case

    # If CS, choose pi and sigmasq according to the Carbonetto-Stephens model
    if (CS){
      idmod <- varbvs::varbvs(X_mat, y, Z = Z[ , 1], verbose = FALSE)
      sigmasq <- mean(idmod$sigma)
      pi_vec <- mean(1 / (1 + exp(-idmod$logodds))) # need to convert to log base 10
    }

    # loop to optimize sigmabeta_sq; for each pair of candidate values of sigma in
    # sigmavec, pi in pi_vec, store the resulting ELBO
    elbo_sigmaXpi <- sigma_loop_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq,
                                  sigmabetasq_vec, pi_vec, tolerance, max_iter)

    # Select the value of sigma_beta that maximizes the ELBO
    sigmabeta_sq <- sigmabetasq_vec[which(elbo_sigmaXpi
                                          == max(elbo_sigmaXpi), T)[,"row"]]

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
              alpha_matrices = incl_probs_asym, ELBO = ELBO_p))
}

# generate data and covariates
discrete_data <- F # true if discrete example is desired
if (discrete_data) {
  dat <- generate_discrete()
  tau_ <- 0.1 # the bandwidth parameter
}else{
  dat <- generate_continuous()
  tau_ <- 0.56
}

data_mat <- dat$data
Z.cov <- dat$covts

package <- F # true if the package version is desired
if (package){
  out <- covdepGE::covdepGE(data_mat, Z.cov, tau_, print_time = T, CS = T,
                            scale = F,
                            sigmabetasq_vec = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10))
}else{
  Rcpp::sourceCpp("c_dev.cpp")
  out <- covdepGE1(data_mat, Z.cov, tau_, print_time = T, CS = T, scale = F,
                   sigmabetasq_vec = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10))
}

# check to see that this modified code produces the same results as the original code
if (discrete_data){
  load("out_original_discrete.Rdata")
}else{
  load("out_original_continuous.Rdata")
}

# check for equality between the alpha matrices
same_alpha <- T
for (j in 1:length(out$alpha_matrices)) {
  if (all.equal(out$alpha_matrices[[j]],
                out_original$original_alpha_matrices[[j]]) != T) {
    same_alpha <- F
    break
  }
}
same_alpha

# check for equality between the inclusion probabilities
same_probs <- T
for (j in 1:length(out$inclusion_probs)) {
  if (all.equal(out$inclusion_probs[[j]],
                out_original$original_incl_probs[[j]]) != T) {
    same_probs <- F
    break
  }
}
same_probs

# check for equality between ELBO
all.equal(unname(unlist(lapply(out$ELBO, `[[`, 3))), out_original$original_ELBO)
