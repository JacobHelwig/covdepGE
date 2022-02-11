source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/cavi_search.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/weights.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/covdepGE_R.R")
Rcpp::sourceCpp("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/src/covdepGE_c.cpp")

# function for generating the data and the covariates
generate_continuous <- function(n1 = 60, n2 = 60, n3 = 60, p = 4){

  # create covariate for individuals in each of the three intervals

  # define the dimensions of the data
  n <- sum(n1, n2, n3)

  # define the limits of the intervals
  limits1 <- c(-3, -1)
  limits2 <- c(-1, 1)
  limits3 <- c(1, 3)

  # define the covariate values within each interval
  z1 <- runif(n1, limits1[1], limits1[2])
  z2 <- runif(n2, limits2[1], limits2[2])
  z3 <- runif(n3, limits3[1], limits3[2])
  # z1 <- seq(limits1[1], limits1[2], length = n1)
  # z2 <- seq(limits2[1], limits2[2], length = n2)
  # z3 <- seq(limits3[1], limits3[2], length = n3)
  Z <- matrix(sort(c(z1, z2, z3)), n, 1)

  # create precision matrices

  # the shared part of the structure for all three intervals is a 2 on the
  # diagonal and a 1 in the (2, 3) position
  common_str <- diag(p + 1)
  common_str[2, 3] <- 1

  # define constants for the structure of interval 2
  beta1 <- diff(limits2)^-1
  beta0 <- -limits2[1] * beta1

  # interval 2 has two different linear functions of Z in the (1, 2) position
  # and (1, 3) positions; define structures for each of these components
  int2_str12 <- int2_str13 <- matrix(0, p + 1, p + 1)
  int2_str12[1, 2] <- int2_str13[1, 3] <- 1

  # define the precision matrices for each of the individuals in interval 2
  int2_prec <- lapply(z2, function(z) common_str +
                        ((1 - beta0 - beta1*z)*int2_str12) +
                        ((beta0 + beta1*z)*int2_str13))

  # interval 1 has a 1 in the (1, 2) and interval 3 has a 1 in the (1, 3) position;
  # define structures for each of these components
  int1_str12 <- int3_str13 <- matrix(0, p + 1, p + 1)
  int1_str12[1, 2] <- int3_str13[1, 3] <- 1

  # define the precision matrices for each of the individuals in interval 1 and interval 3
  int1_prec <- rep(list(common_str + int1_str12), n1)
  int3_prec <- rep(list(common_str + int3_str13), n3)

  # put all of the precision matrices into one list
  prec_mats <- c(int1_prec, int2_prec, int3_prec)

  # symmetrize the precision matrices
  prec_mats <- lapply(prec_mats, function(mat) t(mat) + mat)

  # invert the precision matrices to get the covariance matrices
  cov_mats <- lapply(prec_mats, solve)

  # generate the data using the covariance matrices
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))

  return(list(data = data_mat, covts = Z, true_precision = prec_mats))
}

# colors for plots
set.seed(1)
colors <- c("chartreuse3", "chocolate2", "cornflowerblue", "darkgoldenrod1",
            "darkmagenta", "deepskyblue3", "forestgreen", "darkorchid3",
            "darkred", "darkslategray")
colors <- c(colors, sample(colors()[sapply(colors(), function(color)
  !(substr(color, 1, 4)) %in% c("grey", "gray"))], 180))

cont <- generate_continuous(p = 25)
data_mat <- cont$data
Z <- cont$covts
n <- nrow(data_mat)
p <- ncol(data_mat) - 1


# if the covariates should be centered and scaled, do so ([ , ] for attributes)
Z <- matrix(scale(Z)[ , ], n)

kde <- T; tau <- 0.5; norm <- 2
# get weights
D <- get_weights(Z, norm, kde, tau)
bandwidths <- D$bandwidths
D <- D$D

n_param <- 9
pi_vec <- sigmasq_vec <- sigmabetasq_vec <- update_sigmasq <- update_sigmabetasq <- NULL

# if the user has not provided values for pi_vec, instantiate the grid
if (is.null(pi_vec)){

  pi_vec <- seq(0.05, 0.45, length.out = n_param)
}

# instantiate the matrices of hyperparameters
if (is.null(sigmasq_vec)){

  # if no values have been passed to sigmasq_vec, instantiate as the variance
  # of the data
  sigmasq_vec <- matrix(var(as.vector(data_mat)), n, n_param)

  if (is.null(update_sigmasq)){

    # if the user has also not specified an option for updating sigmasq_vec,
    # update sigmasq
    update_sigmasq <- T
  }
}else{

  # otherwise, the user has passed some values for sigmasq; create a matrix
  # where the j-th column is the j-th value repeated n times
  sigmasq_vec <- matrix(sigmasq_vec, n, n_param, T)

  if (is.null(update_sigmasq)){

    # if the user has not specified an option for updating sigmasq_vec but has
    # provided some values to sigmavec, do not update them
    update_sigmasq <- F
  }
}
if (is.null(sigmabetasq_vec)){

  # if no values have been passed to sigmabetasq_vec, instantiate as 1
  sigmabetasq_vec <- matrix(1, n, n_param)

  if (is.null(update_sigmabetasq)){

    # if the user has also not specified an option for updating
    # sigmabetasq_vec, update sigmabetasq
    update_sigmabetasq <- T
  }
}else{

  # otherwise, the user has passed some values for sigmabetasq; create a
  # matrix where the j-th column is the j-th value repeated n times
  sigmabetasq_vec <- matrix(sigmabetasq_vec, n, n_param, T)

  if (is.null(update_sigmabetasq)){

    # if the user has not specified an option for updating sigmabetasq_vec but
    # has provided some values to sigmabetasq_vec, do not update them
    update_sigmabetasq <- F
  }
}

# main loop over the predictors

# otherwise, CAVI will be executed sequentially


resp_index <- 15
alpha <- 0.2
mu <- 0
tolerance <- 1e-12
max_iter <- 1000
warnings = T
CS <- F
R <- T

# Set variable number `resp_index` as the response
y <- data_mat[, resp_index]

# Set the remaining p variables as predictors
X_mat <- data_mat[, -resp_index]

# perform the grid search and final CAVI; save the results to res
res <- cavi_search(X_mat, Z, D, y, alpha, mu, sigmasq_vec,
                                 update_sigmasq, sigmabetasq_vec,
                                 update_sigmabetasq, pi_vec, tolerance,
                                 max_iter, warnings, resp_index, CS, R)

mu_mat <- matrix(mu, n, p)
alpha_mat <- matrix(alpha, n, p)
sigmasq <- sigmasq_vec[ , 1]
sigmabeta_sq <- sigmabetasq_vec[ , 1]
pi_est <- 0.45
cavi_R(y, D, X_mat, mu_mat, alpha_mat, sigmasq, update_sigmasq,
        sigmabeta_sq, update_sigmabetasq, pi_est, tolerance,
        max_iter, upper_limit = 9)

