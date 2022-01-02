setwd("~/Jacob/covdepGE/dev")
#setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev")
rm(list = ls())

package <- F # true if the package version is desired
if (package){
  rm(list = ls())
  library(covdepGE)
}else{
  rm(list = ls())
  if ("covdepGE" %in% .packages()) detach("package:covdepGE", unload = TRUE)
  source("~/Jacob/covdepGE/R/covdepGE_main.R")
  source("~/Jacob/covdepGE/R/weights.R")
  source("~/Jacob/covdepGE/R/checks.R")
  source("~/Jacob/covdepGE/R/gg_covdepGE.R")
  # source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/weights.R")
  # source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/checks.R")
  # source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/gg_covdepGE.R")
  Rcpp::sourceCpp("~/Jacob/covdepGE/src/covdepGE_c.cpp")
}

source("generate_data.R")

cont <- generate_continuous()
data_mat <- cont$data
Z <- cont$covts


## -----------------------------------------------------------------------------
## -----------------------------bug correction----------------------------------
## -----------------------------------------------------------------------------

# BUG: named matrices to gg_adjMat
mat <- matrix(1:9, 3)
colnames(mat) <- row.names(mat) <- letters[1:3]
gg_adjMat(mat)

mat <- matrix(1:9, 3)
colnames(mat) <- letters[1:3]
gg_adjMat(mat)

# BUG: too many sigma
set.seed(1)
n_neg <- 20
n_pos <- 80
n <- n_neg + n_pos
p <- 4

# generate the extraneous covariate
Z_neg <- sort(-runif(n_neg))
Z_pos <- sort(runif(n_pos))
Z <- c(Z_neg, Z_pos)

# create true covariance structure for 2 groups: positive Z and negative Z
true_graph_pos <- true_graph_neg <- matrix(0, p + 1, p + 1)
true_graph_pos[1, 2] <- true_graph_pos[2, 1] <- true_graph_neg[1, 3] <- true_graph_neg[3, 1] <- 1

# generate the covariance matrices as a function of Z
sigma_mats_neg <- lapply(Z_neg, function(z) z * true_graph_neg + diag(p + 1))
sigma_mats_pos <- lapply(Z_pos, function(z) z * true_graph_pos + diag(p + 1))
sigma_mats <- c(sigma_mats_neg, sigma_mats_pos)

# generate the data using the covariance matrices
data_mat <- t(sapply(sigma_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))

# use varbvs to get the hyperparameter sigma
sigmasq <- rep(NA, p + 1)
for (j in 1:(p + 1)){
  sigmasq[j] <- mean(varbvs::varbvs(data_mat[ , -j], Z, data_mat[ , j], verbose = F)$sigma)
}

# estimate the covariance structure
out <- covdepGE(data_mat,
                Z, # extraneous covariates
                sigmasq = mean(sigmasq), # hyperparameter residual variance
                var_min = 1e-16, # smallest sigmabeta_sq grid value
                var_max = 0.5, # largest sigmabeta_sq grid value
                n_sigma = 50, # length of the sigmabeta_sq grid
                pi_vec = 0.1, # prior inclusion probability grid
                tolerance = 1e-5, # variational parameter exit condition 1
                max_iter = 1e3, # variational parameter exit condition 2
                edge_threshold = 0.75, # minimum inclusion probability
                sym_method = "min", # how to symmetrize the alpha matrices
                print_time = T,
)
out$ELBO

# BUG: constant covariate

# create covariate for individuals in each of the three clusters

# define the dimensions of the data
n1 <- 60; n2 <- 60; n3 <- 60
n <- sum(n1, n2, n3)
p <- 4

# define the limits of the covariate
limits1 <- c(-.990, -.331)
limits2 <- c(-.229, 0.329)
limits3 <- c(0.431, 0.990)

# define the covariate
z1 <- seq(limits1[1], limits1[2], length = n %/% 3)
z2 <- seq(limits2[1], limits2[2], length = n %/% 3)
z3 <- seq(limits3[1], limits3[2], length = n %/% 3)
Z <- matrix(c(z1, z2, z3), n, 1)

# create the precision matrices for each individual

# the shared part of the structure for all three clusters is a 2 on the diagonal and a 1 in the (2, 3) position
common_str <- diag(p + 1)
common_str[2, 3] <- 1

# define constants for the structure of cluster 2
const1 <- 0.23; const2 <- 0.56

# cluster 2 has two different linear functions of Z in the (1, 2) and (1, 3) positions; define structures for each of these components
cl2_str12 <- cl2_str13 <- matrix(0, p + 1, p + 1)
cl2_str12[1, 2] <- cl2_str13[1, 3] <- 1

# define the precision matrices for each of the individuals in cluster 2
cl2_prec <- lapply(z2, function(z) common_str + ((1 - (z + const1) / const2) * cl2_str12) + ((z + const1) / const2 * cl2_str13))

# cluster 1 has a 1 in the (1, 2) and cluster 3 has a 1 in the (1, 3) position; define structures for each of these components
cl1_str12 <- cl3_str13 <- matrix(0, p + 1, p + 1)
cl1_str12[1, 2] <- cl3_str13[1, 3] <- 1

# define the precision matrices for each of the individuals in cluster 1 and cluster 3
cl1_prec <- rep(list(common_str + cl1_str12), n1)
cl3_prec <- rep(list(common_str + cl3_str13), n3)

# put all of the precision matrices into one list
prec_mats <- c(cl1_prec, cl2_prec, cl3_prec)

# symmetrize them and invert to get the covariance matrices
prec_mats <- lapply(prec_mats, function(mat) t(mat) + mat)
cov_mats <- lapply(prec_mats, solve)

# generate the data using the covariance matrices
data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))

# use varbvs to get the hyperparameter sigma
sigmasq <- rep(NA, p + 1)
for (j in 1:(p + 1)){
  sigmasq[j] <- mean(varbvs::varbvs(data_mat[ , -j], Z, data_mat[ , j], verbose = F)$sigma)
}

# estimate the dependence structure independent of the covariate

out_indep <- covdepGE(data_mat,
                      rep(0, n), # extraneous covariates
                      sigmasq = mean(sigmasq), # hyperparameter residual variance
                      var_min = 1e-13, # smallest sigmabeta_sq grid value
                      var_max = 0.4, # largest sigmabeta_sq grid value
                      n_sigma = 100, # length of the sigmabeta_sq grid
                      pi_vec = 0.1, # prior inclusion probability grid
                      tolerance = 1e-12, # variational parameter exit condition 1
                      max_iter = 1e3, # variational parameter exit condition 2
                      edge_threshold = 0.75, # minimum inclusion probability
                      sym_method = "min", # how to symmetrize the alpha matrices
                      print_time = T,
)

out_indep <- covdepGE(data_mat,
                      rep(0, n), # extraneous covariates
                      sigmasq = mean(sigmasq), # hyperparameter residual variance
                      var_min = 1e-13, # smallest sigmabeta_sq grid value
                      var_max = 0.4, # largest sigmabeta_sq grid value
                      n_sigma = 100, # length of the sigmabeta_sq grid
                      pi_vec = 0.1, # prior inclusion probability grid
                      tolerance = 1e-12, # variational parameter exit condition 1
                      max_iter = 1e3, # variational parameter exit condition 2
                      edge_threshold = 0.75, # minimum inclusion probability
                      sym_method = "min", # how to symmetrize the alpha matrices
                      print_time = T,
                      scale = F
)

out_indep <- covdepGE(data_mat,
                      rep(0, n), # extraneous covariates
                      sigmasq = mean(sigmasq), # hyperparameter residual variance
                      var_min = 1e-13, # smallest sigmabeta_sq grid value
                      var_max = 0.4, # largest sigmabeta_sq grid value
                      n_sigma = 100, # length of the sigmabeta_sq grid
                      pi_vec = 0.1, # prior inclusion probability grid
                      tolerance = 1e-12, # variational parameter exit condition 1
                      max_iter = 1e3, # variational parameter exit condition 2
                      edge_threshold = 0.75, # minimum inclusion probability
                      sym_method = "min", # how to symmetrize the alpha matrices
                      print_time = T,
                      kde = F
)

out_indep <- covdepGE(data_mat,
                      rep(1, n), # extraneous covariates
                      sigmasq = mean(sigmasq), # hyperparameter residual variance
                      var_min = 1e-15,
                      var_max = 5,
                      n_sigma = 50,
                      pi_vec = 0.1, # prior inclusion probability grid
                      tolerance = 1e-10, # variational parameter exit condition 1
                      max_iter = 1e3, # variational parameter exit condition 2
                      edge_threshold = 0.75, # minimum inclusion probability
                      sym_method = "min", # how to symmetrize the alpha matrices
                      print_time = T,
                      scale = F,
                      kde = F
)

out_indep$ELBO
unique(out_indep$alpha_matrices)

## -----------------------------------------------------------------------------
## -----------------------------covdepGE----------------------------------------
## -----------------------------------------------------------------------------

## -----------------------------data_mat----------------------------------------

# NA
data_mat_NA <- data_mat
data_mat_NA[1 , 1] <- NA
covdepGE(data_mat_NA, Z)

# Inf
data_mat_Inf <- data_mat
data_mat_Inf[1 , 1] <- Inf
covdepGE(data_mat_Inf, Z)

# wrong input
covdepGE(lm, Z)

## -----------------------------Z-----------------------------------------------

# NA
Z_NA <- Z
Z_NA[1] <- NA
covdepGE(data_mat, Z_NA)

# Inf
Z_Inf <- Z
Z_Inf[1] <- Inf
covdepGE(data_mat, Z_Inf)

# wrong input
covdepGE(data_mat, lm)

# wrong size
covdepGE(data_mat, Z[-1])

## -----------------------------tau---------------------------------------------

# NA
covdepGE(data_mat, Z, tau = NA)

# Inf
covdepGE(data_mat, Z, tau = Inf)

# wrong input
covdepGE(data_mat, Z, tau = "here")

# non-positive
covdepGE(data_mat, Z, tau = 0)

# vector with one non-positive
covdepGE(data_mat, Z, tau = c(0, rep(1, nrow(data_mat))))

# vector of the wrong length
covdepGE(data_mat, Z, tau = rep(1, nrow(data_mat) - 1))

# matrix
covdepGE(data_mat, Z, tau = matrix(rep(1, nrow(data_mat), 9, 20)))

## -----------------------------kde---------------------------------------------

# NA
covdepGE(data_mat, Z, kde = NA)

# non-logical
covdepGE(data_mat, Z, kde = 4)

# bool vect
covdepGE(data_mat, Z, kde = c(T, F))

# bool matrix
covdepGE(data_mat, Z, kde = matrix(c(T, F)))

## -----------------------------alpha-------------------------------------------

# NA
covdepGE(data_mat, Z, alpha = NA)

# Inf
covdepGE(data_mat, Z, alpha = Inf)

# outside of [0,1]
covdepGE(data_mat, Z, alpha = 4)

# non-scalar
covdepGE(data_mat, Z, alpha = c(0.5, 0.5))

# matrix
covdepGE(data_mat, Z, alpha = matrix(rep(0.5, 4), 2))

## -----------------------------mu----------------------------------------------

# NA
covdepGE(data_mat, Z, mu = NA)

# Inf
covdepGE(data_mat, Z, mu = Inf)

# non-scalar
covdepGE(data_mat, Z, mu = c(0.5, 0.5))

# matrix
covdepGE(data_mat, Z, mu = matrix(rep(0.5, 4), 2))

## -----------------------------sigma_sq----------------------------------------

# NA
covdepGE(data_mat, Z, sigmasq = NA)

# Inf
covdepGE(data_mat, Z, sigmasq = Inf)

# non-scalar
covdepGE(data_mat, Z, sigmasq = c(0.5, 0.5))

# matrix
covdepGE(data_mat, Z, sigmasq = matrix(rep(0.5, 4), 2))

# negative
covdepGE(data_mat, Z, sigmasq = -5)

## -----------------------------sigmabetasq_vec---------------------------------

# with an NA
covdepGE(data_mat, Z, sigmabetasq_vec = c(NA, 2))

# with an Inf
covdepGE(data_mat, Z, sigmabetasq_vec = c(Inf, 2))

# with a negative
covdepGE(data_mat, Z, sigmabetasq_vec = c(-0.5, 0.5))

# non-numeric
covdepGE(data_mat, Z, sigmabetasq_vec = c(T, 0.5))

# non-numeric scalar
covdepGE(data_mat, Z, sigmabetasq_vec = F)

# matrix
covdepGE(data_mat, Z, sigmabetasq_vec = matrix(rep(4, 4), 2))

## -----------------------------var_min-----------------------------------------

# NA
covdepGE(data_mat, Z, var_min = NA)

# Inf
covdepGE(data_mat, Z, var_min = Inf)

# negative
covdepGE(data_mat, Z, var_min = -0.5)

# non-numeric
covdepGE(data_mat, Z, var_min = F)

# non-scalar
covdepGE(data_mat, Z, var_min = c(1, 2))

# matrix
covdepGE(data_mat, Z, var_min = matrix(rep(4, 4), 2))

## -----------------------------var_max-----------------------------------------

# NA
covdepGE(data_mat, Z, var_max = NA)

# Inf
covdepGE(data_mat, Z, var_max = Inf)

# negative
covdepGE(data_mat, Z, var_max = -0.5)

# non-numeric
covdepGE(data_mat, Z, var_max = F)

# non-scalar
covdepGE(data_mat, Z, var_max = c(1, 2))

# matrix
covdepGE(data_mat, Z, var_max = matrix(rep(4, 4), 2))

# supply var_max < var_min
covdepGE(data_mat, Z, var_min = 1, var_max = 0.5)

## -----------------------------n_sigma-----------------------------------------

# NA
covdepGE(data_mat, Z, n_sigma = NA)

# Inf
covdepGE(data_mat, Z, n_sigma = Inf)

# negative
covdepGE(data_mat, Z, n_sigma = -1)

# non-integer
covdepGE(data_mat, Z, n_sigma = 1.5)

# non-scalar
covdepGE(data_mat, Z, n_sigma = c(1, 2))

# matrix
covdepGE(data_mat, Z, n_sigma = matrix(rep(4, 4), 2))

## -----------------------------pi_vec------------------------------------------

# with an NA
covdepGE(data_mat, Z, pi_vec = c(NA, 2))

# with an Inf
covdepGE(data_mat, Z, pi_vec = c(Inf, 2))

# with a negative
covdepGE(data_mat, Z, pi_vec = c(-0.5, 0.5))

# non-numeric
covdepGE(data_mat, Z, pi_vec = c(F, 0.5))

# non-numeric scalar
covdepGE(data_mat, Z, pi_vec = F)

# matrix
covdepGE(data_mat, Z, pi_vec = matrix(rep(4, 4), 2))

## -----------------------------norm--------------------------------------------

# NA
covdepGE(data_mat, Z, norm = NA)

# less than 1
covdepGE(data_mat, Z, norm = 0.5)

# non-numeric
covdepGE(data_mat, Z, norm = T)

# non-scalar
covdepGE(data_mat, Z, norm = c(1, 2))

# matrix
covdepGE(data_mat, Z, norm = matrix(rep(4, 4), 2))

## -----------------------------scale-------------------------------------------

# NA
covdepGE(data_mat, Z, scale = NA)

# non-logical
covdepGE(data_mat, Z, scale = 10)

# logicial vector
covdepGE(data_mat, Z, scale = c(T, F))

## -----------------------------tolerance---------------------------------------

# NA
covdepGE(data_mat, Z, tolerance = NA)

# Inf
covdepGE(data_mat, Z, tolerance = Inf)

# negative
covdepGE(data_mat, Z, tolerance = -1)

# non-scalar
covdepGE(data_mat, Z, tolerance = c(1, 2))

## -----------------------------max_iter----------------------------------------

# NA
covdepGE(data_mat, Z, max_iter = NA)

# Inf
covdepGE(data_mat, Z, max_iter = Inf)

# negative
covdepGE(data_mat, Z, max_iter = -1)

# non-integer
covdepGE(data_mat, Z, max_iter = 1.5)

# non-scalar
covdepGE(data_mat, Z, max_iter = c(1, 2))

# matrix
covdepGE(data_mat, Z, max_iter = matrix(rep(4, 4), 2))

## -----------------------------edge_threshold----------------------------------

# NA
covdepGE(data_mat, Z, edge_threshold = NA)

# Inf
covdepGE(data_mat, Z, edge_threshold = Inf)

# outside of (0,1)
covdepGE(data_mat, Z, edge_threshold = 1)

# non-scalar
covdepGE(data_mat, Z, edge_threshold = c(0.5, 0.5))

# matrix
covdepGE(data_mat, Z, edge_threshold = matrix(rep(0.5, 4), 2))

## -----------------------------sym_method--------------------------------------

# NA
covdepGE(data_mat, Z, sym_method = NA)

# Inf
covdepGE(data_mat, Z, sym_method = Inf)

# not valid character
covdepGE(data_mat, Z, sym_method = "Inf")

# non-scalar
covdepGE(data_mat, Z, sym_method = c("max", "min"))

## -----------------------------print_time--------------------------------------

# NA
covdepGE(data_mat, Z, print_time = NA)

# non-logical
covdepGE(data_mat, Z, print_time = 10)

# logicial vector
covdepGE(data_mat, Z, print_time = c(T, F))

## -----------------------------warnings----------------------------------------

# NA
covdepGE(data_mat, Z, warnings = NA)

# non-logical
covdepGE(data_mat, Z, warnings = 10)

# logicial vector
covdepGE(data_mat, Z, warnings = c(T, F))

## -----------------------------------------------------------------------------
## -----------------------------gg_adjMat---------------------------------------
## -----------------------------------------------------------------------------

out <- covdepGE(data_mat, Z)
gg_adjMat(out, 1)

## -----------------------------out---------------------------------------------

# non-list/ matrix
gg_adjMat(7, 1)

# list without proper values
gg_adjMat(list(7), 1)

# matrix with non-numerics
gg_adjMat(matrix(sample(letters, 25), 5))

# matrix with more than two values
gg_adjMat(matrix(rnorm(100), 10), color1 = "dodgerblue", grid_color = "tomato")

## -----------------------------l-----------------------------------------------

# NA
gg_adjMat(out, NA)

# Inf
gg_adjMat(out, Inf)

# Negative
gg_adjMat(out, -5)

# non-integer
gg_adjMat(out, 5.1)

# vector
gg_adjMat(out, c(5, 5))

## -----------------------------prob_shade--------------------------------------

# NA
gg_adjMat(out, 1, prob_shade = NA)

# Inf
gg_adjMat(out, 1, prob_shade = Inf)

# non-logical
gg_adjMat(out, 1, prob_shade = 5)

# logical vector
gg_adjMat(out, 1, prob_shade = c(T, T))

## -----------------------------color0------------------------------------------

# NA
gg_adjMat(out, 1, color0 = NA)

# Inf
gg_adjMat(out, 1, color0 = Inf)

# non-color character
gg_adjMat(out, 1, color0 = "here")

# vector of colors
gg_adjMat(out, 1, color0 = c("red", "blue"))

## -----------------------------color1------------------------------------------

# NA
gg_adjMat(out, 1, color1 = NA)

# Inf
gg_adjMat(out, 1, color1 = Inf)

# non-color character
gg_adjMat(out, 1, color1 = "here")

# vector of colors
gg_adjMat(out, 1, color1 = c("red", "blue"))

## -----------------------------grid_color------------------------------------------

# NA
gg_adjMat(out, 1, grid_color = NA)

# Inf
gg_adjMat(out, 1, grid_color = Inf)

# non-color character
gg_adjMat(out, 1, grid_color = "here")

# vector of colors
gg_adjMat(out, 1, grid_color = c("red", "blue"))

## -----------------------------incl_probs--------------------------------------

# NA
gg_adjMat(out, 1, incl_probs = NA)

# Inf
gg_adjMat(out, 1, incl_probs = Inf)

# non-logical
gg_adjMat(out, 1, incl_probs = 5)

# logical vector
gg_adjMat(out, 1, incl_probs = c(T, T))

## -----------------------------prob_prec--------------------------------------

# NA
gg_adjMat(out, 1, prob_prec = NA)

# Inf
gg_adjMat(out, 1, prob_prec = Inf)

# non-integer
gg_adjMat(out, 1, prob_prec = 5.1)

# integer vector
gg_adjMat(out, 1, prob_prec = c(1, 1))

# negative
gg_adjMat(out, 1, prob_prec = -1)

## -----------------------------font_size---------------------------------------

# NA
gg_adjMat(out, 1, font_size = NA)

# Inf
gg_adjMat(out, 1, font_size = Inf)

# vector
gg_adjMat(out, 1, font_size = c(1, 1))

# negative
gg_adjMat(out, 1, font_size = -1)

# non-numeric
gg_adjMat(out, 1, font_size = "1")

## -----------------------------font_color0-------------------------------------

# NA
gg_adjMat(out, 1, font_color0 = NA)

# Inf
gg_adjMat(out, 1, font_color0 = Inf)

# non-color character
gg_adjMat(out, 1, font_color0 = "here")

# vector of colors
gg_adjMat(out, 1, font_color0 = c("red", "blue"))

## -----------------------------font_color1-------------------------------------

# NA
gg_adjMat(out, 1, font_color1 = NA)

# Inf
gg_adjMat(out, 1, font_color1 = Inf)

# non-color character
gg_adjMat(out, 1, font_color1 = "here")

# vector of colors
gg_adjMat(out, 1, font_color1 = c("red", "blue"))

## -----------------------------------------------------------------------------
## -----------------------------gg_inclusionCurve-------------------------------
## -----------------------------------------------------------------------------

out <- covdepGE(data_mat, Z)
gg_inclusionCurve(out, 1, 2)

## -----------------------------out---------------------------------------------

# non-list
gg_inclusionCurve(7, 1, 2)

# list without proper values
gg_inclusionCurve(list(7), 1, 2)

## -----------------------------col_idx1----------------------------------------

# NA
gg_inclusionCurve(out, NA, 2)

# Inf
gg_inclusionCurve(out, Inf, 2)

# Negative
gg_inclusionCurve(out, -5, 2)

# non-integer
gg_inclusionCurve(out, 5.1, 2)

# vector
gg_inclusionCurve(out, c(5, 5), 2)

## -----------------------------col_idx1----------------------------------------

# NA
gg_inclusionCurve(out, 2, NA)

# Inf
gg_inclusionCurve(out, 2, Inf)

# Negative
gg_inclusionCurve(out, 5, -2)

# non-integer
gg_inclusionCurve(out, 5, 2.1)

# vector
gg_inclusionCurve(out, 2, c(5, 5))

## -----------------------------line_type---------------------------------------

# NA
gg_inclusionCurve(out, 2, 3, line_type = NA)

# Inf
gg_inclusionCurve(out, 2, 3, line_type = Inf)

# non-character
gg_inclusionCurve(out, 2, 5, line_type = -2)

# wrong character
gg_inclusionCurve(out, 5, 2, line_type = "here")

# vector
gg_inclusionCurve(out, 5, 2, line_type = c(1, 2))

## -----------------------------line_size---------------------------------------

# NA
gg_inclusionCurve(out, 1, 2, line_size = NA)

# Inf
gg_inclusionCurve(out, 1, 2, line_size = Inf)

# non-numeric
gg_inclusionCurve(out, 1, 2, line_size = "5")

# vector
gg_inclusionCurve(out, 1, 2, line_size = c(1, 2))

## -----------------------------line_color--------------------------------------

# NA
gg_inclusionCurve(out, 1, 2, line_color = NA)

# Inf
gg_inclusionCurve(out, 1, 2, line_color = Inf)

# non-color character
gg_inclusionCurve(out, 1, 2, line_color = "here")

# vector of colors
gg_inclusionCurve(out, 1, 2, line_color = c("red", "blue"))

## -----------------------------point_shape-------------------------------------

# NA
gg_inclusionCurve(out, 1, 2, point_shape = NA)

# Inf
gg_inclusionCurve(out, 1, 2, point_shape = Inf)

## -----------------------------point_size--------------------------------------

# NA
gg_inclusionCurve(out, 1, 2, point_size = NA)

# Inf
gg_inclusionCurve(out, 1, 2, point_size = Inf)

# non-numeric
gg_inclusionCurve(out, 1, 2, point_size = "here")

# vector
gg_inclusionCurve(out, 1, 2, point_size = c(1, 2))

## -----------------------------point_color-------------------------------------

# NA
gg_inclusionCurve(out, 1, 2, point_color = NA)

# Inf
gg_inclusionCurve(out, 1, 2, point_color = Inf)

# non-color
gg_inclusionCurve(out, 1, 2, point_color = "here")

# vector of colors
gg_inclusionCurve(out, 1, 2, point_color = c("red", "blue"))

## -----------------------------point_fill-------------------------------------

# NA
gg_inclusionCurve(out, 1, 2, point_fill = NA)

# Inf
gg_inclusionCurve(out, 1, 2, point_fill = Inf)

# non-color
gg_inclusionCurve(out, 1, 2, point_fill = "here")

# vector of colors
gg_inclusionCurve(out, 1, 2, point_fill = c("red", "blue"))

## -----------------------------sort--------------------------------------------

# NA
gg_inclusionCurve(out, 1, 2, sort = NA)

# Inf
gg_inclusionCurve(out, 1, 2, sort = Inf)

# non-logical
gg_inclusionCurve(out, 1, 2, sort = 5.1)

# logical vector
gg_inclusionCurve(out, 1, 2, sort = c(T, T))

## -----------------------------------------------------------------------------
## -----------------------------gg_adjMats--------------------------------------
## -----------------------------------------------------------------------------

out <- covdepGE(data_mat, Z)
gg_adjMats(out)

## -----------------------------out---------------------------------------------

# non-list
gg_adjMats(7)

# list without proper values
gg_adjMats(list(7), 1, 2)

## -----------------------------graph_colors------------------------------------

# NA
gg_adjMats(out, rep(NA, 4))

# Inf
gg_adjMats(out, rep(Inf, 4))

# non-color
gg_adjMats(out, rep("here", 4))

# too short vector of colors
gg_adjMats(out, c("red", "blue"))


## -----------------------------seed--------------------------------------------

# NA
gg_adjMats(out, seed = NA)

# Inf
gg_adjMats(out, seed = Inf)

# vector
gg_adjMats(out, seed = c(5, 5))

# non-numeric
gg_adjMats(out, seed = matrix(c(1, "here")))


