setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev")
rm(list = ls())
source("generate_data.R")

R_code <- F # true if R code instead of C++ should be used

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
Z <- dat$covts

package <- F # true if the package version is desired
if (package){
  out <- covdepGE::covdepGE(data_mat, Z, tau_, kde = F, CS = T, scale = F,
                            sigmabetasq_vec = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10),
                            parallel = T, num_workers = 5)
}else{
  if ("covdepGE" %in% .packages()) detach("package:covdepGE", unload = TRUE)
  source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/covdepGE_main.R")
  source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/cavi_search.R")
  source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/weights.R")
  source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/checks.R")
  source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/gg_covdepGE.R")
  source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/covdepGE_R.R")
  Rcpp::sourceCpp("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/src/covdepGE_c.cpp")
  out <- covdepGE(data_mat, Z, tau_, kde = F, CS = T, scale = F,
                  sigmabetasq_vec = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10), R = R_code)
}

# check to see that this modified code produces the same results as the original code
if (discrete_data){
  load("out_original_discrete.Rdata")
}else{
  load("out_original_continuous.Rdata")
}

# check for equality between the alpha matrices
same_alpha <- T
total_diff <- 0
for (j in 1:length(out$alpha_matrices)) {
  if (all.equal(out$alpha_matrices[[j]],
                out_original$original_alpha_matrices[[j]]) != T) {
    total_diff <- total_diff + norm(out$alpha_matrices[[j]] -
                                      out_original$original_alpha_matrices[[j]], "F")
    same_alpha <- F
  }
}
total_diff
same_alpha

# check for equality between the inclusion probabilities
same_probs <- T
total_diff <- 0
for (j in 1:length(out$inclusion_probs)) {
  if (all.equal(out$inclusion_probs[[j]],
                out_original$original_incl_probs[[j]]) != T) {
    same_probs <- F
    total_diff <- total_diff + norm(out$inclusion_probs[[j]] -
                                      out_original$original_incl_probs[[j]], "F")
  }
}
total_diff
same_probs

# check for equality between ELBO
all.equal(unname(unlist(lapply(out$CAVI_details, `[[`, "ELBO"))), out_original$original_ELBO)

