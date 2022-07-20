setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev")
rm(list = ls())
source("generate_data.R")

discrete_data <- F
package <- F
sbq <- c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10)

# generate data and covariates
if (discrete_data) {
  dat <- generate_discrete()
  tau_ <- 0.1 # the bandwidth parameter
}else{
  dat <- generate_continuous()
  tau_ <- 0.56
}

data_mat <- dat$data
Z <- dat$covts

# fit varbvs to get pi and ssq
sigma <- pip_ <- rep(NA, ncol(data_mat))
for (j in 1:ncol(data_mat)){
  vout <- varbvs::varbvs(data_mat[ , -j], NULL, data_mat[ , j], verbose = F)
  sigma[j] <- mean(vout$sigma)
  pip_[j] <- mean(1/(1 + 10^(-vout$logodds)))
}
pip_ <- unique(pip_)

if (package){
  out <- covdepGE::covdepGE(data_mat, Z, tau = tau_, kde = F, CS = T, scale = F,
                            sbsq = c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10),
                            R = R_code, max_iter = 100, warnings = F,
                            alpha_tol = 1e-10, center_data = F,
                            parallel = parallel, num_workers = ncol(data_mat))

}else{
  if ("covdepGE" %in% .packages()) detach("package:covdepGE", unload = TRUE)
  source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/main.R")
  source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/cavi.R")
  source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/weights.R")
  source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/plots.R")
  source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/data.R")
  Rcpp::sourceCpp("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/src/covdepGE_c.cpp")
  out <- covdepGE(data_mat, Z, "grid_search", ssq = sigma, sbsq = sbq,
                  pip = pip_, tau = tau_)
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
for (j in 1:length(out$graphs$inclusion_probs_asym)) {
  if (all.equal(out$graphs$inclusion_probs_asym[[j]],
                out_original$original_alpha_matrices[[j]]) != T) {
    total_diff <- total_diff + norm(out$graphs$inclusion_probs_asym[[j]] -
                                      out_original$original_alpha_matrices[[j]], "F")
    same_alpha <- F
  }
}
total_diff
same_alpha

summary(abs(unlist(out_original$original_alpha_matrices) - unlist(out$graphs$inclusion_probs_sym)))
