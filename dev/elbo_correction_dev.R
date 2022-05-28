# ELBO dev
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/covdepGE_R.R")
Rcpp::sourceCpp("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/src/covdepGE_c.cpp")
Rcpp::sourceCpp()
rm(list = ls())
source("generate_data.R")
set.seed(1)

dat <- generate_continuous()
X <- dat$data
Z <- dat$covts

