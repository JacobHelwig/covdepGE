#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List dev(){
  return(Rcpp::List::create(Rcpp::Named("test") = 1));
}

// [[Rcpp::export]]
double ELBO_calculator_c (const arma::colvec& y, const arma::colvec& D,
                          const arma::mat& X, const arma::colvec& ssq_var,
                          const arma::colvec& mu, const arma::colvec& alpha,
                          double ssq, double sbsq, double pip) {

  // square of mu vector
  arma::colvec mu_sq = arma::pow(mu, 2);

  // calculate the expected value of the log-prior under q
  double eqpr = sum(-alpha / 2 * log(ssq * sbsq) - alpha % (ssq_var + mu_sq) /
                    (2 * ssq * sbsq) + alpha * log(pip) + (1 - alpha) * log(1 - pip));

  // calculate the expected value of the log-likelihood under q
  double eqlik = -sum(
    D % (arma::pow(y - (X * (mu % alpha)), 2) + arma::pow(X, 2) *
      (alpha % (mu_sq + ssq_var) - arma::pow(alpha, 2) % mu_sq))) /
      (2 * ssq) - X.n_rows / 2 * log(2 * M_PI * ssq);

  // calculate the expected value of log q under q
  double eqq = sum(-alpha / 2 % log(ssq_var) - alpha / 2 + alpha % log(alpha +
                   0.000001) + (1 - alpha) % log(1 - alpha + 0.000001));

  // calculate elbo and return
  return(eqpr + eqlik - eqq);
}



/*** R
# ELBO dev
#rm(list = ls())
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/covdepGE_R.R")
source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/weights.R")
source("generate_data.R")
set.seed(1)

dat <- generate_continuous()
X <- dat$data
Z <- dat$covts
n <- nrow(X)
p <- ncol(X)

y <- X[ , 1]
D <- get_weights(Z, 2, T, 0)$D[1 , ]
X <- X[ , -1]
ssq_var <- rexp(p - 1)
mu <- rnorm(p - 1)
alpha <- runif(p - 1)
ssq <- 0.6
sbsq <- 0.5
pip <- 0.2

all.equal(ELBO_calculator_R(y, D, X, ssq_var, mu, alpha, ssq, sbsq, pip), ELBO_calculator_c(y, D, X, ssq_var, mu, alpha, ssq, sbsq, pip))
*/
