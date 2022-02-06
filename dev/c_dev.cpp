#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
arma::mat S_sq_update(const arma::mat& D, const arma::mat& X_mat,
                      arma::colvec& sigmasq, arma::colvec& sigmabeta_sq) {

  int p = X_mat.n_cols;

  arma::mat S_sq = arma::repmat(sigmasq, 1, p) /
    ((arma::pow(X_mat, 2).t() * D).t() + 1 / arma::repmat(sigmabeta_sq, 1, p));

  return S_sq;
}


/***R
set.seed(1)
n <- 10; p <- 3
D <- matrix(abs(rnorm(n^2)), n, n)
X_mat <- matrix(rnorm(n * p), n, p)
sigmasq <- rnorm(n)
sigmabeta_sq <- rnorm(n)

all.equal(S_sq_update(D, X_mat, sigmasq, sigmabeta_sq), matrix(sigmasq, n, p) / (t(t(X_mat^2) %*% D) + 1 / matrix(sigmabeta_sq, n, p)))

*/
