#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
Rcpp::List test(arma::mat X){
  arma::mat Y = X;
  X = X + 1;
  Y = arma::log(Y);
  Y.replace(-arma::datum::inf, log(DBL_MIN));
  return(Rcpp::List::create(Rcpp::Named("X") = X, Rcpp::Named("Y") = Y));
}


/***R
X <- matrix(c(0, 1e-300))
test(X)
*/
