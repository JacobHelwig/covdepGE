#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
void test(arma::vec y, arma::vec x) {

  double var_x = arma::var(x);

  if (arma::any(y > var_x)){
    Rcout << y.elem(arma::find(y > var_x)) << "\n\n" << var_x;
  }

}


/***R

set.seed(2)

x <- rnorm(10); y <- rnorm(10)

test(y, x)
var(x)
*/
