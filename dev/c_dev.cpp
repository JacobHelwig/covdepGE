#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
Rcpp::List test(){

  Rcpp::List res1(200);
  Rcpp::List res2(200);
  Rcpp::List res3(200);

  for (int i = 0; i<200; i++){
    res1(i) = arma::randn(200, 5);
    res2(i) = arma::randn(200, 5);
    res3(i) = arma::randn(200, 5);
  }

  return(Rcpp::List::create(
      Rcpp::Named("mu") = res1, Rcpp::Named("alpha") = res2,
      Rcpp::Named("ssq_var") = res3));
}

//[[Rcpp::export]]
Rcpp::List test2(){

  arma::cube res(600, 5, 200);

  for (int i = 0; i<200; i++){
    res.slice(i) = arma::randn(600, 5);
  }

  return(Rcpp::List::create(Rcpp::Named("mu") = res));
}

//[[Rcpp::export]]
Rcpp::List test3(){

  arma::mat res1 = arma::randn(200, 5);
  arma::mat res2 = arma::randn(200, 5);
  arma::mat res3 = arma::randn(200, 5);

  return(Rcpp::List::create(
      Rcpp::Named("mu") = res1, Rcpp::Named("alpha") = res2,
      Rcpp::Named("ssq_var") = res3));
}

//[[Rcpp::export]]
Rcpp::List test4(){

  Rcpp::List res1(200);
  Rcpp::List res2(200);
  Rcpp::List res3(200);

  for (int i = 0; i<200; i++){
    res1(i) = arma::randn(200, 5);
    res2(i) = arma::randn(200, 5);
    res3(i) = arma::randn(200, 5);
  }

  return(Rcpp::List::create(
      Rcpp::Named("mu") = res1, Rcpp::Named("alpha") = res2,
      Rcpp::Named("ssq_var") = res3));
}


/***R
library(microbenchmark)
test3_R <- function(){
  res <- vector("list", 200)
  for (j in 1:200){
    res[[j]] <- test3()
  }
  return(res)
}
microbenchmark(test(), test2(), test3_R())
*/
