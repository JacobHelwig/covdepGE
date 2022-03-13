// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cavi_c
Rcpp::List cavi_c(const arma::colvec& y, const arma::mat& D, const arma::mat& X, const arma::mat& mu0, const arma::mat& alpha0, double ssq, double sbsq, double pip, double elbo_tol, double alpha_tol, int max_iter, bool grid_search);
RcppExport SEXP _covdepGE_cavi_c(SEXP ySEXP, SEXP DSEXP, SEXP XSEXP, SEXP mu0SEXP, SEXP alpha0SEXP, SEXP ssqSEXP, SEXP sbsqSEXP, SEXP pipSEXP, SEXP elbo_tolSEXP, SEXP alpha_tolSEXP, SEXP max_iterSEXP, SEXP grid_searchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< double >::type ssq(ssqSEXP);
    Rcpp::traits::input_parameter< double >::type sbsq(sbsqSEXP);
    Rcpp::traits::input_parameter< double >::type pip(pipSEXP);
    Rcpp::traits::input_parameter< double >::type elbo_tol(elbo_tolSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_tol(alpha_tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type grid_search(grid_searchSEXP);
    rcpp_result_gen = Rcpp::wrap(cavi_c(y, D, X, mu0, alpha0, ssq, sbsq, pip, elbo_tol, alpha_tol, max_iter, grid_search));
    return rcpp_result_gen;
END_RCPP
}
// grid_search_c
Rcpp::List grid_search_c(const arma::colvec& y, const arma::mat& D, const arma::mat& X, const arma::mat& mu, const arma::mat& alpha, const arma::colvec& ssq, const arma::colvec& sbsq, const arma::colvec& pip, double elbo_tol, double alpha_tol, int max_iter, bool grid_search);
RcppExport SEXP _covdepGE_grid_search_c(SEXP ySEXP, SEXP DSEXP, SEXP XSEXP, SEXP muSEXP, SEXP alphaSEXP, SEXP ssqSEXP, SEXP sbsqSEXP, SEXP pipSEXP, SEXP elbo_tolSEXP, SEXP alpha_tolSEXP, SEXP max_iterSEXP, SEXP grid_searchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type ssq(ssqSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type sbsq(sbsqSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type pip(pipSEXP);
    Rcpp::traits::input_parameter< double >::type elbo_tol(elbo_tolSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_tol(alpha_tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< bool >::type grid_search(grid_searchSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_search_c(y, D, X, mu, alpha, ssq, sbsq, pip, elbo_tol, alpha_tol, max_iter, grid_search));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_covdepGE_cavi_c", (DL_FUNC) &_covdepGE_cavi_c, 12},
    {"_covdepGE_grid_search_c", (DL_FUNC) &_covdepGE_grid_search_c, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_covdepGE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
