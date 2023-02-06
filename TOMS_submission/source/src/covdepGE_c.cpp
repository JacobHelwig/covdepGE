#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// -----------------------------------------------------------------------------
// -----------------------------ELBO_calculator_c-------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// Calculates ELBO for a fixed response j and observation l
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x 1 vector; weights (i-th entry is the weight of the i-th observation
// with respect to the l-th observation using the l-th observation's bandwidth)
// X: n x p matrix; data_mat with the j-th column removed
// ssq_var, mu, alpha: p x 1 vectors; variational parameters. the k-th entry is
// the k-th parameter for the l-th observation
// ssq, sbsq: double; spike and slab variance hyperparameters
// pip: double; spike and slab probability of inclusion
// -----------------------------RETURNS-----------------------------------------
// Returns: double; ELBO for the l-th observation and j-th column fixed as the
// response
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double ELBO_calculator_c (const arma::colvec& y, const arma::colvec& D,
                          const arma::mat& X, const arma::colvec& ssq_var,
                          const arma::colvec& mu, const arma::colvec& alpha,
                          double ssq, double sbsq, double pip) {

  // square of mu vector
  arma::colvec mu_sq = arma::pow(mu, 2);

  // log of alpha and 1 - alpha with infinities replaced
  arma::colvec log_alpha = log(alpha);
  log_alpha.replace(-arma::datum::inf, log(DBL_MIN));
  arma::colvec log_1minalpha = log(1 - alpha);
  log_1minalpha.replace(-arma::datum::inf, log(DBL_MIN));

  // calculate the expected value of the log-prior under q
  double eqpr = sum(-alpha / 2 * log(ssq * sbsq) - alpha % (ssq_var + mu_sq) /
                    (2 * ssq * sbsq) + alpha * log(pip) + (1 - alpha) * log(1 - pip));

  // calculate the expected value of the log-likelihood under q
  double eqlik = -sum(
    D % (arma::pow(y - (X * (mu % alpha)), 2) + arma::pow(X, 2) *
      (alpha % (mu_sq + ssq_var) - arma::pow(alpha, 2) % mu_sq))) /
      (2 * ssq) - X.n_rows / 2 * log(2 * M_PI * ssq);

  // calculate the expected value of log q under q
  double eqq = sum(-alpha / 2 % log(ssq_var) - alpha / 2 + alpha % log_alpha
                     + (1 - alpha) % log_1minalpha);

  // calculate elbo and return
  return(eqpr + eqlik - eqq);
}

// -----------------------------------------------------------------------------
// -----------------------------total_elbo_c------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// for a fixed response, calculate and sum the ELBO for all observations in the
// data
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x n matrix; weights (k,l entry is the weight of the k-th observation
// with respect to the l-th observation using the l-th observation's bandwidth)
// X: n x p matrix; data_mat with the j-th column removed
// ssq_var, mu, alpha: n x p matrices; variational parameters. the l, k
// entry is the k-th parameter for the l-th observation
// ssq, sbsq: doubles; spike and slab variance hyperparameters
// pip: double; spike and slab probability of inclusion
// -----------------------------RETURNS-----------------------------------------
// elbo_tot: double; ELBO for the regression weighted with respect to the l-th
// observation with the j-th column fixed as the response
// -----------------------------------------------------------------------------
double total_ELBO_c (const arma::colvec& y, const arma::mat& D,
                     const arma::mat& X, const arma::mat& ssq_var,
                     const arma::mat& mu, const arma::mat& alpha, double ssq,
                     double sbsq, double pip){

  // get sample size
  int n = X.n_rows;

  // instantiate variable to store the total ELBO
  double elbo_tot = 0;

  // loop over all observations
  for (int l = 0; l < n; l++){

    // calculate the ELBO for l-th observation and add it to the total ELBO
    elbo_tot += ELBO_calculator_c(y, D.col(l), X, ssq_var.row(l).t(),
                                  mu.row(l).t(), alpha.row(l).t(), ssq, sbsq,
                                  pip);
  }

  return elbo_tot;
}

// -----------------------------------------------------------------------------
// -----------------------------mu_update_c-------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// for a fixed response, update mu matrix
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x n matrix; weights (k,l entry is the weight of the k-th observation
// with respect to the l-th observation using the l-th observation's bandwidth)
// X: n x p matrix; data_mat with the j-th column removed
// ssq_var, mu, alpha: n x p matrices; variational parameters. the l, k entry is
// the k-th parameter for the l-th observation
// ssq: double; spike and slab variance hyperparameter
// -----------------------------------------------------------------------------
void mu_update_c (const arma::colvec& y, const arma::mat& D, const arma::mat& X,
                  const arma::mat& ssq_var, arma::mat& mu,
                  const arma::mat& alpha, double ssq){

  // get sample size and number of parameters
  int n = X.n_rows;
  int p = X.n_cols;

  // instantiate matrices for the update loop
  arma::mat mu_stack(n, p), alpha_stack(n, p), X_mu_alpha(n, p);
  arma::mat X_mu_alpha_k(n, p), y_k(n, p), d_x_y(n, p);

  // loop over the observations to update mu row by row
  for (int l = 0; l < n; l++){

    // l-th row of mu, alpha stacked n times
    mu_stack = arma::repmat(mu.row(l), n, 1);
    alpha_stack = arma::repmat(alpha.row(l), n ,1);

    // take the element-wise product of X, mu_stack, and alpha stack
    X_mu_alpha = X % mu_stack % alpha_stack;

    // the k-th column is the rowSums of X_mu_alpha minus the k-th column of
    // X_mu_alpha (accounts for m \neq k in summation)
    X_mu_alpha_k = arma::repmat(arma::sum(X_mu_alpha, 1), 1, p) - X_mu_alpha;

    // the k-th column is y minus the k-th column of X_mu_alpha_k
    y_k = arma::repmat(y, 1, p) - X_mu_alpha_k;

    // the k-th column is d_:,l * x_:,k * y_k_:,k
    d_x_y = arma::repmat(D.col(l), 1, p) % X % y_k;

    // the update of the l-th row of mu
    mu.row(l) = ssq_var.row(l) % sum(d_x_y, 0) / ssq;
  }
}

// -----------------------------------------------------------------------------
// -----------------------------alpha_update_c----------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// for a fixed response, update alpha matrix
// -----------------------------ARGUMENTS---------------------------------------
// ssq_var, mu, alpha: n x p matrices; variational parameters. the l, k entry is
// the k-th parameter for the l-th observation
// alpha1, alpha2_denom, alpha3: double (1), n x p matrices (2, 3): terms used
// to calculate the logit of alpha
// ssq, sbsq: doubles; spike and slab variance hyperparameters
// pip: double; spike and slab probability of inclusion
// -----------------------------------------------------------------------------
void alpha_update_c(const arma::mat& ssq_var, const arma::mat& mu,
                    arma::mat& alpha, double alpha1,
                    const arma::mat& alpha2_denom, const arma::mat& alpha3,
                    double ssq, double sbsq, double pip){

  // calculate the logit of alpha
  arma::mat alpha_logit = (alpha1 + (arma::pow(mu, 2) / alpha2_denom) + alpha3);

  // transform from logit to probabilities of inclusion; update alpha
  arma::mat exp_logit = arma::exp(alpha_logit);
  alpha = exp_logit / (1 + exp_logit);

  // handle nan due to inf/inf resulting from exponentiation of large values;
  // these probabilities are indescernible from 1
  alpha.replace(arma::datum::nan, 1);
}

// -----------------------------------------------------------------------------
// -----------------------------cavi_c------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// The core function that performs CAVI to calculate and return the final
// variational estimates of a single regression for each of the n observations
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x n matrix; weights (k,l entry is the weight of the k-th observation
// with respect to the l-th observation using the l-th observation's bandwidth)
// X: n x p matrix; data_mat with the j-th column removed
// mu0, alpha0: n x p matrices; initial value of variational parameters. the
// l, k entry is the k-th parameter for the l-th observation
// ssq, sbsq: doubles; spike and slab variance hyperparameters
// pip: double; spike and slab probability of inclusion
// alpha_tol: double; when the square root of the sum of squared changes in
// the elements of alpha are within alpha_tol, stop iterating
// max_iter: integer; maximum number of iterations
// -----------------------------RETURNS-----------------------------------------
// mu; n x p matrix; final mu values
// alpha; n x p matrix; final alpha values
// ssq_var; n x p matrix; final ssq_var values
// elbo: double; final value of ELBO summed across all observations
// iter: integer; number of iterations to converge
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List cavi_c(const arma::colvec& y, const arma::mat& D,
                  const arma::mat& X, const arma::mat& mu0,
                  const arma::mat& alpha0, double ssq, double sbsq, double pip,
                  double alpha_tol, int max_iter) {

  // get sample size and number of parameters
  int n = X.n_rows;
  int p = X.n_cols;

  // matrices and vectors for updated variational parameters and hyperparameters
  arma::mat mu = mu0;
  arma::mat alpha = alpha0;

  // matrices for tracking the convergence of alpha parameters
  arma::mat alpha_last(n, p), change_alpha(n, p);

  // variational variance of regression coefficients update
  arma::mat ssq_var = ssq / ((arma::pow(X, 2).t() * D).t() + 1 / sbsq);

  // 1st and 3rd term of the alpha update, denominator of the second term
  double alpha1 = log(pip / (1 - pip));
  arma::mat alpha3 = arma::log(arma::sqrt(ssq_var / (ssq * sbsq)));
  arma::mat alpha2_denom = 2 * ssq_var;

  // track the number of iterations to converge
  int conv_iter = max_iter;

  // CAVI loop (optimize variational parameters)
  for (int k = 0; k < max_iter; k++){

    // mu update
    mu_update_c(y, D, X, ssq_var, mu, alpha, ssq);

    // alpha update

    // save the last value of alpha and update it
    alpha_last = alpha;
    alpha_update_c(ssq_var, mu, alpha, alpha1, alpha2_denom, alpha3, ssq, sbsq,
                   pip);

    // calculate element-wise change in alpha
    change_alpha = alpha - alpha_last;

    // if the Frobenius norm of the change in alpha is within the tolerance,
    // end updates
    if (sqrt(arma::accu(arma::pow(change_alpha, 2))) < alpha_tol){
      conv_iter = k + 1;
      break;
    }
  }

  // calculate ELBO across n observations
  double ELBO = total_ELBO_c(y, D, X, ssq_var, mu, alpha, ssq, sbsq, pip);

  // return final mu matrix, alpha matrix, the final ELBO, the number of
  // iterations to converge, and the fitted variance hyperparameters
  return(Rcpp::List::create(
      Rcpp::Named("mu") = mu, Rcpp::Named("alpha") = alpha,
      Rcpp::Named("ssq_var") = ssq_var, Rcpp::Named("elbo") = ELBO,
      Rcpp::Named("iter") = conv_iter));
}

// -----------------------------------------------------------------------------
// -----------------------------grid_search_c-----------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// for a fixed response, run CAVI for each observation across a grid of
// hyperparameters
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x n matrix; weights (k,l entry is the weight of the k-th observation
// with respect to the l-th observation using the l-th observation's bandwidth)
// X: n x p matrix; data_mat with the j-th column removed
// mu, alpha: n x p; matrices of variational parameters. the l, k entry
// is the k-th parameter for the l-th observation
// ssq: nssq x 1 vector; error term variance candidates
// sbsq: nsbsq x 1 vector; slab variance candidates
// pip: npip x 1 vector; candidate spike and slab probabilities of inclusion
// alpha_tol: double; when the square root of the sum of squared changes in
// the elements of alpha are within alpha_tol, stop iterating
// max_iter: integer; maximum number of iterations
// -----------------------------RETURNS-----------------------------------------
// elbo: double; final value of ELBO summed across all observations that achieved
// the maximum accross the hyperparameter grid
// mu; n x p matrix; final mu values corresponding to the hyperparameters that
// maximized the ELBO
// alpha; n x p matrix; final alpha values corresponding to the hyperparameters
// that maximized the ELBO
// ssq_var; n x p matrix; final ssq_var values corresponding to the
// hyperparameters that maximized the ELBO
// ssq: double; value of ssq that maximized the ELBO
// sbsq: double; value of sbsq that maximized the ELBO
// pip: double; value of pip that maximized the ELBO
// elbo_vec: (nssq * nsbsq * npip) x 1 vector; ELBO for each point in the
// hyperparameter grid
// iter: (nssq * nsbsq * npip) x 1 vector; number of iterations to converge for
// each point in the hyperparameter grid
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List grid_search_c(const arma::colvec& y, const arma::mat& D,
                         const arma::mat& X, const arma::mat& mu,
                         const arma::mat& alpha, const arma::colvec& ssq,
                         const arma::colvec& sbsq, const arma::colvec& pip,
                         double alpha_tol, int max_iter){

  // get sample size and number of parameters
  int n = X.n_rows;
  int p = X.n_cols;

  // get number of grid points
  int n_param = pip.n_elem;

  // storage for the best ELBO, new elbo, and converged iter
  double elbo_best, elbo_new;
  int iter_new;

  // storage for the hyperparameters and variational parameters corresponding to
  // the current best ELBO
  arma::mat mu_best(n, p), alpha_best(n, p), ssqv_best(n, p);
  double ssq_best, sbsq_best, pip_best;

  // instantiate a list to store the result of cavi_c
  Rcpp::List out(5);

  // vector to store ELBOs and converged iter
  arma::colvec elbo_store(n_param), iter_store(n_param);

  // perform CAVI for each grid point
  for (int j = 0; j < n_param; j++){

    // run CAVI; store elbo and number of iterations to converge
    out = cavi_c(y, D, X, mu, alpha, ssq(j), sbsq(j), pip(j), alpha_tol,
                 max_iter);
    elbo_new = as<double>(out["elbo"]);
    iter_new = as<double>(out["iter"]);
    elbo_store(j) = elbo_new;
    iter_store(j) = iter_new;

    // if the new ELBO is greater than the current best, update the best ELBO and
    // parameters and whether or not the model converged
    if (elbo_best <  elbo_new or j == 0){

      elbo_best = as<double>(out["elbo"]);
      mu_best = as<arma::mat>(out["mu"]);
      alpha_best = as<arma::mat>(out["alpha"]);
      ssqv_best = as<arma::mat>(out["ssq_var"]);
      ssq_best = ssq(j);
      sbsq_best = sbsq(j);
      pip_best = pip(j);
    }
  }

  return(Rcpp::List::create(
      Rcpp::Named("elbo") = elbo_best, Rcpp::Named("mu") = mu_best,
      Rcpp::Named("alpha") = alpha_best, Rcpp::Named("ssq_var") = ssqv_best,
      Rcpp::Named("ssq") = ssq_best, Rcpp::Named("sbsq") = sbsq_best,
      Rcpp::Named("pip") = pip_best, Rcpp::Named("elbo_vec") = elbo_store,
      Rcpp::Named("iter") = iter_store));
}
