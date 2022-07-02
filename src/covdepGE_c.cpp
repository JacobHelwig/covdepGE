#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// -----------------------------------------------------------------------------
// -----------------------------ELBO_calculator_c-------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// Calculates ELBO for a fixed response j and individual l
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x 1 vector; weights (i-th entry is the weight of the i-th
// individual with respect to the l-th individual using the l-th individual's
// bandwidth)
// X: n x p matrix; data_mat with the j-th column removed
// ssq_var, mu, alpha: p x 1 vectors; variational parameters. the k-th entry is the
// k-th parameter for the l-th individual
// ssq, sbsq: double; spike and slab variance hyperparameters
// pip: double; spike and slab probability of inclusion
// -----------------------------RETURNS-----------------------------------------
// Returns: double; ELBO for the l-th individual and j-th column fixed as the
// response
// -----------------------------------------------------------------------------
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

// -----------------------------------------------------------------------------
// -----------------------------total_elbo_c------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// for a fixed response, calculate and sum the ELBO for all individuals in the
// data
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x n matrix; weights (k,l entry is the weight of the k-th individual
// with respect to the l-th individual using the l-th individual's bandwidth)
// X: n x p matrix; data_mat with the j-th column removed
// ssq_var, mu, alpha: n x p matrices; variational parameters. the l, k
// entry is the k-th parameter for the l-th individual
// ssq, sbsq: doubles; spike and slab variance hyperparameters
// pip: double; spike and slab probability of inclusion
// -----------------------------RETURNS-----------------------------------------
// elbo_tot: double; ELBO for the l-th individual with j-th column fixed as the
// response
// -----------------------------------------------------------------------------
arma::colvec ELBO_c (const arma::colvec& y, const arma::mat& D, const arma::mat& X,
               const arma::mat& ssq_var, const arma::mat& mu,
               const arma::mat& alpha, double ssq, double sbsq, double pip){

  // get sample size
  int n = X.n_rows;

  // instantiate variable to store the ELBO for each individual
  arma::colvec elbo(n);

  // loop over all individuals
  for (int l = 0; l < n; l++){

    // calculate and store the ELBO for l-th individual
    elbo(l) = ELBO_calculator_c(y, D.col(l), X, ssq_var.row(l).t(),
         mu.row(l).t(), alpha.row(l).t(), ssq, sbsq, pip);
  }

  return elbo;
}

// -----------------------------------------------------------------------------
// -----------------------------mu_update_c-------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// for a fixed response, update mu matrix
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x n matrix; weights (k,l entry is the weight of the k-th individual
// with respect to the l-th individual using the l-th individual's bandwidth)
// X: n x p matrix; data_mat with the j-th column removed
// ssq_var, mu, alpha: n x p matrices; variational parameters. the l, k entry is
// the k-th parameter for the l-th individual
// ssq: double; spike and slab variance hyperparameter
// -----------------------------------------------------------------------------
void mu_update_c (const arma::colvec& y, const arma::mat& D, const arma::mat& X,
                  const arma::mat& ssq_var, arma::mat& mu,
                  const arma::mat& alpha, double ssq){

  // get sample size and number of parameters
  int n = X.n_rows;
  int p = X.n_cols;

  // instantiate matrices for the update loop
  arma::mat mu_stack(n, p);
  arma::mat alpha_stack(n, p);
  arma::mat X_mu_alpha(n, p);
  arma::mat X_mu_alpha_k(n, p);
  arma::mat y_k(n, p);
  arma::mat d_x_y(n, p);

  // loop over the individuals to update mu row by row
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
// mu, alpha: n x p matrices; variational parameters. the l, k entry is
// the k-th parameter for the l-th individual
// alpha_logit_term 1, 2, 3: double (1), n x p matrices (2, 3): terms used to
// calculate the logit of alpha
// upper_limit: double; during the alpha update, values of logit(alpha) greater
// than upper_limit will be assigned a probability of 1; this avoids issues
// with exponentiation of large numbers creating Infinity divided by Infinity
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

  // handle NA's due to division by infinity resulting from exponentiation of
  // large values; these probabilities are indescernible from 1

  // find large values
  double upper_limit = 9;
  arma::uvec index1 = arma::find(alpha_logit > upper_limit);

  // replace them
  alpha.elem(index1) = arma::vec(index1.n_rows, arma::fill::ones);
}

// -----------------------------------------------------------------------------
// -----------------------------cavi_c------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// The core function that performs CAVI to calculate and return the final
// variational estimates of a single regression for each of the n individuals
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x n matrix; weights (k,l entry is the weight of the k-th individual
// with respect to the l-th individual using the l-th individual's bandwidth)
// X: n x p matrix; data_mat with the j-th column removed
// mu, alpha: n x p matrices; variational parameters. the l, k entry is
// the k-th parameter for the l-th individual
// ssq, sbsq: doubles; spike and slab variance hyperparameters
// pip: double; spike and slab probability of inclusion
// tolerance: double; when the square root of the sum of squared changes in
// the elements of alpha are within tolerance, stop iterating
// max_iter: integer; maximum number of iterations
// upper_limit: double; during the alpha update, values of logit(alpha) greater
// than upper_limit will be assigned a probability of 1
// -----------------------------RETURNS-----------------------------------------
// var_alpha: n x p matrix; final alpha values
// var_ELB: double; final value of ELBO summed across all individuals
// converged_iter: integer; number of iterations to reach convergence
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List cavi_c(const arma::colvec& y, const arma::mat& D,
                  const arma::mat& X, const arma::mat& mu0,
                  const arma::mat& alpha0, double ssq, double sbsq, double pip,
                  double elbo_tol, double alpha_tol, int max_iter) {

  // matrices and vectors for updated variational parameters and hyperparameters
  arma::mat mu = mu0;
  arma::mat alpha = alpha0;

  // matrices for tracking the convergence of alpha parameters
  arma::mat alpha_last;
  arma::mat change_alpha;

  // variational variance of regression coefficients update
  arma::mat ssq_var = ssq / ((arma::pow(X, 2).t() * D).t() + 1 / sbsq);

  // 1st and 3rd term of the alpha update, denominator of the second term
  double alpha1 = log(pip / (1 - pip));
  arma::mat alpha3 = arma::log(arma::sqrt(ssq_var / (ssq * sbsq)));
  arma::mat alpha2_denom = 2 * ssq_var;

  // define last_elbo and cur_elbo for monitoring ELBO progress
  double cur_elbo;
  double last_elbo;

  // elbo_tol less than 0 indicates that elbo should not be monitored during CAVI
  if (elbo_tol >= 0){
    cur_elbo = 0;
    last_elbo = sum(ELBO_c(y, D, X, ssq_var, mu, alpha, ssq, sbsq, pip));
  }

  // CAVI loop (optimize variational parameters)
  for (int k = 0; k < max_iter; k++){

    // mu update
    mu_update_c(y, D, X, ssq_var, mu, alpha, ssq);

    // alpha update

    // save the last value of alpha and update it
    alpha_last = alpha;
    alpha_update_c(ssq_var, mu, alpha, alpha1, alpha2_denom, alpha3, ssq, sbsq,
                   pip);

    // check for convergence of ELBO
    if (elbo_tol >= 0){
      cur_elbo = sum(ELBO_c(y, D, X, ssq_var, mu, alpha, ssq, sbsq, pip));

      // if new elbo is within the tolerance or has increased, end updates
      if ((cur_elbo - last_elbo < elbo_tol)){
        break;
      }

      // otherwise, the last_elbo becomes the current elbo
      last_elbo = cur_elbo;
    }

    // calculate element-wise change in alpha
    change_alpha = alpha - alpha_last;

    // if the Frobenius norm of the change in alpha is within the tolerance,
    // end updates
    if (sqrt(arma::accu(arma::pow(change_alpha, 2))) < alpha_tol){
      break;
    }
  }

  // calculate ELBO across n individuals
  arma::colvec ELBO = ELBO_c(y, D, X, ssq_var, mu, alpha, ssq, sbsq, pip);

  // return final variational parameters and ELBO
  return(Rcpp::List::create(
      Rcpp::Named("mu") = mu, Rcpp::Named("alpha") = alpha,
      Rcpp::Named("ssq_var") = ssq_var, Rcpp::Named("elbo") = ELBO));
}

// -----------------------------------------------------------------------------
// -----------------------------grid_search_c-----------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// for a fixed response, run CAVI for each individual across a grid of
// hyperparameters and return the ELBO for each of the grid points
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x n matrix; weights (k,l entry is the weight of the k-th individual
// with respect to the l-th individual using the l-th individual's bandwidth)
// X: n x p matrix; data_mat with the j-th column removed
// mu, alpha: n x p; matrices of variational parameters. the l, k entry
// is the k-th parameter for the l-th individual
// ssq_vec: n_param x 1 vector; error term variance candidates
// sbsq_vec: n_param x 1 vector; slab variance candidates
// pi_vec: n_param x 1 vector; candidate spike and slab probabilities of inclusion
// tolerance: double; when the square root of the sum of squared changes in
// the elements of alpha are within tolerance, stop iterating
// -----------------------------RETURNS-----------------------------------------
// returns list with 2 values:
// elbo: n_param x 1 vector; ELBO for each hyperparameter candidate
// num_converged: integer; the number of hyperparameter candidates for which
// the CAVI converged
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List grid_search_c(const arma::colvec& y, const arma::mat& D,
                         const arma::mat& X, const arma::mat& mu,
                         const arma::mat& alpha, const arma::colvec& ssq,
                         const arma::colvec& sbsq, const arma::colvec& pip,
                         double elbo_tol, double alpha_tol, int max_iter){

  // get number of grid points
  int n_param = pip.n_elem;

  // create storage for the hyperparameters, variational parameters, and ELBO
  Rcpp::List elbo_list(n_param);
  Rcpp::List mu_list(n_param), alpha_list(n_param), ssqv_list(n_param);

  // instantiate a list to store the result of cavi_c
  Rcpp::List out;

  // perform CAVI for each grid point
  for (int j = 0; j < n_param; j++){

    // run CAVI
    out = cavi_c(y, D, X, mu, alpha, ssq(j), sbsq(j), pip(j), elbo_tol,
                 alpha_tol, max_iter);

    // store the final ELBO and variational parameters
    elbo_list(j) = as<arma::colvec>(out["elbo"]);
    mu_list(j) = as<arma::mat>(out["mu"]);
    alpha_list(j) = as<arma::mat>(out["alpha"]);
    ssqv_list(j) = as<arma::mat>(out["ssq_var"]);
  }

  return(Rcpp::List::create(
      Rcpp::Named("elbo") = elbo_list, Rcpp::Named("mu") = mu_list,
      Rcpp::Named("alpha") = alpha_list, Rcpp::Named("ssq_var") = ssqv_list));
}
