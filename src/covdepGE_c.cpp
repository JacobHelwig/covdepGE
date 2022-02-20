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
// X_mat: n x p matrix; data_mat with the j-th column removed
// S_sq, mu, alpha: p x 1 vectors; variational parameters. the k-th entry is the
// k-th parameter for the l-th individual
// sigmasq, sigmabeta_sq: double; spike and slab variance hyperparameters
// pi_est: double; spike and slab probability of inclusion
// -----------------------------RETURNS-----------------------------------------
// Returns: double; ELBO for the l-th individual and j-th column fixed as the
// response
// -----------------------------------------------------------------------------
double ELBO_calculator_c (const arma::colvec& y, const arma::colvec& D,
                          const arma::mat& X_mat, const arma::colvec& S_sq,
                          const arma::colvec& mu, const arma::colvec& alpha,
                          double sigmasq, double sigmabeta_sq, double pi_est) {

  // square of: y minus X matrix-multiplied with element-wise product of mu and
  // alpha
  arma::colvec y_Xmualpha_sq = arma::pow(y - (X_mat * (mu % alpha)), 2);

  // square of mu vector
  arma::colvec mu_sq = arma::pow(mu, 2);

  // calculate each of the terms that sum to ELBO
  double t1 = -sum(D % y_Xmualpha_sq) / (2 * sigmasq);
  double t2 = (-sum(D % (arma::pow(X_mat, 2) * (alpha % (mu_sq + S_sq) -
               arma::pow(alpha, 2) % mu_sq))) /(2 * sigmasq));
  double t3 = sum(alpha % (1 + arma::log(S_sq))) / 2;
  double t4 = (-sum(alpha % arma::log((alpha + 0.000001) / pi_est) +
               (1 - alpha) % log((1 - alpha + 0.000001) / (1 - pi_est))));
  double t5 = (-sum(alpha % ((mu_sq + S_sq) / (2 * sigmasq * sigmabeta_sq) +
               log(sigmasq * sigmabeta_sq) / 2)));
  double t6 = 0.5 * log(1 / (2 *  + M_PI * sigmasq));

  return(t1 + t2 + t3 + t4 + t5 + t6);
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
// X_mat: n x p matrix; data_mat with the j-th column removed
// S_sq, mu, alpha: n x p matrices; variational parameters. the l, k
// entry is the k-th parameter for the l-th individual
// sigmasq, sigmabeta_sq: doubles; spike and slab variance hyperparameters
// pi_est: double; spike and slab probability of inclusion
// -----------------------------RETURNS-----------------------------------------
// elbo_tot: double; ELBO for the l-th individual with j-th column fixed as the
// response
// -----------------------------------------------------------------------------
double total_ELBO_c (const arma::colvec& y, const arma::mat& D,
                     const arma::mat& X_mat, const arma::mat& S_sq,
                     const arma::mat& mu_mat, const arma::mat& alpha_mat,
                     const arma::colvec& sigmasq,
                     const arma::colvec& sigmabeta_sq, double pi_est){

  // get sample size
  int n = X_mat.n_rows;

  // instantiate variable to store the total ELBO
  double elbo_tot = 0;

  // loop over all individuals
  for (int l = 0; l < n; l++){

    // calculate the ELBO for l-th individual and add it to the total ELBO
    elbo_tot += ELBO_calculator_c(y, D.col(l), X_mat, S_sq.row(l).t(),
                                  mu_mat.row(l).t(), alpha_mat.row(l).t(),
                                  sigmasq(l), sigmabeta_sq(l), pi_est);
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
// D: n x n matrix; weights (k,l entry is the weight of the k-th individual
// with respect to the l-th individual using the l-th individual's bandwidth)
// X_mat: n x p matrix; data_mat with the j-th column removed
// S_sq, mu, alpha: n x p matrices; variational parameters. the l, k entry is
// the k-th parameter for the l-th individual
// sigmasq: double; spike and slab variance hyperparameter
// -----------------------------------------------------------------------------
void mu_update_c (const arma::colvec& y, const arma::mat& D,
                  const arma::mat& X_mat, const arma::mat& S_sq, arma::mat& mu,
                  const arma::mat& alpha, const arma::colvec& sigmasq){

  // get sample size and number of parameters
  int n = X_mat.n_rows;
  int p = X_mat.n_cols;

  // instantiate matrices for the update loop
  arma::mat mu_stack(n, p);
  arma::mat alpha_stack(n, p);
  arma::mat X_mu_alpha(n, p);
  arma::mat X_mu_alpha_k(n, p);
  arma::mat y_k(n, p);
  arma::mat d_x_y(n, p);

  // loop over the individuals to update mu row by row
  for (int l = 0; l < n; l++){

    // l-th row of mu_mat, alpha_mat stacked n times
    mu_stack = arma::repmat(mu.row(l), n, 1);
    alpha_stack = arma::repmat(alpha.row(l), n ,1);

    // take the element-wise product of X_mat, mu_stack, and alpha stack
    X_mu_alpha = X_mat % mu_stack % alpha_stack;

    // the k-th column is the rowSums of X_mu_alpha minus the k-th column of
    // X_mu_alpha (accounts for m \neq k in summation)
    X_mu_alpha_k = arma::repmat(arma::sum(X_mu_alpha, 1), 1, p) - X_mu_alpha;

    // the k-th column is y minus the k-th column of X_mu_alpha_k
    y_k = arma::repmat(y, 1, p) - X_mu_alpha_k;

    // the k-th column is d_:,l * x_:,k * y_k_:,k
    d_x_y = arma::repmat(D.col(l), 1, p) % X_mat % y_k;

    // the update of the l-th row of mu
    mu.row(l) = S_sq.row(l) % sum(d_x_y, 0) / sigmasq(l);
  }
}

// -----------------------------------------------------------------------------
// -----------------------------alpha_update_c----------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// for a fixed response, update alpha matrix
// -----------------------------ARGUMENTS---------------------------------------
// mu_mat, alpha_mat: n x p matrices; variational parameters. the l, k entry is
// the k-th parameter for the l-th individual
// alpha_logit_term 1, 2, 3: double (1), n x p matrices (2, 3): terms used to
// calculate the logit of alpha
// upper_limit: double; during the alpha update, values of logit(alpha) greater
// than upper_limit will be assigned a probability of 1; this avoids issues
// with exponentiation of large numbers creating Infinity divided by Infinity
// -----------------------------------------------------------------------------
void alpha_update_c(const arma::mat& S_sq, const arma::mat& mu, arma::mat& alpha,
                    const arma::colvec& sigmasq,
                    const arma::colvec& sigmabeta_sq, double pi_est,
                    double upper_limit = 9){

  // get dimensions of the data
  int p = S_sq.n_cols;

  // 1st and 3rd term of the alpha update, denominator of the second term
  double alpha_logit_term1 = log(pi_est / (1 - pi_est));
  arma::mat alpha_logit_term3 = arma::log(arma::sqrt(S_sq /
    arma::repmat(sigmasq % sigmabeta_sq, 1, p)));
  arma::mat alpha_logit_term2_denom = 2 * S_sq;

  // calculate the logit of alpha
  arma::mat alpha_logit = (alpha_logit_term1
                             + (arma::pow(mu, 2) / alpha_logit_term2_denom)
                             + alpha_logit_term3);

  // transform from logit to probabilities of inclusion; update alpha_mat
  arma::mat exp_logit = arma::exp(alpha_logit);
  alpha = exp_logit / (1 + exp_logit);

  // handle NA's due to division by infinity resulting from exponentiation of
  // large values; these probabilities are indescernible from 1

  // find large values
  arma::uvec index1 = arma::find(alpha_logit > upper_limit);

  // replace them
  alpha.elem(index1) = arma::vec(index1.n_rows, arma::fill::ones);
}

// -----------------------------------------------------------------------------
// -----------------------------sigma_update_c----------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// Update the residual and slab variance term for each individual using MAPE
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x n matrix; weights (k,l entry is the weight of the k-th individual
// with respect to the l-th individual using the l-th individual's bandwidth)
// X_mat: n x p matrix; data_mat with the j-th column removed
// mu_mat, alpha_mat, S_sq: n x p matrices; variational parameters. the l, k
// entry is the k-th parameter for the l-th individual
// sigmasq, sigmabeta_sq: doubles; spike and slab variance hyperparameters
// -----------------------------RETURNS-----------------------------------------
// returns list with two values:
// sigmasq: n x 1 vector; updated residual variance
// sigmabeta_sq: n x 1 vector; updated slab variance
// -----------------------------------------------------------------------------
void sigma_update_c(const arma::colvec& y, const arma::mat& D,
                    const arma::mat& X_mat, arma::mat& S_sq,
                    const arma::mat& mu_mat, const arma::mat& alpha_mat,
                    arma::colvec& sigmasq, arma::colvec& sigmabeta_sq,
                    bool update_sigmasq, bool update_sigmabetasq) {

  // calculate terms needed for both updates

  // get dimensions of the data
  int n = X_mat.n_rows;
  int p = X_mat.n_cols;

  // calculate mu^2
  arma::mat mu_sq = arma::pow(mu_mat, 2);

  // calculate the third value in the numerator of the sigmasq update for each
  // individual; the l-th value is for individual l
  arma::colvec num_term3 = arma::sum(alpha_mat % (S_sq + mu_sq), 1) / sigmabeta_sq;

  // calculate the denominator for each individual; l-th value is for individual l
  arma::colvec denom = arma::sum(alpha_mat, 1) + n;

  // if sigmasq is to be updated, update it
  if (update_sigmasq){

    // terms for sigmasq update:

    // calulate alpha^2
    arma::mat alpha_sq = arma::pow(alpha_mat, 2);

    // calculate expected value of beta for each individual; l-th row is
    // E(beta) for individual l
    arma::mat rho = mu_mat % alpha_mat;

    // find fitted values using expected value of beta for each individual; l-th
    // column is fitted values for individual l
    arma::mat fitted = X_mat * rho.t();

    // calculate the squared residuals for each of the fitted values for each
    // individual; l-th column is residuals for individual l
    arma::mat resid2 = arma::pow(arma::repmat(y, 1, n) - fitted, 2);

    // calculate the sum of the weighted squared residuals for each individual;
    // l-th value is the SWSR for individual l
    arma::colvec resid_w = arma::sum(resid2 % D, 0).t();

    // calculate the second values in the numerator for each individual; the l-th
    // row is for individual l
    arma::mat num_term2 = alpha_mat % S_sq + alpha_mat % mu_sq - alpha_sq % mu_sq;

    // vector for storing weights
    arma::colvec weights;

    // matrix for storing the weighted version of X
    arma::mat X_w;

    // vector for storing the diagonal of the weighted XtX
    arma::colvec XtX_w;

    // double for storing the second numerator term for the l-th individual
    double num_term2_l;

    // iterate over the individuals to update each error variance
    for (int l = 0; l < n; l++){

      // sigmasq update for individual l:

      // calculate weighted version of X
      weights = arma::sqrt(D.col(l));
      X_w = X_mat % arma::repmat(weights, 1, p);

      // diagonal elements of X transpose X weighted
      XtX_w = arma::diagvec(X_w.t() * X_w);

      // second numerator term
      num_term2_l = arma::as_scalar(XtX_w.t() * num_term2.row(l).t());

      // apply update
      sigmasq(l) = (resid_w(l) + num_term2_l + num_term3(l)) / denom(l);
    }
  }

  // if sigmabeta_sq is to be updated, update it
  if(update_sigmabetasq){

    // if sigmasq was updated, update S_sq
    if(update_sigmasq){
      S_sq = arma::repmat(sigmasq, 1, p) /
        ((arma::pow(X_mat, 2).t() * D).t() + 1 / arma::repmat(sigmabeta_sq, 1, p));
    }

    // terms for sigmabeta_sq update

    // numerator is num_term3 without the division by sigmabeta_sq
    num_term3 = num_term3 % sigmabeta_sq;

    // denominator is denom without summing n, and also scaled by sigmasq
    denom = sigmasq % (denom - n);

    // update the slab variance
    sigmabeta_sq = num_term3 / denom;
  }
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
// X_mat: n x p matrix; data_mat with the j-th column removed
// mu_mat, alpha_mat: n x p matrices; variational parameters. the l, k entry is
// the k-th parameter for the l-th individual
// sigmasq, sigmabeta_sq: doubles; spike and slab variance hyperparameters
// pi_est: double; spike and slab probability of inclusion
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
Rcpp::List cavi_c(const arma::colvec& y, const arma::mat& D,
                  const arma::mat& X_mat, const arma::mat& mu_mat,
                  const arma::mat& alpha_mat, const arma::colvec& sigmasq_vec,
                  bool update_sigmasq, const arma::colvec& sigmabeta_sq_vec,
                  bool update_sigmabetasq, double pi_est, double tolerance,
                  int max_iter, double upper_limit = 9) {

  // get dimensions of the data
  int p = X_mat.n_cols;

  // matrices and vectors for updated variational parameters and hyperparameters
  arma::mat mu = mu_mat;
  arma::mat alpha = alpha_mat;
  arma::colvec sigmasq = sigmasq_vec;
  arma::colvec sigmabeta_sq = sigmabeta_sq_vec;

  // matrices for tracking the convergence of alpha parameters
  arma::mat alpha_last;
  arma::mat change_alpha;

  // integer for tracking the iteration at which convergence is reached
  int converged_iter = max_iter;

  // matrix for storing S_sq
  arma::mat S_sq;

  // CAVI loop (optimize variational parameters)
  for (int k = 0; k < max_iter; k++){

    // S_sq update
    S_sq = arma::repmat(sigmasq, 1, p) /
      ((arma::pow(X_mat, 2).t() * D).t() + 1 / arma::repmat(sigmabeta_sq, 1, p));

    // mu update
    mu_update_c(y, D, X_mat, S_sq, mu, alpha, sigmasq);

    // alpha update

    // save the last value of alpha and update it
    alpha_last = alpha;
    alpha_update_c(S_sq, mu, alpha, sigmasq, sigmabeta_sq, pi_est);

    // calculate change in alpha
    change_alpha = alpha - alpha_last;

    // if the square root of the sum of squared changes in alpha is within the
    // tolerance, break from the for loop
    if (sqrt(arma::accu(arma::pow(change_alpha, 2))) < tolerance){
      converged_iter = k;
      break;
    }

    // update the variance terms using MAPE
    if (update_sigmasq | update_sigmabetasq){
      sigma_update_c(y, D, X_mat, S_sq, mu, alpha, sigmasq, sigmabeta_sq,
                     update_sigmasq, update_sigmabetasq);
    }
  }

  // calculate ELBO across n individuals
  double ELBO = total_ELBO_c(y, D, X_mat, S_sq, mu, alpha, sigmasq,
                             sigmabeta_sq, pi_est);

  // return final mu matrix, alpha matrix, the final ELBO, the number of
  // iterations to converge, and the fitted variance hyperparameters
  return(Rcpp::List::create(Rcpp::Named("var_mu") = mu,
                            Rcpp::Named("var_alpha") = alpha,
                            Rcpp::Named("var_elbo") = ELBO,
                            Rcpp::Named("converged_iter") = converged_iter,
                            Rcpp::Named("sigmasq") = sigmasq,
                            Rcpp::Named("sigmabeta_sq") = sigmabeta_sq));
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
// X_mat: n x p matrix; data_mat with the j-th column removed
// mu_mat, alpha_mat: n x p; matrices of variational parameters. the l, k entry
// is the k-th parameter for the l-th individual
// sigmasq_vec: n_param x 1 vector; error term variance candidates
// sigmabeta_sq_vec: n_param x 1 vector; slab variance candidates
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
                         const arma::mat& X_mat, const arma::mat& mu_mat,
                         const arma::mat& alpha_mat,
                         const arma::mat& sigmasq_vec, bool update_sigmasq,
                         const arma::mat& sigmabeta_sq_vec,
                         bool update_sigmabetasq, const arma::colvec& pi_vec,
                         double tolerance, int max_iter, double upper_limit = 9){

  // get the number of grid points
  int n_param = pi_vec.n_elem;

  // storage for the best ELBO and new elbo
  double elbo_best = -1;
  double elbo_new;

  // storage for the hyperparameters and variational parameters corresponding to
  // the current best ELBO
  arma::mat mu_best;
  arma::mat alpha_best;
  arma::colvec sigmasq_best;
  arma::colvec sigmabetasq_best;
  double pi_best = -1;

  // storage for best iterations
  int iterations_best = -1;

  // instantiate a list to store the result of cavi_c
  Rcpp::List out;

  // perform CAVI for each grid point
  for (int j = 0; j < n_param; j++){

    // run CAVI
    out = cavi_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq_vec.col(j),
                 update_sigmasq, sigmabeta_sq_vec.col(j), update_sigmabetasq,
                 pi_vec(j), tolerance, max_iter, upper_limit);
    elbo_new = as<double>(out["var_elbo"]);

    // if the new ELBO is greater than the current best, update the best ELBO and
    // parameters and whether or not the model converged
    if ((elbo_best <  elbo_new or (std::isnan(elbo_best) and !std::isnan(elbo_new))) or j == 0){

      elbo_best = as<double>(out["var_elbo"]);
      mu_best = as<arma::mat>(out["var_mu"]);
      alpha_best = as<arma::mat>(out["var_alpha"]);
      sigmasq_best = as<arma::colvec>(out["sigmasq"]);
      sigmabetasq_best = as<arma::colvec>(out["sigmabeta_sq"]);
      pi_best = pi_vec(j);
      iterations_best = as<int>(out["converged_iter"]);
    }
  }

  return(Rcpp::List::create(Rcpp::Named("elbo") = elbo_best,
                            Rcpp::Named("mu") = mu_best,
                            Rcpp::Named("alpha") = alpha_best,
                            Rcpp::Named("sigmasq") = sigmasq_best,
                            Rcpp::Named("sigmabeta_sq") = sigmabetasq_best,
                            Rcpp::Named("pi") = pi_best,
                            Rcpp::Named("iterations") = iterations_best));
}
