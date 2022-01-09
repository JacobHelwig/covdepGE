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
                     double sigmasq, double sigmabeta_sq, double pi_est){

  // get sample size
  int n = X_mat.n_rows;

  // instantiate variable to store the total ELBO
  double elbo_tot = 0;

  // loop over all individuals
  for (int l = 0; l < n; l++){

    // calculate the ELBO for l-th individual and add it to the total ELBO
    elbo_tot += ELBO_calculator_c(y, D.col(l), X_mat, S_sq.row(l).t(),
                                  mu_mat.row(l).t(), alpha_mat.row(l).t(),
                                  sigmasq, sigmabeta_sq, pi_est);
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
                  const arma::mat& alpha, const double sigmasq){

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
    mu.row(l) = S_sq.row(l) % sum(d_x_y, 0) / sigmasq;
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
void alpha_update_c(const arma::mat& mu, arma::mat& alpha,
                    double alpha_logit_term1,
                    const arma::mat& alpha_logit_term2_denom,
                    const arma::mat& alpha_logit_term3, double upper_limit = 9){

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
// monitor_elbo: bool; if true, elbo will be tracked and returned
// monitor_period: double; dictates the periodicity of the elbo recording
// upper_limit: double; during the alpha update, values of logit(alpha) greater
// than upper_limit will be assigned a probability of 1
// -----------------------------RETURNS-----------------------------------------
// var_alpha: n x p matrix; final alpha values
// var_ELB: double; final value of ELBO summed across all individuals
// converged_iter: integer; number of iterations to reach convergence
// elbo_history: 2 x d matrix OR 0 x 0 matrix; if monitor_elbo is false, 0 x 0
// matrix; otherwise, d = floor(coverged_iter / monitor_period) +
// (coverged_iter % monitor_period == 1). The first row is the iterations at
// which the ELBO was recorded, and the second row is the ELBO measurements
// -----------------------------------------------------------------------------
//[[Rcpp::export]]
Rcpp::List cavi_c(const arma::colvec& y, const arma::mat& D,
                  const arma::mat& X_mat, const arma::mat& mu_mat,
                  const arma::mat& alpha_mat, double sigmasq,
                  double sigmabeta_sq, double pi_est, double tolerance,
                  int max_iter, bool monitor_elbo, int monitor_period,
                  double upper_limit = 9) {

  // instantiate matrices for updated variational parameters with starting
  // values dictated by the matrices passed as arguments
  arma::mat mu = mu_mat;
  arma::mat alpha = alpha_mat;

  // matrices for tracking the convergence of alpha parameters
  arma::mat alpha_last;
  arma::mat change_alpha;

  // integer for tracking the iteration at which convergence is reached
  int converged_iter = max_iter;

  // matrix for tracking elbo; construct this matrix if elbo is to be tracked
  arma::mat elbo_history;
  if (monitor_elbo){

    // create a vector for indexing the iterations at which ELBO is recorded
    arma::rowvec elbo_idx = arma::regspace<arma::rowvec>(1, monitor_period,
                                                         max_iter);

    // create the matrix for storing ELBO
    elbo_history = arma::mat(2, elbo_idx.n_elem);

    // in the first row, put the iteration indexing vector
    elbo_history.row(0) = elbo_idx;
  }

  // S_sq update
  arma::mat S_sq = (sigmasq * arma::pow(arma::pow(X_mat, 2).t() *
    D + 1 / sigmabeta_sq, -1)).t();

  // 1st and 3rd term of the alpha update, denominator of the second term
  double alpha_logit_term1 = log(pi_est / (1 - pi_est));
  arma::mat alpha_logit_term3 = arma::log(arma::sqrt(S_sq /
    (sigmasq * sigmabeta_sq)));
  arma::mat alpha_logit_term2_denom = 2 * S_sq;

  // CAVI loop (optimize variational parameters)
  for (int k = 0; k < max_iter; k++){

    // mu update
    mu_update_c(y, D, X_mat, S_sq, mu, alpha, sigmasq);

    // alpha update

    // save the last value of alpha and update it
    alpha_last = alpha;
    alpha_update_c(mu, alpha, alpha_logit_term1, alpha_logit_term2_denom,
                   alpha_logit_term3);

    // if the iteration aligns with the monitor period, record the elbo
    if ((k % monitor_period == 0) & monitor_elbo){
      elbo_history(1, k / monitor_period) = total_ELBO_c(y, D, X_mat, S_sq, mu,
                   alpha, sigmasq, sigmabeta_sq, pi_est);
    }

    // calculate change in alpha
    change_alpha = alpha - alpha_last;

    // if the square root of the sum of squared changes in alpha is within the
    // tolerance, break from the for loop
    if (sqrt(arma::accu(arma::pow(change_alpha, 2))) < tolerance){
      converged_iter = k;
      break;
    }
  }

  // calculate ELBO across n individuals
  double ELBO = total_ELBO_c(y, D, X_mat, S_sq, mu, alpha, sigmasq,
                             sigmabeta_sq, pi_est);

  // if ELBO is being tracked, remove the excess columns from elbo_history, and
  // record the ELBO value at the last iteration if it is not already recorded
  if (monitor_elbo){

    // remove excess entries
    elbo_history = elbo_history.cols(arma::find(elbo_history.row(0) <=
      converged_iter));

    // assess whether the elbo value at the last iteration is in elbo_history
    if((converged_iter - 1) % monitor_period != 0){

      // the final ELBO is not in elbo_history; add another column to
      // elbo_history with the final ELBO
      elbo_history = arma::join_rows(
        elbo_history, (arma::colvec){(double)converged_iter, ELBO});
    }
  }

  // return final alpha matrix, the final ELBO, the number of iterations to
  // converge, and the elbo history matrix
  return(Rcpp::List::create(Rcpp::Named("var_alpha") = alpha,
                            Rcpp::Named("var_elbo") = ELBO,
                            Rcpp::Named("converged_iter") = converged_iter,
                            Rcpp::Named("elbo_history") = elbo_history));
}

// -----------------------------------------------------------------------------
// -----------------------------grid_search_c-----------------------------------
// -----------------------------------------------------------------------------
// -----------------------------DESCRIPTION-------------------------------------
// for a fixed response, run CAVI for each individual across a grid of
// sigmabeta_sq and pi values and return the ELBO for each of the grid points
// -----------------------------ARGUMENTS---------------------------------------
// y: n x 1 vector; responses (j-th column of the data)
// D: n x n matrix; weights (k,l entry is the weight of the k-th individual
// with respect to the l-th individual using the l-th individual's bandwidth)
// X_mat: n x p matrix; data_mat with the j-th column removed
// mu_mat, alpha_mat: n x p; matrices of variational parameters. the l, k entry
// is the k-th parameter for the l-th individual
// sigmasq: double; spike and slab variance hyperparameter
// sigmabeta_sq_vec: n_sigma x 1 vector; spike-and-slab hyperparameter
// candidates
// pi_vec: n_pi x 1 vector; candidate spike and slab probabilities of inclusion
// tolerance: double; when the square root of the sum of squared changes in
// the elements of alpha are within tolerance, stop iterating
// max_iter: integer; maximum number of iterations
// monitor_elbo: bool; if true, elbo will be tracked and returned
// monitor_period: integer; dictates the periodicity of the elbo recording
// upper_limit: double, during the alpha update, values of logit(alpha) greater
// than upper_limit will be assigned a probability of 1
// -----------------------------RETURNS-----------------------------------------
// returns n_sigma x n_pi matrix of ELBOs
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List grid_search_c(const arma::colvec& y, const arma::mat& D,
                         const arma::mat& X_mat, const arma::mat& mu_mat,
                         const arma::mat& alpha_mat, double sigmasq,
                         const arma::colvec& sigmabeta_sq_vec,
                         const arma::colvec& pi_vec, double tolerance,
                         int max_iter, bool monitor_elbo, int monitor_period,
                         double upper_limit = 9){

  // get the number of sigmas that are being considered
  int n_sigma = sigmabeta_sq_vec.n_elem;

  // get the number of pi's that are being considered
  int n_pi = pi_vec.n_elem;

  // instantiate a matrix for storing the ELBO corresponding to each grid point
  arma::mat elbo_sigmaXpi(n_sigma, n_pi, arma::fill::zeros);

  // instantiate a list to store the result of cavi_c
  Rcpp::List out;

  // double to store the ELBO of the estimated graphs
  double elbo_graph;

  // integer to store the number of iterations it took for cavi_c to converge
  int converged_iter;

  // count the number of grid points for which convergence was attained
  int converged_ct = 0;

  // instantiate a matrix for storing pi, sigmabeta_sq values that result in
  // convergence not being attained
  arma::mat non_converged_pairs;

  // create a vector for storing the elbo_history vectors when convergence is
  // not attained, and a matrix for storing these vectors
  arma::rowvec nc_elbo_new_row;
  arma::mat non_converged_elbo;

  // for each sigma, pi pair, estimate n graphs and record the total ELBO
  // summed across these n graphs
  for (int j = 0; j < n_sigma; j++){
    for (int k = 0; k < n_pi; k++){

      // get the ELBO and convergence status for cavi_c
      out = cavi_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq_vec(j),
                   pi_vec(k), tolerance, max_iter, monitor_elbo, monitor_period,
                   upper_limit);
      elbo_graph = out["var_elbo"];

      // if cavi_c has converged, increment the convergence count
      converged_iter = out["converged_iter"];
      if (converged_iter != max_iter){
        converged_ct++;
      }
      else{

        // otherwise, convergence was not attained; add the sigabeta_sq and pi
        // values to non_converged_pairs

        // first, check to see if non_converged_pairs has any entries
        if (non_converged_pairs.n_elem == 0){

          // non_converged_pairs does not have any entries; the current pair
          // will be the first entry
          non_converged_pairs = (arma::rowvec){sigmabeta_sq_vec(j), pi_vec(k)};

          // if elbo is being tracked, then the non_converged_elbo matrix also
          // must be constructed
          if (monitor_elbo){
            non_converged_elbo = as<arma::mat>(out["elbo_history"]);
          }
        }
        else{

          // otherwise, non_converged_pairs already has some entries; add the
          // latest non-converged pair as the last row
          non_converged_pairs = arma::join_cols(
            non_converged_pairs, (arma::rowvec){sigmabeta_sq_vec(j), pi_vec(k)});

          // if elbo is being tracked, then the latest elbo history should be
          // added as the last row of non_converged_elbo
          if (monitor_elbo){
            nc_elbo_new_row = (as<arma::mat>(out["elbo_history"])).row(1);
            non_converged_elbo = arma::join_cols(non_converged_elbo,
                                                 nc_elbo_new_row);
          }
        }
      }

      // add the ELBO to the matrix of ELBOs
      elbo_sigmaXpi(j, k) = elbo_graph;
    }
  }

  return (Rcpp::List::create(Rcpp::Named("elbo_grid") = elbo_sigmaXpi,
                            Rcpp::Named("num_converged") = converged_ct,
                            Rcpp::Named("non_converged_pairs") =
                              non_converged_pairs,
                            Rcpp::Named("non_converged_elbo") =
                              non_converged_elbo));
}
