#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// Calculates ELBO for a fixed response and individual
// [[Rcpp::export]]
double ELBO_calculator_c (const arma::colvec& y, const arma::colvec& D,
                        const arma::mat& X_mat, const arma::colvec& S_sq,
                        const arma::colvec& mu, const arma::colvec& alpha,
                        double sigmasq, double sigmabeta_sq, double pi_est) {

  // square of: y minus X matrix-multiplied with element-wise product of mu and alpha
  arma::colvec y_Xmualpha_sq = arma::pow(y - (X_mat * (mu % alpha)), 2);

  // square of mu vector
  arma::colvec mu_sq = arma::pow(mu, 2);

  // calculate each of the terms that sum to ELBO
  double t1 = -sum(D % y_Xmualpha_sq) / (2 * sigmasq);
  double t2 = (-sum(D % (arma::pow(X_mat, 2) * (alpha % (mu_sq + S_sq) - arma::pow(alpha, 2) % mu_sq))) /
               (2 * sigmasq));
  double t3 = sum(alpha % (1 + arma::log(S_sq))) / 2;
  double t4 = (-sum(alpha % arma::log((alpha + 0.000001) / pi_est) +
               (1 - alpha) % log((1 - alpha + 0.000001) / (1 - pi_est))));
  double t5 = (-sum(alpha % ((mu_sq + S_sq) / (2 * sigmasq * sigmabeta_sq) +
               log(sigmasq * sigmabeta_sq) / 2)));
  double t6 = 0.5 * log(1 / (2 *  + M_PI * sigmasq));

  return(t1 + t2 + t3 + t4 + t5 + t6);
}

// for a fixed response, calculate and sum the ELBO for all individuals in the data
// [[Rcpp::export]]
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

// mu update: for a fixed response, update mu matrix
//[[Rcpp::export]]
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

// alpha update: for a fixed response, update alpha matrix
//[[Rcpp::export]]
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
  arma::uvec index1 = arma::find(alpha_logit > upper_limit); // find large values
  alpha.elem(index1) = arma::vec(index1.n_rows, arma::fill::ones); // replace them
}

// The core function that calculates the variational parameter updates and
// returns the final variational estimates for a single regression.
// Arguments:
// y: n by 1 response, i.e the j th variable whose conditional distribution given the
// remaining variables is being calculated.
// X_mat: n by p data matrix with the j-th variable removed
// sigmasq: variance of response given the parameters (homoscedastic part,
// actual variance sigma_sq/w_i)
// sigmabeta_sq: prior variance of coefficient parameter
// pi_est: estimate of spike and slab mixture proportion.
// mu_mat & alpha_mat: n by p matrices of variational parameters mu and alpha;
// the i,j-th entry corresponds to the j-th parameter for the i-th individual
//[[Rcpp::export]]
Rcpp::List cov_vsvb_c(const arma::colvec& y, const arma::mat& D,
                      const arma::mat& X_mat, const arma::mat& mu_mat,
                      const arma::mat& alpha_mat, double sigmasq,
                      double sigmabeta_sq, double pi_est, double tolerance = 1e-9,
                      int max_iter = 100, double upper_limit = 9) {

  // instantiate matrices for updated variational parameters with starting
  // values dictated by the matrices passed as arguments
  arma::mat mu = mu_mat;
  arma::mat alpha = alpha_mat;

  // matrices for tracking the convergence of alpha parameters
  arma::mat alpha_last;
  arma::mat change_alpha;

  // S_sq update
  arma::mat S_sq = (sigmasq * arma::pow(arma::pow(X_mat, 2).t() * D + 1 / sigmabeta_sq, -1)).t();

  // 1st and 3rd term of the alpha update, denominator of the second term
  double alpha_logit_term1 = log(pi_est / (1 - pi_est));
  arma::mat alpha_logit_term3 = arma::log(arma::sqrt(S_sq / (sigmasq * sigmabeta_sq)));
  arma::mat alpha_logit_term2_denom = 2 * S_sq;

  // loop to optimize variational parameters
  for (int k = 1; k < max_iter; k++){ // start counting at 1 to match

    // mu update
    mu_update_c(y, D, X_mat, S_sq, mu, alpha, sigmasq);

    // alpha update

    // save the last value of alpha and update it
    alpha_last = alpha;
    alpha_update_c(mu, alpha, alpha_logit_term1, alpha_logit_term2_denom, alpha_logit_term3);

    // calculate change in alpha
    change_alpha = alpha - alpha_last;

    // if the sum of squared changes in alpha is within the tolerance, break
    // from the for loop
    if (sqrt(arma::accu(arma::pow(change_alpha, 2))) < tolerance){
      break;
    }
  }

  // calculate ELBO across n individuals
  double ELBO_LB = total_ELBO_c(y, D, X_mat, S_sq, mu, alpha, sigmasq, sigmabeta_sq, pi_est);

  // return matrices of variational parameters and the ELBO
  return(Rcpp::List::create(Rcpp::Named("var.alpha") = alpha,
                            Rcpp::Named("var.mu") = mu,
                            Rcpp::Named("var.S_sq") = S_sq,
                            Rcpp::Named("var.elbo") = ELBO_LB));
}



/*** R
source("generate_data.R")
dat <- generate_continuous()
data_mat <- dat$data
Z <- dat$covts
n <- nrow(data_mat); p <- ncol(data_mat) - 1

# generate data
# set.seed(3)
# n <- 100
# p <- 4
# data_mat <- matrix(sample(-3:3, n * (p + 1), T), n, p + 1)
# Z <- matrix(sample(1:n, n), n, 1)
tau <- 0.1

# D is an n by n matrix of weights; the i, j entry is the similarity between
# individuals i and j
D <- matrix(1, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    D[j, i] <- dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, tau)
  }
}

# Scale weights to sum to n
for (i in 1:n) {
  D[, i] <- n * (D[, i] / sum(D[, i]))
}

resp_index <- 2

# Set variable number `resp_index` as the response
y <- data_mat[, resp_index]

# Set the remaining p variables as predictor
X_mat <- data_mat[, -resp_index]

# instantiate initial values of variational parameters
set.seed(1)
alpha_mat <- matrix(runif(n * p), n, p)
sigmabeta_sq <- 3
sigmasq <- 1
S_sq <- t(sigmasq * (t(X_mat^2) + 1 / sigmabeta_sq)^(-1))
mu_mat <- matrix(rnorm(n * p), n, p)

# Setting hyperparameter values for sigmasq and the probability of inclusion
# according to the Carbonetto Stephens model
idmod <- varbvs::varbvs(X_mat, y, Z = Z[, 1], verbose = FALSE)
sigmasq <- mean(idmod$sigma)
pi_est <- mean(1 / (1 + exp(-idmod$logodds)))

#-------------------------------------------------------------------------------
#----------------------------------MAIN FUNCTION--------------------------------
#-------------------------------------------------------------------------------

cov_vsvb <- function(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq, pi_est) {

  mu_mat2 <- rlang::duplicate(mu_mat)
  alpha_mat2 <- rlang::duplicate(alpha_mat)

  n <- nrow(X_mat); p <- ncol(X_mat)

  # threshold for calculating the reverse logit of alpha
  upper_limit <- 9

  # exit condition tolerance
  tol <- 1e-9

  change_alpha <- matrix(Inf, n, p)
  alpha_last <- matrix(Inf, n, p)

  max_iter <- 100
  iter <- 1

  # S_sq update
  S_sq2 <- t(sigmasq * (t(X_mat^2) %*% D + 1 / sigmabeta_sq)^(-1))

  # 1st and 3rd term of the alpha update, denominator of the second term
  alpha_logit_term1 <- log(pi_est / (1 - pi_est))
  alpha_logit_term3 <- log(sqrt(S_sq2 / (sigmasq * sigmabeta_sq)))
  alpha_logit_term2_denom <- (2 * S_sq2)

  # loop to optimize variational parameters
  while (sqrt(sum(change_alpha^2)) > tol & iter < max_iter) {

    # mu update
    mu_update_c(y, D, X_mat, S_sq2, mu_mat2, alpha_mat2, sigmasq)

    # alpha update

    # save the last value of alpha
    alpha_last <- rlang::duplicate(alpha_mat2)
    alpha_update_c(mu_mat2, alpha_mat2, alpha_logit_term1,
                   alpha_logit_term2_denom, alpha_logit_term3)

    # calculate change in alpha
    change_alpha <- alpha_mat2 - alpha_last

    iter <- iter + 1
  }
  print(iter)
  # calculate ELBO across n individuals
  # want to maximize this by optimizing sigma beta
  ELBO_LB <- total_ELBO_c(y, D, X_mat, S_sq2, mu_mat2, alpha_mat2, sigmasq,
                          sigmabeta_sq, pi_est)

  # return matrices of variational parameters and the ELBO
  list(var.alpha = alpha_mat2, var.mu = mu_mat2, var.S_sq = S_sq2, var.elbo = ELBO_LB)
}

# R
set.seed(1)
alpha_mat <- matrix(runif(n * p), n, p)
S_sq <- t(sigmasq * (t(X_mat^2) + 1 / sigmabeta_sq)^(-1))
mu_mat <- matrix(rnorm(n * p), n, p)
h1 <- cov_vsvb(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq, pi_est)

# C++
set.seed(1)
alpha_mat <- matrix(runif(n * p), n, p)
S_sq <- t(sigmasq * (t(X_mat^2) + 1 / sigmabeta_sq)^(-1))
mu_mat <- matrix(rnorm(n * p), n, p)
h2 <- cov_vsvb_c(y, D, X_mat, mu_mat, alpha_mat, sigmasq, sigmabeta_sq, pi_est)

all.equal(h1$var.alpha, h2$var.alpha)

*/


