#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// Calculates ELBO for a fixed response and individual
// [[Rcpp::export]]
double ELBO_calculator_c (const arma::colvec& y, const arma::colvec& W,
                        const arma::mat& X_mat, const arma::colvec& S_sq,
                        const arma::colvec& mu, const arma::colvec& alpha,
                        double sigmasq, double sigmabeta_sq, double pi_est) {

  // square of: y minus X matrix-multiplied with element-wise product of mu and alpha
  arma::colvec y_Xmualpha_sq = arma::pow(y - (X_mat * (mu % alpha)), 2);

  // square of mu vector
  arma::colvec mu_sq = arma::pow(mu, 2);

  // calculate each of the terms that sum to ELBO
  double t1 = -sum(W % y_Xmualpha_sq) / (2 * sigmasq);
  double t2 = (-sum(W % (arma::pow(X_mat, 2) * (alpha % (mu_sq + S_sq) - arma::pow(alpha, 2) % mu_sq))) /
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
  double n = X_mat.n_rows;

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

/*** R
# generate data
set.seed(3)
n <- 5
p <- 3
data_mat <- matrix(sample(-3:3, n * (p + 1), T), n, p + 1)
Z <- matrix(sample(1:n, n), n, 1)
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
alpha_mat <- matrix(0.2, n, p)
sigmabeta_sq <- 3
sigmasq <- 1
S_sq <- t(sigmasq * (t(X_mat^2) + 1 / sigmabeta_sq)^(-1))
mu_mat <- matrix(0, n, p, byrow = TRUE)

# Setting hyperparameter values for sigmasq and the probability of inclusion
# according to the Carbonetto Stephens model
idmod <- varbvs(X_mat, y, Z = Z[, 1], verbose = FALSE)
sigmasq <- mean(idmod$sigma)
pi_est <- mean(1 / (1 + exp(-idmod$logodds)))

ELBO_calculator <- function(y, W, X_mat, S_sq, mu, alpha, sigmasq, sigmabeta_sq, pi_est) {
  n <- nrow(X_mat)
  p <- ncol(X_mat)
  mu <- matrix(mu, p, 1)
  alpha <- matrix(alpha, p, 1)

  s <- matrix(S_sq, p, 1)
  mu_alpha <- matrix(mu * alpha, p, 1)
  W <- matrix(W, n, 1)
  t1 <- -sum(W * (y - X_mat %*% mu_alpha)^2) / (2 * sigmasq)
  t2 <- -sum(W * ((X_mat)^2 %*% (alpha * (mu^2 + s) - alpha^2 * mu^2))) / (2 * sigmasq)
  t3 <- sum(alpha * ((1 + log(s)))) / 2
  t4 <- -sum(alpha * log((alpha + 0.000001) / pi_est) + (1 - alpha) * log((1 - alpha + 0.000001) / (1 - pi_est)))
  t5 <- -sum(alpha * ((mu^2 + s) / (2 * sigmasq * sigmabeta_sq) + log(sigmasq * sigmabeta_sq) / 2))
  t6 <- sum(0.5 * log(1 / (2 * pi * sigmasq)))
  t <- t1 + t2 + t3 + t4 + t5 + t6
  return(t1 + t2 + t3 + t4 + t5 + t6)
}

elbo_R <- 0
elbo_C <- 0
for (l in 1:n){
  elbo_R <- elbo_R + ELBO_calculator(y, D[, l], X_mat, S_sq[l, ], mu_mat[l, ],
                                     alpha_mat[l, ], sigmasq, sigmabeta_sq, pi_est)
  elbo_C <- elbo_C + ELBO_calculator_c(y, D[, l], X_mat, S_sq[l, ], mu_mat[l, ],
                                       alpha_mat[l, ], sigmasq, sigmabeta_sq, pi_est)
}
elbo_R
elbo_C

tot_ELBO_R <- function(){
  for (l in 1:n){
    ELBO_calculator(y, D[, l], X_mat, S_sq[l, ], mu_mat[l, ], alpha_mat[l, ],
                    sigmasq, sigmabeta_sq, pi_est)
  }
}

microbenchmark::microbenchmark(total_ELBO(y, D, X_mat, S_sq, mu_mat, alpha_mat, sigmasq, sigmabeta_sq, pi_est),
                               tot_ELBO_R())

*/

