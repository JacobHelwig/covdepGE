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
void mu_update_c (const arma::colvec& y, const arma::mat& W,
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
    d_x_y = arma::repmat(W.col(l), 1, p) % X_mat % y_k;

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
set.seed(1)
alpha_mat <- matrix(runif(n * p), n, p)
sigmabeta_sq <- 3
sigmasq <- 1
S_sq <- t(sigmasq * (t(X_mat^2) + 1 / sigmabeta_sq)^(-1))
mu_mat <- matrix(rnorm(n * p), n, p)

# Setting hyperparameter values for sigmasq and the probability of inclusion
# according to the Carbonetto Stephens model
idmod <- varbvs(X_mat, y, Z = Z[, 1], verbose = FALSE)
sigmasq <- mean(idmod$sigma)
pi_est <- mean(1 / (1 + exp(-idmod$logodds)))

#-------------------------------------------------------------------------------
#------------------------------ELBO CALCULATION---------------------------------
#-------------------------------------------------------------------------------

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

total_ELBO_c(y, D, X_mat, S_sq, mu_mat, alpha_mat, sigmasq, sigmabeta_sq, pi_est) == elbo_R

microbenchmark::microbenchmark(total_ELBO_c(y, D, X_mat, S_sq, mu_mat, alpha_mat, sigmasq, sigmabeta_sq, pi_est),
                               tot_ELBO_R())


#-------------------------------------------------------------------------------
#----------------------------------MU UPDATE------------------------------------
#-------------------------------------------------------------------------------

# original mu update
R_mupdate <- function(){
  mu_mat2 <- matrix(NA, n, p)
  for (l in 1:n){

    # the l-th row of mu_mat, alpha_mat stacked n times
    mu_stack <- matrix(mu_mat[l, ], n, p, T)
    alpha_stack <- matrix(alpha_mat[l, ], n, p, T)

    # the element-wise product of X_mat, mu_stack, and alpha stack;
    # the i,j entry is x_i,j * mu_l,j * alpha_l,j
    X_mu_alpha <- X_mat * mu_stack * alpha_stack

    # the k-th column is the rowSums of X_mu_alpha minus the k-th column of
    # X_mu_alpha (accounts for m \neq k in summation)
    X_mu_alpha_k <- matrix(rowSums(X_mu_alpha), n, p) - X_mu_alpha

    # the k-th column is y minus the k-th column of X_mu_alpha_k
    y_k <- matrix(y, n, p) - X_mu_alpha_k

    # the k-th column is d_:,l * x_:,k * y_k_:,k
    d_x_y <- D[ , l] * X_mat * y_k

    # the update of the l-th row of mu
    mu_mat2[l, ] <- S_sq[l, ] / sigmasq * colSums(d_x_y)

  }
  mu_mat2
}


# c update
mu_update_c(y, D, X_mat, S_sq, mu_mat, alpha_mat, sigmasq)

microbenchmark::microbenchmark(mu_update_c(y, D, X_mat, S_sq, mu_mat, alpha_mat, sigmasq),
                               R_mupdate())

#-------------------------------------------------------------------------------
#----------------------------------ALPHA UPDATE---------------------------------
#-------------------------------------------------------------------------------

alpha_mat2 <- matrix(NA, n, p)

# R alpha update
upper_limit <- 9;

# 1st and 3rd term of the alpha update, denominator of the second term
alpha_logit_term1 <- logit(pi_est)
alpha_logit_term3 <- log(sqrt(S_sq / (sigmasq * sigmabeta_sq)))
alpha_logit_term2_denom <- (2 * S_sq)

# calculate the logit of alpha
alpha_logit <- (alpha_logit_term1 +
                  (mu_mat^2 / alpha_logit_term2_denom) +
                  alpha_logit_term3)

# transform from logit to probabilities of inclusion; update alpha_mat
exp_logit <- exp(alpha_logit)
alpha_mat2 <- exp_logit / (1 + exp_logit)

# handle NA's due to division by infinity resulting from exponentiation of
# large values; these probabilities are indescernible from 1
#alpha_mat[is.infinite(exp_logit)] <- 1
alpha_mat2[alpha_logit > upper_limit] <- 1

alpha_mat2
# c update
alpha_mat
alpha_update_c(mu_mat, alpha_mat, alpha_logit_term1, alpha_logit_term2_denom,
               alpha_logit_term3)
alpha_mat
*/

