source("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/R/covdepGE_R.R")

# function for performing cavi
cavi0 <- function(y, D, X_mat, mu_mat, alpha_mat, sigmasq, update_sigmasq,
                  sigmabeta_sq, update_sigmabetasq, pi_est, tolerance,
                  max_iter, upper_limit = 9) {

  mean_alpha_tracker <- mean_mu_tracker <- mean_ssq_tracker <- t1_tracker <-
    t2_tracker <- t3_tracker <- denom_tracker <- rep(NA, max_iter)

  ELBO_tracker <- rep(NA, ceiling(max_iter / 5))

  n <- nrow(X_mat)
  p <- ncol(X_mat)

  # instantiate matrices for updated variational parameters with starting
  # values dictated by the matrices passed as arguments
  mu <- matrix(mu_mat, n, p)
  alpha <- matrix(alpha_mat, n, p)

  # matrices for tracking the convergence of alpha parameters
  alpha_last <- matrix(NA, n, p)
  change_alpha <- matrix(NA, n, p)

  # integer for tracking the iteration at which convergence is reached
  converged_iter <- max_iter

  # CAVI loop (optimize variational parameters)
  for (k in 1:max_iter){

    mean_alpha_tracker[k] <- mean(alpha)
    mean_ssq_tracker[k] <- mean(sort(sigmasq, T)[1:20])
    mean_mu_tracker[k] <- mean((abs(mu[order(rowMeans(abs(mu)), decreasing = T)[1:20], ])))

    if (k %% 5 == 0) ELBO_tracker[k] <- total_ELBO_R(y, D, X_mat, S_sq, mu,
                                                     alpha, sigmasq, sigmabeta_sq,
                                                     pi_est)

    # S_sq update
    S_sq <- matrix(sigmasq, n, p) / (t(t(X_mat^2) %*% D) + 1 /
                                       matrix(sigmabeta_sq, n, p))

    # mu update
    mu <- mu_update_R(y, D, X_mat, S_sq, mu, alpha, sigmasq);
    #mu_update_c(y, D, X_mat, S_sq, mu, alpha, sigmasq);

    # alpha update

    # save the last value of alpha and update it
    alpha_last <- matrix(alpha, n, p)
    alpha <- alpha_update_R(S_sq, mu, alpha, sigmasq, sigmabeta_sq, pi_est)
    #alpha_update_c(S_sq, mu, alpha, sigmasq, sigmabeta_sq, pi_est)

    # calculate change in alpha
    change_alpha <- alpha - alpha_last;

    # if the square root of the sum of squared changes in alpha is within the
    # tolerance, break from the for loop
    if (sqrt(sum((change_alpha^2))) < tolerance){
      converged_iter <- k
      break;
    }

    # update the variance terms using MAPE
    if (update_sigmasq | update_sigmabetasq){
      sigma_update <- sig_upd(y, D, X_mat, S_sq, mu, alpha, sigmasq,
                                     sigmabeta_sq)
      # sigma_update_c(y, D, X_mat, S_sq, mu, alpha, sigmasq,
      #                sigmabeta_sq, update_sigmasq, update_sigmabetasq)
      if (update_sigmasq) sigmasq <- sigma_update$sigmasq
      if (update_sigmabetasq) sigmabeta_sq <- sigma_update$sigmabeta_sq
      t1_tracker[k] <- sigma_update$t1
      t2_tracker[k] <- sigma_update$t2
      t3_tracker[k] <- sigma_update$t3
      denom_tracker[k] <- sigma_update$den
    }

    if (any(is.na(sigmasq))){
      warning(k); break
    }
  }

  # calculate ELBO across n individuals
  ELBO <- total_ELBO_R(y, D, X_mat, S_sq, mu, alpha, sigmasq, sigmabeta_sq, pi_est)
  #ELBO <- total_ELBO_c(y, D, X_mat, S_sq, mu, alpha, sigmasq, sigmabeta_sq, pi_est)

  # return final alpha matrix, the final ELBO, the number of iterations to
  # converge, and the elbo history matrix
  return(list(var_mu = mu, var_alpha = alpha, var_elbo = ELBO,
              converged_iter = converged_iter, sigmasq = sigmasq,
              sigmabeta_sq = sigmabeta_sq, mat = na.omit(mean_alpha_tracker),
              mmt = na.omit(mean_mu_tracker), mst = na.omit(mean_ssq_tracker),
              elb = na.omit(ELBO_tracker), t1 = na.omit(t1_tracker),
              t2 = na.omit(t2_tracker), t3 = na.omit(t3_tracker),
              den = na.omit(denom_tracker)))
}

sig_upd <- function(y, D, X_mat, S_sq, mu_mat, alpha_mat, sigmasq,
                    sigmabeta_sq) {

  # get dimensions of the data
  n <- nrow(X_mat)
  p <- ncol(X_mat)

  # terms for sigmasq update:

  # calulate alpha^2 and mu^2
  alpha_sq <- alpha_mat^2
  mu_sq <- mu_mat^2

  # calculate expected value of beta for each individual; l-th row is
  # E(beta) for individual l
  rho <- mu_mat * alpha_mat

  # find fitted values using expected value of beta for each individual; l-th
  # column is fitted values for individual l
  fitted <- X_mat %*% t(rho)

  # calculate the squared residuals for each of the fitted values for each
  # individual; l-th column is residuals for individual l
  resid2 <- (matrix(y, n, n) - fitted)^2

  # calculate the sum of the weighted squared residuals for each individual;
  # l-th value is the SWSR for individual l
  resid_w <- colSums(resid2 * D)

  # calculate the second values in the numerator for each individual; the l-th
  # row is for individual l
  num_term2 <- alpha_mat * S_sq + alpha_mat * mu_sq - alpha_sq * mu_sq

  # calculate the third values in the numerator for each individual; the l-th
  # value is for individual l
  num_term3 <- rowSums(alpha_mat * (S_sq + mu_sq)) / sigmabeta_sq

  # calculate the denominator for each individual; l-th value is for individual l
  denom <- rowSums(alpha_mat) + n

  nt_l <- rep(NA, n)
  # iterate over the individuals to update each error variance
  for (l in 1:n){

    # sigmasq update for individual l:

    # calculate weighted version of X
    weights <- sqrt(D[ , l])
    X_w <- X_mat * matrix(weights, n, p)

    # diagonal elements of X transpose X weighted
    XtX_w <- diag(t(X_w) %*% X_w)

    # second numerator term
    num_term2_l <- t(XtX_w) %*% num_term2[l , ]
    nt_l[l] <- num_term2_l
    # apply update
    sigmasq[l] <- (resid_w[l] + num_term2_l + num_term3[l]) / denom[l]
  }

  # terms for sigmabeta_sq update

  # numerator is num_term3 without the division by sigmabeta_sq
  num <- num_term3 * sigmabeta_sq


  den_old <- denom
  # denominator is denom without summing n, and also scaled by sigmasq
  denom <- sigmasq * (denom - n)

  # update the slab variance
  sigmabeta_sq <- num / denom

  return(list(sigmasq = sigmasq, sigmabeta_sq = sigmabeta_sq,
              t1 = mean(sort(abs(resid_w), T)[1:20]),
              t2 = mean(sort(abs(nt_l), T)[1:20]),
              t3 = mean(sort(abs(num_term3), T)[1:20]),
              den = mean(sort(abs(den_old), T)[1:20])))
}

# prep the inputs for cavi_R
#save(input, file = "dev/analyses_demos_experiments/hyperparameter_specification/blowup_input.Rda")

load("dev/analyses_demos_experiments/hyperparameter_specification/blowup_input.Rda")
n <- nrow(input$X); p <- ncol(input$X) - 1
y <- input$X[ , 15]
D <- input$D
X_mat <- input$X[ , -15]
mu_mat <- matrix(0, n, p)
alpha_mat <- matrix(0.2, n, p)
sigmasq <- rep(var(y), n)
update_sigma <- update_sigmabetasq <- T
sigmabeta_sq <- rep(1, n)
pi_est <- 0.45
tolerance <- 1e-12
max_iter <- 1e3

# blow up
out1 <- cavi0(y, D, X_mat, mu_mat, alpha_mat, sigmasq, update_sigma,
              sigmabeta_sq, update_sigmabetasq, pi_est, tolerance, max_iter)

library(ggplot2)
library(latex2exp)
burn <- 10

# visualize sigma versus mu
sig_mu <- reshape2::melt(data.frame(index = 1:burn, mu = scale(out1$mmt[1:burn]),
                                    sigmasq = scale(out1$mst[1:burn])),
                         measure.vars = c("mu", "sigmasq"),
                         variable.name = "parameter", value.name = "z-score")

ggplot(sig_mu, aes(index, `z-score`, color = parameter)) +
  geom_line() + theme_bw() +
  ggtitle(TeX(paste("Values of $\\sigma^2$ and $\\mu$, first ten iterations, p + 1 = 25"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(labels = unname(TeX(c("$\\mu$", "$\\sigma^2$"))))

# not blow up
out2 <- cavi0(input$X[ , 1], D, input$X[ , 2:5], mu_mat[ , 1:4],
              alpha_mat[ , 1:4], sigmasq, update_sigma, sigmabeta_sq,
              update_sigmabetasq, pi_est, tolerance, max_iter)

# visualize sigma versus mu
sig_mu <- reshape2::melt(data.frame(index = 1:burn,
                                    mu = scale(out2$mmt[1:burn]),
                                    sigmasq = scale(out2$mst[1:burn])),
                         measure.vars = c("mu", "sigmasq"),
                         variable.name = "parameter", value.name = "z-score")

ggplot(sig_mu, aes(index, `z-score`, color = parameter)) +
  geom_line() + theme_bw() +
  ggtitle(TeX(paste("Values of $\\sigma^2$ and $\\mu$, first ten iterations, p + 1 = 5"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(labels = unname(TeX(c("$\\mu$", "$\\sigma^2$"))))

sig_sigterms <- reshape2::melt(data.frame(index = 1:burn,
                                          sigmasq = scale(out1$mst[1:burn]),
                                          term1 = scale(out1$t1[1:burn]),
                                          term2 = scale(out1$t2[1:burn]),
                                          term3 = scale(out1$t3[1:burn]),
                                          denom = scale(out1$den[1:burn])),
                               id.vars = "index", variable.name = "parameter",
                               value.name = "z-score")
sig_sigterms_un <- reshape2::melt(data.frame(index = 1:burn,
                                             sigmasq = out1$mst[1:burn],
                                             term1 = out1$t1[1:burn],
                                             term2 = out1$t2[1:burn],
                                             term3 = out1$t3[1:burn],
                                             denom = out1$den[1:burn]),
                                  id.vars = "index", variable.name = "parameter")

ggplot(sig_sigterms_un, aes(index, value, color = parameter)) +
  geom_line() + theme_bw() +
  ggtitle(TeX(paste("Values of $\\sigma^2$ and terms in $\\sigma^2$ update, first ten iterations, p + 1 = 25"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(labels = unname(TeX(c("$\\sigma^2$", "term 1", "term 2", "term 3", "denom"))))

ggplot(sig_sigterms, aes(index, `z-score`, color = parameter)) +
  geom_line() + theme_bw() +
  ggtitle(TeX(paste("Values of $\\sigma^2$ and terms in $\\sigma^2$ update, first ten iterations, p + 1 = 25"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(labels = unname(TeX(c("$\\sigma^2$", "term 1", "term 2", "term 3", "denom"))))

ggplot(sig_sigterms[sig_sigterms$parameter %in% c("sigmasq", "term1"), ],
       aes(index, `z-score`, color = parameter)) +
  geom_line() + theme_bw() +
  ggtitle(TeX(paste("Values of $\\sigma^2$ and terms in $\\sigma^2$ update, first ten iterations, p + 1 = 25"))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_discrete(labels = unname(TeX(c("$\\sigma^2$", "term 1"))))
