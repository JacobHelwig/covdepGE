## The core function that calculates the variational parameter updates and returns the final variational estimates
## for a single regression. So in the arguments, y plays the role of the response, i.e the j th variable whose
## conditional distribution given the remaining variables is being calculated.

## This calls the ELBO_calculator function

## X: the data matrix except the j th variable
## XtX is x transpose times x.
## DXtX:diagonal elements of XtX
## Diff_mat: XtX-diag(DXtX)
## Xty: x transpose times y
## sigmasq: variance of response given the parameters (homoscedastic part, actual variance sigma_sq/w_i)
## sigmabeta_sq: prior variance of coefficient parameter
## true_pi: estimate of spike and slab mixture proportion.


cov_vsvb <- function(y, X, Z, XtX, DXtX, Diff_mat, Xty, sigmasq, sigmabeta_sq, true_pi) {
  thres <- 1e-7
  tol <- 1e-9

  msg <- function(s, ...) {
    time <- format(Sys.time(), "%X")
    cat(sprintf("%s %s\n", time, s))
  }


  change_alpha <- rep(0.001, n * p) # alpha_new - alpha_int

  max_iter <- 100
  iter <- 1
  Mu_vec <- matrix(rep(mu, n), n * p, 1)
  while (sqrt(sum(change_alpha^2)) > tol & iter < max_iter) { # The max_iter controls the max number of iterations until convergence

    alpha_int <- alpha ## Initialization of inclusion probability parameter.


    alpha_mat <- matrix(alpha, n, p, byrow = TRUE)

    alpha_vec <- matrix(alpha, n * p, 1, byrow = TRUE)
    for (i in 1:n) {
      S_sq[i, ] <- sigmasq * (t(DXtX_Big_ind) %*% D_long[, i] + 1 / sigmabeta_sq)^(-1) ## variance parameter
    }

    S_sq_vec <- matrix(t(S_sq), n * p, 1)

    for (i in 1:n) {
      y_XW <- y_long_vec * X_vec * D_long[, i]
      y_XW_mat <- matrix(y_XW, n, p, byrow = TRUE)

      X_mu_alpha <- X_vec * Mu_vec * alpha_vec
      xmualpha_mat <- t(matrix(X_mu_alpha, p, n)) %*% (matrix(1, p, p) - diag(rep(1, p)))
      XW_mat <- matrix(X_vec * D_long[, i], n, p, byrow = TRUE) * xmualpha_mat






      mu_mat[i, ] <- (t(y_XW_mat) %*% matrix(1, n, 1) - (t(XW_mat) %*% matrix(1, n, 1))) * (S_sq[i, ] / sigmasq) ### ### CAVI updation of mean variational parameter mu
    }
    Mu_vec <- matrix(t(mu_mat), n * p, 1)

    vec_1 <- log(true_pi / (1 - true_pi)) ## term 1 of the update of inclusion probability alpha

    vec_2 <- as.matrix(0.5 * log(S_sq_vec / (sigmasq * sigmabeta_sq))) ## term 2 of update of alpha
    vec_3 <- as.matrix(Mu_vec^2 / (2 * S_sq_vec)) ## term 3 of update of alpha
    # (vec_1+vec_2+vec_3)[1:p]
    unlogitalpha <- vec_1 + vec_2 + vec_3 # Sum of 3 terms for alpha update
    # thres=10^{-9}
    lthres <- logit(thres)
    uthres <- logit(1 - thres)
    indlarge <- which(unlogitalpha > uthres)
    indsmall <- which(unlogitalpha < lthres)
    unlogitalpha[indlarge] <- uthres
    unlogitalpha[indsmall] <- lthres
    alpha[which(unlogitalpha > 9)] <- 1 # thresholding very large values to 1 for computational stability
    alpha[which(unlogitalpha <= 9)] <- 1 / (1 + exp(-unlogitalpha[which(unlogitalpha <= 9)])) ### ### CAVI updation of variational parameter alpha


    e <- 0
    for (i in 1:n) { ## calculates ELBO for the j th variable by adding the contribution of the parameter
      ## corresponding to every individual in study. i th iteration takes the contribution of the variational
      ## parameters corresponding to  the i th individual in study, but the information is borrowed from
      ## all the n individuals depending on the weights coded in D[,i]
      e <- e + ELBO_calculator(y, X_mat, S_sq[i, ], mu_mat[i, ], alpha_mat[i, ], sigmasq, sigmabeta_sq, true_pi, D[, i], n, p)
    }
    ELBO_LB <- e


    alpha_new <- alpha
    change_alpha <- alpha_new - alpha_int



    ELBO_LBit[iter] <- ELBO_LB
    iter <- iter + 1
  }
  ELBO_LBit <- ELBO_LBit[1:(iter - 1)]
  list(var.alpha = alpha, var.mu = mu_mat, var.S_sq = S_sq, var.elbo = ELBO_LB, var.elboit = ELBO_LBit)
}
