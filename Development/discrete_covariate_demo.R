set.seed(1)

# rm(list=ls())
rm(list = ls())
source("cov_vsvb.R")
source("ELBO_calculator.R")
library(ggplot2)
library(reshape2)
library(MASS)
library(varbvs)
library(ks)

logit <- function(x) {
  if ((x == 0) | (x == 1)) {
    return(0)
  } else {
    return(log(x / (1 - x)))
  }
}

# Data generation 
n <- 100
p <- 10

# generating the precision matrix: Assume two discrete covariate levels
Lam1 <- matrix(0, p + 1, 1)
Lam2 <- Lam1
Lam1 <- c(3, 3, 3, 3, rep(0, p - 3)) * 5 # For Z[i]=-0.1

# Same lambda for both covariate levels, corresponds to covariate 
Lam2 <- Lam1 


# independent levels
#corresponds to covariate dependent model, uncomment to try this out.
# Lam2=c(3,3,3,3,3,3,3,0,0,0,0)*5#For Z[i]= 0.1
# Lam2=c(rep(0,p-3),3,3,3,3)*5#For Z[i]= 0.1 

# covariance matrix for covariate level 1
Var1 <- solve(Lam1 %*% t(Lam1) + diag(rep(10, p + 1))) 

# covariance matrix for covariate level 2
Var2 <- solve(Lam2 %*% t(Lam2) + diag(rep(10, p + 1))) 

X1 <- mvrnorm(n / 2, rep(0, p + 1), Var1)
X2 <- mvrnorm(n / 2, rep(0, p + 1), Var2)

# Generating the covariates

# Initializing the covariate matrix
Z <- matrix(-1, n, p) 

# covariate creation
for (i in 1:n) {
  for (j in 1:p) {
    #Uncomment to lay with other covariate values
    # Z[i,j]=rnorm(1,-2,.1)*(i<50) +rnorm(1,2,0.1)*(i>=50) 
    Z[i, j] <- -.1 * (i <= n / 2) + .1 * (i > n / 2)
  }
}

# D is an n by n matrix of weights 
D <- matrix(1, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    D[i, j] <- dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, .1)
  }
}
for (i in 1:n) {
  # Scale weights to sum to n
  D[, i] <- n * (D[, i] / sum(D[, i])) 
  # D[,i]=1 # When there is no covariate information, set the weights to 
  # be 1 throughout.
}

# the i-th column of D_long is the i-th column of D with the elements repeated
# p times
D_long <- matrix(rep(D, each = p), n*p)

beta <- matrix(0, n, p) # Ground truth of the dependence structure
resp_index <- 1 # The index we consider as response

# The variable specific inclusion probability matrix: 
# i-th row corresponds to the dependence structure for the i-th subject, 
# j-th matrix corresponds to the j th variable as response and the remaining as 
# predictors.
mylist <- rep(list(beta), p + 1) 

data_mat <- rbind(X1, X2)

Adj_Mat_vb <- array(0, dim = c(p + 1, p + 1))

# big ind matrix is a matrix of p stacked I_p identities  
Big_ind <- matrix(rep(diag(p), p), n * p, p, T)

# the big loop
for (resp_index in 1:(p + 1)) {
  # for (i in 1:n) {
  #   beta[i, ] <- (t(Lam1[-resp_index]) > 0) * (i <= n / 2) + 
  #     (t(Lam2[-resp_index]) > 0) * (i > n / 2) # Ground truth
  # 
  #   for (j in 1:p) {
  #     #Uncomment to lay with other covariate values
  #     # Z[i,j]=rnorm(1,-2,.1)*(i<50) +rnorm(1,2,0.1)*(i>=50) 
  #     Z[i, j] <- -.1 * (i <= n / 2) + .1 * (i > n / 2)
  #   }
  # }

  # Set variable number `resp_index` as the response
  y <- data_mat[, resp_index] 
  
  # Set the remaining p variables as predictor.
  X_mat <- data_mat[, -resp_index] 
  X_vec <- matrix(0, n * p, 1)
  # X<- matrix(rep(0,n^2*p),nrow=n,ncol=n*p)

  X <- matrix(0, nrow = n, ncol = n * p)
  
  # X_vec is a vector of length n*p that is the rows of X_mat "unravelled" by row;
  # that is, the first element of X_vec is the 1,1 of X_mat; the second element
  # is the 1,2; third is 1,3, ect.
  
  # X is a n by n*p matrix; it consists of rbinding n n by p matrices together
  # the j-th matrix is the j-th row of X_mat in the j-th row, and 0's o.w.
  for (i in 1:n) {
    for (j in 1:p) {
      k <- p * (i - 1) + 1
      X[i, k + j - 1] <- X_mat[i, j]
      X_vec[k + j - 1] <- X[i, k + j - 1]
    }
  }
  ELBO_LBit <- rep(0, 10000)
  
  # big diag mat looks exactly the same as X, except it replaces all non-zero
  # entries with 1
  Big_diag_mat <- (X != 0) * 1
  # Big_diag_mat <- matrix(0, nrow = n, ncol = n * p)
  # for (i in 1:n) {
  #   k <- p * (i - 1)
  #   for (j in 1:p) {
  #     Big_diag_mat[i, k + j] <- 1
  #   }
  # }

  q <- matrix(2, n, 1)

  sigmasq <- 1 # Initialization of the hyperparameter value
  E <- rnorm(n, 0, sigmasq)

  # a block diagonal matrix; the j-th block is the transpose of the j-th row of 
  # X times the j-th row of X; it is an n*p by n*p matrix
  XtX <- t(X) %*% X
  
  # the diagonal values of XtX
  DXtX <- diag(XtX)
  DXtX_rep <- rep(DXtX, p)
  DXtX_mat <- matrix(DXtX_rep, n * p, p)
  
  # XtX with its diagonal removed and replaced with 0
  Diff_mat <- XtX - diag(DXtX)

  # # D is a matrix of weights
  # D <- matrix(1, n, n)
  # for (i in 1:n) {
  #   for (j in 1:n) {
  #     D[i, j] <- dnorm(norm(Z[i, ] - Z[j, ], "2"), 0, .1)
  #   }
  # }
  # for (i in 1:n) {
  #   # Scale weights to sum to n
  #   D[, i] <- n * (D[, i] / sum(D[, i])) 
  #   # D[,i]=1 # When there is no covariate information, set the weights to 
  #   # be 1 throughout.
  # }

  # Initialization of the inclusion probability matrix for a fixed variable
  # with i-th row corresponding to i th subject.
  alpha <- rep(0.2, n * p) 
  
  sigmabeta_sq <- 3 # Initialization for hyperparameter
  mu <- rep(0, p) # Variational parameter
  true_pi <- 0.5 # Hyperparameter
  
  # y_long_vec is each element of y repeated p times 
  #y_long_vec <- as.vector(t(y %*% matrix(1, 1, p)))
  y_long_vec <- rep(y, each = p)
  
  Xty <- t(X) %*% y
  beta_mat <- matrix(beta, n, p, byrow = TRUE)
  mu_mat <- beta_mat

  # # the i-th column of D_long is the i-th column of D with the elements repeated
  # # p times
  # D_long <- matrix(rep(D, each = p), n*p)
  # D_long <- matrix(0, n * p, n)
  # for (i in 1:n) {
  #   D_long[ , i] <- rep(D[ , i], each = p)
  #   #D_long[, i] <- matrix(t(D[, i] %*% matrix(1, 1, p)), n * p, 1)
  # }

  S_sq <- matrix(sigmasq * (DXtX + 1 / sigmabeta_sq)^(-1), n, p)

  iter <- 1

  # ind_vec <- seq(0, (n - 1) * p, by = p)
  # Ind_mat <- matrix(0, n, p)
  # for (j in 1:p) {
  #   Ind_mat[, j] <- ind_vec + j
  # }
  # Ind_mat <- matrix(1:(n * p), n, p, T)
  # Big_ind <- matrix(rep(diag(p), p), n * p, p, T)
  # Big_ind <- matrix(0, n * p, p)
  # Big_ind_1 <- matrix(0, n * p, p)
  
   
  # for (j in 1:p) {
  #   Big_ind[Ind_mat[, j], j] <- 1
  #   Big_ind_1[Ind_mat[, j], j] <- 0
  # }

  DXtX_Big_ind <- DXtX_mat * Big_ind
  
  candL <- seq(0.1, 0.9, .2) # Different values of hyperparameter true_pi
  # candL=0.5
  elb <- rep(0, length(candL))

  est_q <- rep(0, n)
  beta_matr <- matrix(0, n, p)

  #################### tuning hyperparameters ##################################
  
  # Setting hyperparameter value as in Carbonetto Stephens model
  idmod <- varbvs(X_mat, y, Z = Z[, 1], verbose = FALSE) 
  inprob <- idmod$pip
  rest_index_set <- setdiff(c(1:(p + 1)), resp_index)

  sigmasq <- mean(idmod$sigma)
  pi_est <- mean(1 / (1 + exp(-idmod$logodds)))
  sigmavec <- c(0.01, 0.05, 0.1, 0.5, 1, 3, 7, 10)
  
  # vector for storing the ELBO for each value of sigma in sigmavec
  elb1 <- matrix(0, length(sigmavec), 1)
  
  # loop to optimize sigma
  for (j in 1:length(sigmavec)) {
    res <- cov_vsvb(y, X, Z, XtX, DXtX, Diff_mat, Xty, sigmasq, sigmavec[j], pi_est)
    elb1[j] <- res$var.elbo
  }
  
  # Select the value of sigma_beta that maximizes the elbo
  sigmabeta_sq <- sigmavec[which.max(elb1)] 
  
  # fit another model using this value of sigma_beta
  result <- cov_vsvb(y, X, Z, XtX, DXtX, Diff_mat, Xty, sigmasq, sigmabeta_sq, pi_est)
  
  # vector of length n * p of inclusion probabilities
  incl_prob <- result$var.alpha
  
  # n by p matrix; the i,j-th entry is the probability of inclusion for the 
  # i-th individual for the j-th variable according to the regression on y
  heat_alpha <- matrix(incl_prob, n, p, byrow = TRUE)
  mylist[[resp_index]] <- heat_alpha
}

mylist2 <- mylist
load("original_discrete_alpha_matrices.Rdata")
same <- T
for (j in 1:length(mylist)){
  if (all.equal(mylist[[j]], mylist2[[j]]) != T){
    same <- F
    break
  }
}
same

# beyond here is demonstrations

# alph <- matrix(0, p + 1, p + 1)
# 
# alph <- matrix(0, p + 1, p + 1)
# 
# SUBJECT <- 1
# for (i in 1:(p + 1)) {
#   alph[i, -i] <- mylist[[i]][SUBJECT, ] # Individual specific inclusion probability matrix
# }
# beta <- matrix(0, p + 1, p + 1)
# for (i in 1:(p + 1)) {
#   for (j in 1:(p + 1)) {
#     beta[i, j] <- (Lam1[i] != 0 & Lam1[j] != 0)
#   }
# }
# diag(beta) <- 0
# 
# heat_alpha <- alph
# 
# a <- heat_alpha
# for (i in 1:(p + 1)) {
#   for (j in i:(p + 1)) {
#     #  a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
#     a[i, j] <- mean(c(heat_alpha[i, j], heat_alpha[j, i]))
#     a[j, i] <- a[i, j]
#   }
# }
# 
# heat_alpha <- a
# 
# 
# 
# alphvec <- sort(as.vector(heat_alpha[which(heat_alpha != 0)]))
# 
# selection1 <- 1 * (heat_alpha > 0.5)
# 
# 
# 
# SUBJECT <- 100
# 
# for (i in 1:(p + 1)) {
#   alph[i, -i] <- mylist[[i]][SUBJECT, ]
# }
# beta <- matrix(0, p + 1, p + 1)
# for (i in 1:(p + 1)) {
#   for (j in 1:(p + 1)) {
#     beta[i, j] <- (Lam2[i] != 0 & Lam2[j] != 0)
#   }
# }
# diag(beta) <- 0
# 
# heat_alpha <- alph
# 
# 
# a <- heat_alpha
# for (i in 1:(p + 1)) {
#   for (j in i:(p + 1)) {
#     #  a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
#     a[i, j] <- mean(c(heat_alpha[i, j], heat_alpha[j, i]))
#     a[j, i] <- a[i, j]
#   }
# }
# 
# heat_alpha <- a
# 
# 
# alphvec <- sort(as.vector(heat_alpha[which(heat_alpha != 0)]))
# 
# selection1 <- 1 * (heat_alpha > 0.5)
# 
# 
# 
# SUBJECT <- 1
# for (i in 1:(p + 1)) {
#   alph[i, -i] <- mylist[[i]][SUBJECT, ]
# }
# beta <- matrix(0, p + 1, p + 1)
# for (i in 1:(p + 1)) {
#   for (j in 1:(p + 1)) {
#     beta[i, j] <- (Lam1[i] != 0 & Lam1[j] != 0)
#   }
# }
# diag(beta) <- 0
# 
# data <- melt(t(beta))
# fig <- ggplot(data, aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile(color = "brown") +
#   scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0), guide = "legend")
# 
# 
# 
# 
# fig <- fig + scale_x_continuous(expand = c(0, 0))
# fig <- fig + scale_y_continuous(expand = c(0, 0))
# 
# 
# fig <- fig + labs(x = expression(bold(Variables)), y = expression(bold(Variables)), title = expression(bold("True Dependence")))
# 
# fig <- fig + theme(
#   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#   panel.background = element_blank()
# )
# fig <- fig + theme(plot.title = element_text(hjust = 0.5))
# 
# 
# fig <- fig + theme(
#   axis.text = element_text(size = 15, face = "bold", colour = "black"),
#   axis.title = element_text(size = 30, face = "bold")
# )
# 
# fig <- fig + theme(
#   legend.title = element_text(face = "bold", size = 25),
#   legend.text = element_text(face = "bold", size = 25),
#   legend.key.size = unit(2, "lines")
# )
# fig <- fig + coord_equal()
# 
# fig
# 
# 
# heat_alpha <- alph
# 
# 
# a <- heat_alpha
# for (i in 1:(p + 1)) {
#   for (j in i:(p + 1)) {
#     #  a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
#     a[i, j] <- 0.5 * (heat_alpha[i, j] + heat_alpha[j, i])
#     a[j, i] <- a[i, j]
#   }
# }
# 
# # heat_alpha=a
# 
# data <- melt(t(a))
# fig <- ggplot(data, aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile(color = "brown") +
#   scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0), guide = "colorbar")
# 
# fig <- fig + scale_x_continuous(expand = c(0, 0))
# fig <- fig + scale_y_continuous(expand = c(0, 0))
# 
# 
# fig <- fig + labs(x = expression(bold(Variables)), y = expression(bold(Variables)), title = expression(bold("Inclusion Probability for Covariate Level 1")))
# 
# fig <- fig + theme(
#   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#   panel.background = element_blank()
# )
# fig <- fig + theme(plot.title = element_text(hjust = 0.5))
# 
# 
# fig <- fig + theme(
#   axis.text = element_text(size = 15, face = "bold", colour = "black"),
#   axis.title = element_text(size = 30, face = "bold")
# )
# 
# fig <- fig + theme(
#   legend.title = element_text(face = "bold", size = 25),
#   legend.text = element_text(face = "bold", size = 25),
#   legend.key.size = unit(2, "lines")
# )
# fig <- fig + coord_equal()
# 
# fig
# 
# alphvec <- setdiff(as.vector(heat_alpha), diag(heat_alpha))
# 
# selection1 <- 1 * (a > 0.5)
# data <- melt(t(selection1))
# fig <- ggplot(data, aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile(color = "brown") +
#   scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0), guide = "legend")
# 
# fig <- fig + scale_x_continuous(expand = c(0, 0))
# fig <- fig + scale_y_continuous(expand = c(0, 0))
# 
# 
# fig <- fig + labs(x = expression(bold(Variables)), y = expression(bold(Variables)), title = expression(bold("Graph Estimate For Covariate Level 1")))
# 
# fig <- fig + theme(
#   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#   panel.background = element_blank()
# )
# fig <- fig + theme(plot.title = element_text(hjust = 0.5))
# 
# 
# fig <- fig + theme(
#   axis.text = element_text(size = 15, face = "bold", colour = "black"),
#   axis.title = element_text(size = 30, face = "bold")
# )
# 
# fig <- fig + theme(
#   legend.title = element_text(face = "bold", size = 25),
#   legend.text = element_text(face = "bold", size = 25),
#   legend.key.size = unit(2, "lines")
# )
# fig <- fig + coord_equal()
# 
# fig
# 
# SUBJECT <- 100
# 
# for (i in 1:(p + 1)) {
#   alph[i, -i] <- mylist[[i]][SUBJECT, ]
# }
# beta <- matrix(0, p + 1, p + 1)
# for (i in 1:(p + 1)) {
#   for (j in 1:(p + 1)) {
#     beta[i, j] <- (Lam2[i] != 0 & Lam2[j] != 0)
#   }
# }
# diag(beta) <- 0
# data <- melt(t(beta))
# fig <- ggplot(data, aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile(color = "brown") +
#   scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0), guide = "legend")
# 
# fig <- fig + scale_x_continuous(expand = c(0, 0))
# fig <- fig + scale_y_continuous(expand = c(0, 0))
# 
# 
# fig <- fig + labs(x = expression(bold(Variables)), y = expression(bold(Variables)), title = expression(bold("True Dependence")))
# 
# fig <- fig + theme(
#   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#   panel.background = element_blank()
# )
# fig <- fig + theme(plot.title = element_text(hjust = 0.5))
# 
# 
# fig <- fig + theme(
#   axis.text = element_text(size = 15, face = "bold", colour = "black"),
#   axis.title = element_text(size = 30, face = "bold")
# )
# 
# fig <- fig + theme(
#   legend.title = element_text(face = "bold", size = 25),
#   legend.text = element_text(face = "bold", size = 25),
#   legend.key.size = unit(2, "lines")
# )
# fig <- fig + coord_equal()
# 
# fig
# 
# heat_alpha <- alph
# 
# 
# a <- heat_alpha
# for (i in 1:(p + 1)) {
#   for (j in i:(p + 1)) {
#     #  a[i,j]=max(heat_alpha[i,j],heat_alpha[j,i])
#     a[i, j] <- 0.5 * (heat_alpha[i, j] + heat_alpha[j, i])
#     a[j, i] <- a[i, j]
#   }
# }
# 
# # heat_alpha=a
# data <- melt(t(a))
# fig <- ggplot(data, aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile(color = "brown") +
#   scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0), guide = "colorbar")
# 
# fig <- fig + scale_x_continuous(expand = c(0, 0))
# fig <- fig + scale_y_continuous(expand = c(0, 0))
# 
# 
# fig <- fig + labs(x = expression(bold(Variables)), y = expression(bold(Variables)), title = expression(bold("Inclusion Probability For Covariate Level 2")))
# 
# fig <- fig + theme(
#   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#   panel.background = element_blank()
# )
# fig <- fig + theme(plot.title = element_text(hjust = 0.5))
# 
# 
# fig <- fig + theme(
#   axis.text = element_text(size = 15, face = "bold", colour = "black"),
#   axis.title = element_text(size = 30, face = "bold")
# )
# 
# fig <- fig + theme(
#   legend.title = element_text(face = "bold", size = 25),
#   legend.text = element_text(face = "bold", size = 25),
#   legend.key.size = unit(2, "lines")
# )
# fig <- fig + coord_equal()
# fig
# 
# alphvec <- setdiff(as.vector(heat_alpha), diag(heat_alpha))
# 
# 
# selection1 <- 1 * (a > 0.5)
# data <- melt(t(selection1))
# fig <- ggplot(data, aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile(color = "brown") +
#   scale_fill_gradient(low = "white", high = "steelblue", breaks = c(1, 0), guide = "legend")
# 
# fig <- fig + scale_x_continuous(expand = c(0, 0))
# fig <- fig + scale_y_continuous(expand = c(0, 0))
# 
# 
# fig <- fig + labs(x = expression(bold(Variables)), y = expression(bold(Variables)), title = expression(bold("Graph Estimate For Covariate Level 2")))
# 
# fig <- fig + theme(
#   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#   panel.background = element_blank()
# )
# 
# fig <- fig + coord_equal()
# fig <- fig + theme(
#   axis.text = element_text(size = 15, face = "bold", colour = "black"),
#   axis.title = element_text(size = 30, face = "bold")
# )
# fig <- fig + theme(plot.title = element_text(hjust = 0.5))
# fig <- fig + theme(
#   legend.title = element_text(face = "bold", size = 25),
#   legend.text = element_text(face = "bold", size = 25),
#   legend.key.size = unit(2, "lines")
# )
# 
# fig
