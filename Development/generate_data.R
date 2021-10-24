#-------------------------------------------------------------------------------
#-------------------FUNCTION TO GENERATE CONTINUOUS DATA------------------------
#-------------------------------------------------------------------------------

library(MASS)

generate_continuous <- function(seed = 1){

  set.seed(seed)

  # Data generation
  n <- 180
  p <- 4

  # function to create sigma matrix for a p-dimensional gaussian given the value
  # of an extraneous covariate
  Var_cont <- function(z) {
    STR <- 1
    pr <- matrix(0, p + 1, p + 1)
    diag(pr) <- 2
    pr[2, 3] <- STR
    pr[1, 2] <- STR * ((z > -1) && (z < -.33)) + (STR - STR * ((z + .23) / .56)) * ((z > -0.23) && (z < 0.33)) + (0) * ((z > 0.43) && (z < 1))
    pr[1, 3] <- 0 * ((z > -1) && (z < -.33)) + (STR * ((z + .23) / .56)) * ((z > -0.23) && (z < 0.33)) + (STR) * ((z > 0.43) && (z < 1))
    pr[2, 1] <- pr[1, 2]
    pr[3, 1] <- pr[1, 3]
    pr[3, 2] <- pr[2, 3]
    Var <- solve(pr)
    return(Var)
  }

  # creation of covariate
  Z <- c(
    seq(-0.99, -0.331, (-.331 + .99) / 59),
    seq(-0.229, 0.329, (.329 + .229) / 59),
    seq(0.431, .99, (.99 - .431) / 59)
  )
  Z <- matrix(Z, n, 1)

  # creation the data matrix; each individual is generated from a MVN with 0 mean
  # and covariance matrix determined by their corresponding extraneous covariate
  data_mat <- matrix(0, n, p + 1)
  for (i in 1:n) {
    data_mat[i, ] <- mvrnorm(1, rep(0, p + 1), Var_cont(Z[i]))
  }

  return(list(data = data_mat, covts = Z))

}

#-------------------------------------------------------------------------------
#-------------------FUNCTION TO GENERATE DISCRETE DATA--------------------------
#-------------------------------------------------------------------------------

generate_discrete <- function(seed = 1){

  set.seed(seed)

  # Data generation
  n <- 100
  p <- 10

  # generating the precision matrix: Assume two discrete covariate levels
  Lam1 <- c(3, 3, 3, 3, rep(0, p - 3)) * 5 # For Z[i]=-0.1

  # Same lambda for both covariate levels, corresponds to covariate
  Lam2 <- Lam1

  # covariance matrix for covariate level 1
  Var1 <- solve(Lam1 %*% t(Lam1) + diag(rep(10, p + 1)))

  # covariance matrix for covariate level 2
  Var2 <- solve(Lam2 %*% t(Lam2) + diag(rep(10, p + 1)))

  # Initializing the covariate matrix
  Z <- matrix(-1, n, p)

  # covariate creation; half the individuals get a 0.1, while the others get -0.1
  for (i in 1:n) {
    for (j in 1:p) {
      Z[i, j] <- -.1 * (i <= n / 2) + .1 * (i > n / 2)
    }
  }

  # create the data matrix; half of the individuals are generated from a MVN with
  # 0 mean vector and covariance matrix corresponding to covariate level 1, while
  # the other are from an MVN with covariance matrix corresponding to covariate
  # level 2
  X1 <- mvrnorm(n / 2, rep(0, p + 1), Var1)
  X2 <- mvrnorm(n / 2, rep(0, p + 1), Var2)
  data_mat <- rbind(X1, X2)

  return(list(data = data_mat, covts = Z))

}
