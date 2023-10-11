# function for generating continuous covariate dependent data
cont_cov_dep_data <- function(p, n1, n2, n3){

  # create covariate for observations in each of the three intervals

  # define number of samples
  n <- sum(n1, n2, n3)

  # define the intervals
  limits1 <- c(-3, -1)
  limits2 <- c(-1, 1)
  limits3 <- c(1, 3)

  # define the interval labels
  interval <- c(rep(1, n1), rep(2, n2), rep(3, n3))

  # draw the covariate values within each interval
  z1 <- sort(stats::runif(n1, limits1[1], limits1[2]))
  z2 <- sort(stats::runif(n2, limits2[1], limits2[2]))
  z3 <- sort(stats::runif(n3, limits3[1], limits3[2]))
  Z <- matrix(c(z1, z2, z3), n, 1)

  # the shared part of the structure for all three intervals is a 2 on the
  # diagonal and a 1 in the (2, 3) position
  common_str <- diag(p)
  common_str[2, 3] <- 1

  # define constants for the structure of interval 2
  beta1 <- diff(limits2)^-1
  beta0 <- -limits2[1] * beta1

  # ggplot() + geom_function(fun=function(x) (x < 1) * pmin(1, 1 - 0.5 - 0.5 * x), color='red',n=1000)+ geom_function(fun=function(x) (x > -1) * pmin(1, 0.5 + 0.5 * x), color='blue',n=100)+ xlim(-3,3)

  # define omega12 and omega 13
  omega12 <- (Z < 1) * pmin(1, 1 - beta0 - beta1 * Z)
  omega13 <- (Z > -1) * pmin(1, beta0 + beta1 * Z)

  # interval 2 has two different linear functions of Z in the (1, 2) position
  # and (1, 3) positions; define structures for each of these components
  str12 <- str13 <- matrix(0, p, p)
  str12[1, 2] <- str13[1, 3] <- 1

  # create the precision matrices
  prec_mats <- vector("list", n)
  for (j in 1:n){
    prec_mats[[j]] <- common_str + omega12[j] * str12 + omega13[j] * str13
  }

  # symmetrize the precision matrices
  true_precision <- lapply(prec_mats, function(mat) t(mat) + mat)

  # invert the precision matrices to get the covariance matrices
  cov_mats <- lapply(true_precision, solve)

  # generate the data using the covariance matrices
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p)))

  return(list(X = data_mat, Z = Z, true_precision = true_precision,
              interval = interval))
}

# function for generating continuous sinusoidal covariate dependent data
cont_cov_dep_sine_data <- function(p, n1, n2, n3){

  # create covariate for observations in each of the three intervals

  # define number of samples
  n <- sum(n1, n2, n3)

  # define the intervals
  limits1 <- c(-3, -2)
  limits2 <- c(-2, -1)
  limits3 <- c(-1, 1)

  # define the interval labels
  interval <- c(rep(1, n1), rep(2, n2), rep(3, n3))

  # draw the covariate values within each interval
  z1 <- stats::runif(n1, limits1[1], limits1[2])
  z1 <- c(z1[1:(floor(n1/2))], -z1[(floor(n1/2) + 1):n1])
  z2 <- stats::runif(n2, limits2[1], limits2[2])
  z2 <- c(z2[1:(ceiling(n2/2))], -z2[(ceiling(n2/2) + 1):n2])
  z3 <- stats::runif(n3, limits3[1], limits3[2])
  Z <- matrix(sort(c(z1, z2, z3)), n, 1)

  # the shared part of the structure for all three intervals is a 2 on the
  # diagonal and a 1 in the (2, 3) position
  common_str <- diag(p)
  common_str[2, 3] <- 1

  # ggplot() + geom_function(fun=function(x) pmax(cospi((x+3)/4),0)+pmax(cospi((x-3)/4),0), color='red',n=1000)+ geom_function(fun=function(x) pmax(0,cospi(x/4)), color='blue',n=100)+ xlim(-3,3)

  # define omega12 and omega 13
  omega12 <- pmax(cospi((Z+3)/4),0)+pmax(cospi((Z-3)/4),0)
  omega13 <- pmax(0,cospi(Z/4))

  # interval 2 has two different linear functions of Z in the (1, 2) position
  # and (1, 3) positions; define structures for each of these components
  str12 <- str13 <- matrix(0, p, p)
  str12[1, 2] <- str13[1, 3] <- 1

  # create the precision matrices
  prec_mats <- vector("list", n)
  for (j in 1:n){
    prec_mats[[j]] <- common_str + omega12[j] * str12 + omega13[j] * str13
  }

  # symmetrize the precision matrices
  true_precision <- lapply(prec_mats, function(mat) t(mat) + mat)

  # invert the precision matrices to get the covariance matrices
  cov_mats <- lapply(true_precision, solve)

  # generate the data using the covariance matrices
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p)))

  return(list(X = data_mat, Z = Z, true_precision = true_precision,
              interval = interval))
}

# function for generating multivariate continuous covariate dependent data
cont_multi_cov_dep_data <- function(p, n){

  # create covariate for observations in each of the three intervals

  # define the intervals
  limits1 <- c(-3, -1)
  limits2 <- c(-1, 1)
  limits3 <- c(1, 3)
  intervals <- list(limits1, limits2, limits3)

  # draw the covariate values within each interval
  Z <- matrix(NA, 0, 2)
  for (int_x in intervals){
    for (int_y in intervals){
      x <- runif(n, int_x[1], int_x[2])
      y <- runif(n, int_y[1], int_y[2])
      Z <- rbind(Z, cbind(x, y))
    }
  }

  # the shared part of the structure for all three intervals is a 2 on the
  # diagonal and a 1 in the (2, 3) position
  common_str <- diag(p)
  common_str[2, 3] <- 1

  # define constants for the structure of interval 2
  beta1 <- diff(limits2)^-1
  beta0 <- -limits2[1] * beta1

  # define omega12 and omega 13
  omega12 <- (Z[ , 1] < 1) * pmin(1, 1 - beta0 - beta1 * Z[ , 1])
  omega13 <- (Z[ , 2] > -1) * pmin(1, beta0 + beta1 * Z[ , 2])

  # interval 2 has two different linear functions of Z in the (1, 2) position
  # and (1, 3) positions; define structures for each of these components
  str12 <- str13 <- matrix(0, p, p)
  str12[1, 2] <- str13[1, 3] <- 1

  # create the precision matrices
  prec_mats <- vector("list", n)
  for (j in 1:(9 * n)){
    prec_mats[[j]] <- common_str + omega12[j] * str12 + omega13[j] * str13
  }

  # symmetrize the precision matrices
  true_precision <- lapply(prec_mats, function(mat) t(mat) + mat)

  # invert the precision matrices to get the covariance matrices
  cov_mats <- lapply(true_precision, solve)

  # generate the data using the covariance matrices
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p)))

  return(list(X = data_mat, Z = Z, true_precision = true_precision))
}

# function for generating 4D continuous covariate dependent data
cont_4_cov_dep_data <- function(p, n, Z=NULL){

    # create covariate for observations in each of the three intervals

    # # define the intervals
    # limits1 <- c(-3, -9/5)
    # limits2 <- c(-9/5, -3/5)
    # limits3 <- c(-3/5, 3/5)
    # limits4 <- c(3/5, 9/5)
    # limits5 <- c(9/5, 3)
    # intervals <- list(limits1,limits2, limits3,limits4,limits5)
    #
    # # draw the covariate values within each interval
    # Z <- matrix(NA, 0, 4)
    # for (int_x in intervals){
    #   for (int_y in intervals){
    #     for (int_z1 in intervals){
    #       for (int_z2 in intervals){
    #         x <- runif(n, int_x[1], int_x[2])
    #         y <- runif(n, int_y[1], int_y[2])
    #         z1 <- runif(n, int_z1[1], int_z1[2])
    #         z2 <- runif(n, int_z2[1], int_z2[2])
    #         Z <- rbind(Z, cbind(x, y, z1, z2))
    #       }
    #     }
    #   }
    # }

    # define the intervals
    limits <- c(-3, 3)

    # draw the covariate values
    if (is.null(Z)){
      Z <- matrix(stats::runif(4 * n, limits[1], limits[2]), n, 4)
    }else{
      n <- nrow(Z)
    }

    # the shared part of the structure is a 2 on the diagonal
    common_str <- diag(p)
    common_str[2, 3] <- common_str[5, 6] <- 1

    # define constants for the covariate-dependent structure
    beta1 <- 0.5

    # define precision value for the first 5 values dependent on the first 5
    # covariates, and the next 5 values dependent on the next 5 covariates
    # p <- 5; by <- 0.1; Z <- seq(-3, 3, by); n <- length(Z); omega <- matrix(Z, length(Z), 4);     common_str <- diag(p)
    omega <- Z

    # 1,2 and 1,3 entries
    omega[,1] <-  pmax(pmin(1, -beta1 * (omega[,1] - 9/5)), 0)
    omega[,2] <-  pmax(pmin(1, beta1 * (omega[,2] + 3/5)), 0)


    # 4,5 and 4,6 entries
    omega[,3] <-  pmax(pmin(1, -beta1 * (omega[,3] - 3/5)), 0)
    omega[,4] <-  pmax(pmin(1, beta1 * (omega[,4] + 9/5)), 0)

    # g <- ggplot()
    # for (i in 1:2){
    #   g <- local({
    #     j <- i
    #     g + geom_line(aes(Z, omega[,j]), col=j) + geom_line(aes(Z, omega[,j+2]), col=j+2)
    #     # g + geom_line(aes(Z, omega[,j+2]), col=j+2)
    #
    #   })
    # }
    # g + xlim(-3, 3) + scale_x_continuous(breaks=seq(-3, 3, 0.2))

    # ggplot() + geom_function(fun=function(x)pmax(pmin(-0.5*(x-0.6),1),0))+ geom_function(fun=function(x)pmax(pmin(-0.5*(x-1.8),1),0)) + geom_function(fun=function(x)pmax(pmin(0.5*(x+0.6),1),0))+ geom_function(fun=function(x)pmax(pmin(0.5*(x+1.8),1),0)) + xlim(-3,3)

    # create the precision matrices
    prec_mats <- vector("list", n)
    for (i in 1:n){
      prec_mats[[i]] <- common_str
      for (j in 1:4){
        shift <- 3 * (j > 2)
        prec_mats[[i]][1 + shift, j + 1 + (j > 2)] <- omega[i, j]
      }
    }

    # symmetrize the precision matrices
    true_precision <- lapply(prec_mats, function(mat) t(mat) + mat)
    # return(true_precision)
    # library(covdepGE)
    # for (i in 1:n){
    #   print(i)
    #   print(matViz(true_precision[[i]], incl_val = T) + ggtitle(i))
    # }

    # invert the precision matrices to get the covariance matrices
    cov_mats <- lapply(true_precision, solve)

    # generate the data using the covariance matrices
    data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p)))

    return(list(X = data_mat, Z = Z, true_precision = true_precision))
}


# extraneous covariate visualizations
if (F){
  rm(list=ls())
  library(kableExtra)
  library(ggplot2)
  library(ggpubr)
  library(extrafont)
  library(latex2exp)
  windowsFonts(Times=windowsFont("TT Times New Roman"))
  colors <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF")
  plots <- list(
  ggplot() + xlim(-3, 3) + ylim(0,1) +
    geom_function(fun=function(x)pmax(0, pmin(1, -0.5 * (x - 1))), n=1e4, aes(col="1"), size=1) +
    geom_function(fun=function(x)pmax(0, pmin(1, 0.5 * (x + 1))), n=1e4, aes(col="2"), size=1) +
    geom_function(fun=function(x)5, n=1e4, aes(col="3"), size=1) +
    geom_function(fun=function(x)5, n=1e4, aes(col="4"), size=1) +
    theme_pubclean() +
    theme(text = element_text(family = "Times", size = 18),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colors, labels = unname(TeX(c("$\\textit{l}=1$", "$\\textit{l}=2$", "$\\textit{l}=3$", "$\\textit{l}=4$"))), name = "") +
    labs(x = " ", y = TeX("$\\Omega(z_l)$")) + theme(legend.key = element_rect(fill = "white"), legend.position = "right") +
    ggtitle(TeX("$\\textit{\\Omega(z_l), q}\\in\\{1,2\\}$")),
  ggplot() + xlim(-3, 3) +
    geom_function(fun=function(x) pmax(pmin(1, -0.5 * (x - 9/5)), 0), n=1e4, aes(col="1"), size=1) +
    geom_function(fun=function(x) pmax(pmin(1, 0.5 * (x + 3/5)), 0), n=1e4, aes(col="2"), size=1) +
    geom_function(fun=function(x) pmax(pmin(1, -0.5 * (x - 3/5)), 0), n=1e4, aes(col="3"), size=1) +
    geom_function(fun=function(x) pmax(pmin(1, 0.5 * (x + 9/5)), 0), n=1e4, aes(col="4"), size=1) +
    theme_pubclean() +
    theme(text = element_text(family = "Times", size = 18),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colors, labels = unname(TeX(c("$\\textit{l}=1$", "$\\textit{l}=2$", "$\\textit{l}=3$", "$\\textit{l}=4$"))), name = "") +    labs(x = TeX("$\\textit{z_l}$"), y = TeX("$\\Omega(z_l)$")) + theme(legend.key = element_rect(fill = "white"), legend.position = "right") +
    ggtitle(TeX("$\\textit{\\Omega(z_l), q=4}$")),
  ggplot() + xlim(-3, 3) +
    geom_function(fun=function(x)pmax(cospi((x+3)/4),0)+pmax(cospi((x-3)/4),0), n=1e4, aes(col="1"), size=1) +
    geom_function(fun=function(x)pmax(0,cospi(x/4)), n=1e4, aes(col="2"), size=1) +
    theme_pubclean() +
    theme(text = element_text(family = "Times", size = 18),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colors, labels = unname(TeX(c("$\\textit{l}=1$", "$\\textit{l}=2$", "$\\textit{l}=3$", "$\\textit{l}=4$"))), name = "") +
    labs(x = " ", y = TeX("$\\Omega(z_l)$")) + theme(legend.key = element_rect(fill = "white"), legend.position = "right") +
    ggtitle(TeX("$\\textit{\\Omega(z_l), q=1}$ (non-linear)"))
  )
  plots <- lapply(plots, function(plot) plot + rremove("ylab"))
  fig <- ggarrange(plotlist = plots, nrow=1, common.legend = TRUE,  legend="bottom")
  fig <- annotate_figure(fig, left = text_grob(TeX("$\\textit{\\Omega(z_l)_{j,k}}$"), family="Times", size = 18, rot=90))
  ggsave("plots/omega_plots.pdf", fig, height = 4, width = 14)
}


# # function for generating 10D continuous covariate dependent data
# cont_4_cov_dep_data <- function(p, n=NULL, Z=NULL){
#
#   # create covariate for observations in each of the three intervals
#
#   # define the intervals
#   limits <- c(-3, 3)
#
#   # draw the covariate values
#   if (is.null(Z)){
#     Z <- matrix(stats::runif(4 * n, limits[1], limits[2]), n, 4)
#   }else{
#     n <- nrow(Z)
#   }
#
#   # the shared part of the structure is a 2 on the diagonal
#   common_str <- diag(p)
#
#   # define constants for the covariate-dependent structure
#   beta1 <- 0.5
#
#   # define precision value for the first 5 values dependent on the first 5
#   # covariates, and the next 5 values dependent on the next 5 covariates
#   # p <- 5; by <- 0.1; Z <- seq(-3, 3, by); n <- length(Z); omega <- matrix(Z, length(Z), 4);     common_str <- diag(p)
#   omega <- Z
#   for (i in 1:2){
#     shift <- -1.5 * (i - 1)
#     omega[,i] <-  pmax(pmin(1, -beta1 * (omega[,i] + shift)), 0)
#     j <- 5-i
#     omega[,j] <-  pmax(pmin(1, beta1 * (omega[,j] - shift)), 0)
#   }
#   # g <- ggplot()
#   # for (i in 1:2){
#   #   g <- local({
#   #     j <- i
#   #     g + geom_line(aes(Z, omega[,j]), col=j) + geom_line(aes(Z, omega[,j+2]), col=j+2)
#   #     # g + geom_line(aes(Z, omega[,j+2]), col=j+2)
#   #
#   #   })
#   # }
#   # g + xlim(-3, 3) + scale_x_continuous(breaks=seq(-3, 3, 0.2))
#
#   # g = ggplot() + xlim(-3, 3)
#   # for (i in 1:2){
#   #   g <- local({
#   #     j <- i
#   #     shift <- (j - 1)
#   #     g + stat_function(fun = function(x) ((x - shift) < 1) * pmin(1, 1 - 0.5 - 0.5 * (x - shift)), col = j) +
#   #       stat_function(fun = function(x) ((x + shift) > -1) * pmin(1, 0.5 + 0.5 * (x + shift)), col = j+5)
#   #   })
#   # }
#   # g
#
#   # create the precision matrices
#   prec_mats <- vector("list", n)
#   for (i in 1:n){
#     prec_mats[[i]] <- common_str
#     for (j in 1:4){
#       prec_mats[[i]][j, j + 1] <- omega[i, j]
#     }
#   }
#
#   # symmetrize the precision matrices
#   true_precision <- lapply(prec_mats, function(mat) t(mat) + mat)
#
#   # library(covdepGE)
#   # for (i in 1:n){
#   #   print(i)
#   #   print(matViz(true_precision[[i]], incl_val = T) + ggtitle(i))
#   # }
#
#   # invert the precision matrices to get the covariance matrices
#   cov_mats <- lapply(true_precision, solve)
#
#   # generate the data using the covariance matrices
#   data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p)))
#
#   return(list(X = data_mat, Z = Z, true_precision = true_precision))
# }
