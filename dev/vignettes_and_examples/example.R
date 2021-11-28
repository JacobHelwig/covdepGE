library(covdepGE)
set.seed(1)
n <- 100
p <- 4

# generate the extraneous covariate
Z_neg <- sort(runif(n / 2) * -1)
Z_pos <- sort(runif(n / 2))
Z <- c(Z_neg, Z_pos)
summary(Z)

# create true covariance structure for 2 groups: positive Z and negative Z
true_graph_pos <- true_graph_neg <- matrix(0, p + 1, p + 1)
true_graph_pos[1, 2] <- true_graph_pos[2, 1] <- 1
true_graph_neg[1, 3] <- true_graph_neg[3, 1] <- 1

# visualize the true covariance structures
gg_adjMat(true_graph_neg) + ggplot2::ggtitle("True graph for individuals with negative Z")
gg_adjMat(true_graph_pos, color1 = "steelblue") + ggplot2::ggtitle("True graph for individuals with positive Z")

# generate the covariance matrices as a function of Z
sigma_mats_neg <- lapply(Z_neg, function(z) z * true_graph_neg + diag(p + 1))
sigma_mats_pos <- lapply(Z_pos, function(z) z * true_graph_pos + diag(p + 1))
sigma_mats <- c(sigma_mats_neg, sigma_mats_pos)

# generate the data using the covariance matrices
data_mat <- t(sapply(sigma_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))

# visualize the sample correlation
gg_adjMat(abs(cor(data_mat[1:(n / 2), ])) - diag(p + 1))
gg_adjMat(abs(cor(data_mat[(n / 2 + 1):n, ])) - diag(p + 1), color1 = "dodgerblue")

# fit the model
out <- covdepGE(data_mat, Z, max_iter = 1000)

# analyze results
gg_adjMat(out, 1)
gg_adjMat(out, 50, color1 = "tomato")
gg_adjMat(out, 54, color1 = "steelblue")
gg_adjMat(out, 100, color1 = "dodgerblue")

gg_inclusionCurve(out, 1, 2)
gg_inclusionCurve(out, 1, 3, point_color = "dodgerblue")

# find sensitivity, specificity, and accuracy

# true positives
TP_neg <- sum(sapply(out$graphs[1:(n / 2)], function(graph) sum(graph == 1 & true_graph_neg == 1)))
TP_pos <- sum(sapply(out$graphs[(n / 2 + 1):n], function(graph) sum(graph == 1 & true_graph_pos == 1)))
TP <- TP_neg + TP_pos

# total positives
num_pos <- sum(true_graph_pos) * n / 2 + sum(true_graph_neg) * n / 2

# true negatives
TN_neg <- sum(sapply(out$graphs[1:(n / 2)], function(graph) sum(graph == 0 & true_graph_neg == 0)))
TN_pos <- sum(sapply(out$graphs[(n / 2 + 1):n], function(graph) sum(graph == 0 & true_graph_pos == 0)))
TN <- TN_neg + TN_pos

# total negatives
num_neg <- length(true_graph_pos) * n - num_pos

(sensitivity <- TP / num_pos)
(specificity <- TN / num_neg)
(accuracy <- (TN + TP) / (num_pos + num_neg))
