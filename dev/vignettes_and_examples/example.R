set.seed(1)
n <- 100
p <- 4

# generate the extraneous covariate
Z_neg <- runif(n / 2, -2, 0)
Z_pos <- runif(n / 2, 0, 2)
Z <- sort(c(Z_neg, Z_pos))

# true covariance structure for individuals with a positive extraneous covariate
true_graph_pos <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, rep(0, 13)),
                         p + 1, p + 1)

# true covariance structure for individuals with a negative extraneous covariate
true_graph_neg <- matrix(c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, rep(0, 13)),
                         p + 1, p + 1)

# generate the covariance matrix for both groups
noise <- runif(length(sigma), 0.25, 0.5) * sample(c(-1, 1), length(sigma), T)
sigma_neg <- true_graph_pos * noise
sigma_pos <- true_graph_neg * noise

# generate the covariance matrix for each individual by multiplying their
# extraneous covariate
sigma_mats <- lapply(Z, `*`, sigma)
sigma_mats <- lapply(sigma_mats, `+`, diag(p + 1))
data_mat <- t(sapply(sigma_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))
out <- covdepGE(data_mat, Z)
gg_adjMat(out, 8)
gg_inclusionCurve(out, 1, 2) + ggplot2::ylim(c(0, 1))

total_true_positives <- sum(sapply(out$graphs, function(graph) sum(graph == 1 & true_graph == 1)))
total_true_negatives <- sum(sapply(out$graphs, function(graph) sum(graph == 0 & true_graph == 0)))
total_positives <-

