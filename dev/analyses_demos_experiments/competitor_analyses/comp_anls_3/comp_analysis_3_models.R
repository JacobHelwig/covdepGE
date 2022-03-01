setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_demos_experiments/competitor_analyses/comp_anls_3")
start <- Sys.time()

# function for generating the data and the covariates
generate_continuous <- function(n1 = 60, n2 = 60, n3 = 60, p = 4){

  # create covariate for individuals in each of the three intervals

  # define the dimensions of the data
  n <- sum(n1, n2, n3)

  # define the limits of the intervals
  limits1 <- c(-3, -1)
  limits2 <- c(-1, 1)
  limits3 <- c(1, 3)

  # define the covariate values within each interval
  z1 <- runif(n1, limits1[1], limits1[2])
  z2 <- runif(n2, limits2[1], limits2[2])
  z3 <- runif(n3, limits3[1], limits3[2])
  Z <- matrix(sort(c(z1, z2, z3)), n, 1)

  # create precision matrices

  # the shared part of the structure for all three intervals is a 2 on the
  # diagonal and a 1 in the (2, 3) position
  common_str <- diag(p + 1)
  common_str[2, 3] <- 1

  # define constants for the structure of interval 2
  beta1 <- diff(limits2)^-1
  beta0 <- -limits2[1] * beta1

  # interval 2 has two different linear functions of Z in the (1, 2) position
  # and (1, 3) positions; define structures for each of these components
  int2_str12 <- int2_str13 <- matrix(0, p + 1, p + 1)
  int2_str12[1, 2] <- int2_str13[1, 3] <- 1

  # define the precision matrices for each of the individuals in interval 2
  int2_prec <- lapply(z2, function(z) common_str +
                        ((1 - beta0 - beta1*z)*int2_str12) +
                        ((beta0 + beta1*z)*int2_str13))

  # interval 1 has a 1 in the (1, 2) and interval 3 has a 1 in the (1, 3) position;
  # define structures for each of these components
  int1_str12 <- int3_str13 <- matrix(0, p + 1, p + 1)
  int1_str12[1, 2] <- int3_str13[1, 3] <- 1

  # define the precision matrices for each of the individuals in interval 1 and interval 3
  int1_prec <- rep(list(common_str + int1_str12), n1)
  int3_prec <- rep(list(common_str + int3_str13), n3)

  # put all of the precision matrices into one list
  prec_mats <- c(int1_prec, int2_prec, int3_prec)

  # symmetrize the precision matrices
  prec_mats <- lapply(prec_mats, function(mat) t(mat) + mat)

  # invert the precision matrices to get the covariance matrices
  cov_mats <- lapply(prec_mats, solve)

  # generate the data using the covariance matrices
  data_mat <- t(sapply(cov_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))

  return(list(data = data_mat, covts = Z, true_precision = prec_mats))
}

library(covdepGE)
library(mclust)

# define experiment parameters
set.seed(1)
n <- 180; p <- 4
n_trials <- 100
sparsity_levels <- c(5, 10, 15)

# results storage for each sparsity level
results <- vector("list", 3)
names(results) <- paste0("sparsity_level", sparsity_levels)

doParallel::registerDoParallel(8)

pb <- txtProgressBar(min = 0, max = n_trials * 3, style = 3)

for (sparsity_level in sparsity_levels){

  # results for sparsity_level
  results_spl <- vector("list", n_trials)
  names(results_spl) <- paste0("trial", 1:n_trials)

  for (j in 1:n_trials){

    # storage for the j-th trial
    results_j <- vector("list", 3)
    names(results_j) <- c("data", "dependent", "independent")

    # generate the data
    cont <- generate_continuous(n2 = sparsity_level)
    X <- cont$data
    Z <- cont$covts

    # save the data
    results_j$data <- cont

    # run the covariate dependent method
    out <- tryCatch(covdepGE(X, Z, max_iter = 1e2, parallel = T, warnings = F,
                             stop_cluster = F),
                    error = function(msg) as.character(msg))

    # save the results
    results_j$dependent <- out

    # run the covariate independent method
    # apply Gaussian Mixture model clustering; selects number of clusters based on
    # the model that results in the best BIC
    gmm <- Mclust(Z, verbose = F)

    # find number of clusters in final clustering
    num_clusters <- length(unique(gmm$classification))

    # storage for the graph for each of the clusters and clustering
    out <- vector("list", num_clusters + 1)
    names(out) <- c("clusters", paste0("cluster", 1:num_clusters))
    out$clusters <- gmm$classification

    # iterate over each of the clusters identified by GMM
    for (k in 1:num_clusters) {

      # fix the datapoints in the k-th cluster
      data_mat_k <- X[gmm$classification == k, ]

      # apply the GGM using covdepGE with constant Z, save the resulting graph
      out[[paste0("cluster", k)]] <- tryCatch(
        covdepGE(data_mat_k, rep(0, nrow(data_mat_k)), max_iter = 1e2,
                 scale = F, kde = F, parallel = T, warnings = F, stop_cluster = F),
        error = function(msg) as.character(msg))
    }

    # save the independent results
    results_j$independent <- out

    # save the results for trial_j
    results_spl[[paste0("trial", j)]] <- results_j

    # update the progress bar
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
  }

  # save the results for the sparsity level
  results[[paste0("sparsity_level", sparsity_level)]] <- results_spl

}

close(pb)

doParallel::stopImplicitCluster()

Sys.time() - start

#save(results, file = "competitor_models3.Rda")

library(ggplot2)
library(latex2exp)
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_demos_experiments/competitor_analyses/comp_anls_3/competitor_models3.Rda")

# storage for performance by sparsity level
sparsity_levels <- names(results)
trials <- names(results$sparsity_level5)
performance <- vector("list", length(sparsity_levels))
names(performance) <- sparsity_levels

for (sparsity_level in sparsity_levels){

  # storage for performance by trial
  perf_trial <- vector("list", length(trials))
  names(perf_trial) <- trials

  for (trial in trials){

    # get dimensions of the data
    data_dim <- dim(results[[sparsity_level]][[trial]]$data$data)
    n <- data_dim[1]; p <- data_dim[2]

    # get the true precision structure
    tru_prec <- results[[sparsity_level]][[trial]]$data$true_precision
    tru_prec <- lapply(tru_prec, function(mat) (mat - diag(diag(mat)) != 0) * 1)

    # get the dependent results
    dep_graphs <- results[[sparsity_level]][[trial]]$dependent$graphs

    # get independent results
    indep_res <- results[[sparsity_level]][[trial]]$independent
    num_clust <- length(unique(indep_res$clusters))
    indep_graphs <- vector("list", n)
    for (clust in 1:num_clust){

      # fix the result for cluster clust
      clust_res <- indep_res[[paste0("cluster", clust)]]

      # find the graph for this cluster and assign it to the individuals in this cluster
      if (length(clust_res$unique_graphs) != 1) stop("Unique graph error")
      clust_gr <- clust_res$unique_graphs$graph1$graph
      clust_indiv <- which(indep_res$clusters == clust)
      indep_graphs[clust_indiv] <- list(clust_gr)
    }

    # calculate the number of true positives/ negatives
    tru1 <- sum(unlist(tru_prec))
    tru0 <- sum(unlist(tru_prec) == 0)

    # true pos/ neg for dep/ indep
    tru1_dep <- sum(sapply(1:n, function(j) sum((tru_prec[[j]] == 1) & (dep_graphs[[j]] == 1))))
    tru0_dep <- sum(sapply(1:n, function(j) sum((tru_prec[[j]] == 0) & (dep_graphs[[j]] == 0))))
    tru1_indep <- sum(sapply(1:n, function(j) sum((tru_prec[[j]] == 1) & (indep_graphs[[j]] == 1))))
    tru0_indep <- sum(sapply(1:n, function(j) sum((tru_prec[[j]] == 0) & (indep_graphs[[j]] == 0))))

    # sensitivity/ specificity/ accuracy for dep/ indep
    perf <- vector("list", 6)
    metrics <- c("sensitivity_", "specificity_", "accuracy_")
    names(perf) <- c(paste0(metrics, c("dependent")), paste0(metrics, "independent"))
    perf[["sensitivity_dependent"]] <- tru1_dep / tru1
    perf[["specificity_dependent"]] <- tru0_dep / tru0
    perf[["accuracy_dependent"]] <- (tru1_dep + tru0_dep) / (tru1 + tru0)
    perf[["sensitivity_independent"]] <- tru1_indep / tru1
    perf[["specificity_independent"]] <- tru0_indep / tru0
    perf[["accuracy_independent"]] <- (tru1_indep + tru0_indep) / (tru1 + tru0)

    # save the trial performance
    perf_trial[[trial]] <- perf
  }

  # save the sparsity level performance
  performance[sparsity_level] <- list(perf_trial)
}

summary(sapply(performance$sparsity_level5, `[[`, "sensitivity_dependent") - sapply(performance$sparsity_level5, `[[`, "sensitivity_independent"))
summary(sapply(performance$sparsity_level15, `[[`, "sensitivity_dependent") - sapply(performance$sparsity_level15, `[[`, "sensitivity_independent"))
summary(sapply(performance$sparsity_level10, `[[`, "sensitivity_dependent") - sapply(performance$sparsity_level10, `[[`, "sensitivity_independent"))
