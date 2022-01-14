rm(list = ls())
start <- Sys.time()
library(covdepGE)
source("generate_data.R")

num_rep <- 500

# store the resulting metrics for each model
elbo1 <- elbo2 <- sens1 <- sens2 <- spec1 <- spec2 <- accu1 <- accu2 <- rep(NA, num_rep)

# store final prior prob of inclusion for varbvs tuple model
pi1 <- matrix(NA, 5, num_rep)

# store sigma and sigmabetasq candidates and ranges from the varbvs tuple model
sigmasq_cands <- sigmabetasq_cands <- matrix(NA, 100, num_rep)
sigmasq_ranges <- sigmabetasq_ranges <- matrix(NA, 5, num_rep)

doParallel::registerDoParallel(7)

pb <- utils::txtProgressBar(0, num_rep, style = 3)
utils::setTxtProgressBar(pb, 0)
for (j in 1:num_rep){

  # generate the data with seed = j
  cont <- generate_continuous(seed = j)
  X <- cont$data
  Z <- cont$covts

  # run CAVI with constant values of sigmasq and pi for each variable
  out1 <- covdepGE(X, Z, var_min = 1e-3, var_max = 1, n_sigma = 20,
                   max_iter_grid = 100, max_iter_final = 1000, CS = T,
                   parallel = T, stop_cluster = F, warnings = F)

  # run CAVI with hyperparameters tuples estimated by varbvs
  out2 <- covdepGE(X, Z, max_iter_grid = 100, max_iter_final = 1000,
                   parallel = T, stop_cluster = F, warnings = F)

  # save the resulting pi values from out2
  pi1[ , j] <- sapply(out2$CAVI_details, `[[`, "pi")

  # save the sigma and sigmabetasq candidates from out2
  sigma_cands <- lapply(out2$CAVI_details, function(variable)
    variable$hyperparameters[ , c("sigmasq", "sigmabetasq")])
  sigma_cands_mat <- do.call(rbind, sigma_cands)
  sigmasq_cands[ , j] <- sigma_cands_mat[ , "sigmasq"]
  sigmabetasq_cands[ , j] <- sigma_cands_mat[ , "sigmabetasq"]

  # find the ranges of the sigma and sigmabetasq candidates from out2
  sigmasq_ranges[ , j] <- sapply(sigma_cands, function(variable) range(
    variable[ , "sigmasq"]) %*% c(-1, 1))
  sigmabetasq_ranges[ , j] <- sapply(sigma_cands, function(variable) range(
    variable[ , "sigmabetasq"]) %*% c(-1, 1))

  # save elbo for each
  elbo1[j] <- out1$model_details$ELBO
  elbo2[j] <- out2$model_details$ELBO

  # get the true graphs
  true_graphs <- lapply(cont$true_precision, function(prec_mat) ((prec_mat - diag(diag(prec_mat))) != 0) * 1)

  # number of 1's and 0's
  num1 <- sum(unlist(true_graphs))
  num0 <- sum(-unlist(true_graphs) + 1)
  num1 + num0 == length(unlist(true_graphs))

  # get sensitivity (true 1)
  true1_1 <- sum(sapply(1:nrow(X), function(gr_idx) sum(
    true_graphs[[gr_idx]] == out1$graphs[[gr_idx]] & true_graphs[[gr_idx]] == 1)))
  true1_2 <- sum(sapply(1:nrow(X), function(gr_idx) sum(
    true_graphs[[gr_idx]] == out2$graphs[[gr_idx]] & true_graphs[[gr_idx]] == 1)))
  sens1[j] <- true1_1 / num1
  sens2[j] <- true1_2 / num1

  # get specificity (true 0)
  true0_1 <- sum(sapply(1:nrow(X), function(gr_idx) sum(
    true_graphs[[gr_idx]] == out1$graphs[[gr_idx]] & true_graphs[[gr_idx]] == 0)))
  true0_2 <- sum(sapply(1:nrow(X), function(gr_idx) sum(
    true_graphs[[gr_idx]] == out2$graphs[[gr_idx]] & true_graphs[[gr_idx]] == 0)))
  spec1[j] <- true0_1 / num0
  spec2[j] <- true0_2 / num0

  # get accuracy
  accu1[j] <- (true0_1 + true1_1) / (num0 + num1)
  accu2[j] <- (true0_2 + true1_2) / (num0 + num1)

  utils::setTxtProgressBar(pb, j)
}

close(pb)

res <- list(
  model1 = list(elbo = elbo1, sensitivity = sens1, specificity = spec1, accuracy = accu1),
  model2 = list(elbo = elbo2, sensitivity = sens2, specificity = spec2, accuracy = accu2),
  hyperparameters = list(pi_final = pi1, sigmasq_cands = sigmasq_cands,
                         sigmabetasq_cands = sigmabetasq_cands,
                         sigmasq_ranges = sigmasq_ranges,
                         sigmabetasq_ranges = sigmabetasq_ranges)
  )
#save(res, file = "~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_and_demos/hyp_spec_exper_results.Rda")
doParallel::stopImplicitCluster()
end <- Sys.time()
end - start

