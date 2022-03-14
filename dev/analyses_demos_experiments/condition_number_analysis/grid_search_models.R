setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_demos_experiments/condition_number_analysis")
start <- Sys.time()

# find the data for the trials resulting in instability
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_demos_experiments/condition_number_analysis/cond_number_models.Rda")
results0 <- results
# get dimensions of the data
dt_dim <- dim(results$trial1$data$data)
n <- dt_dim[1]; p <- dt_dim[2]
n_trials <- length(results)

# get all of the summaries
summs <- lapply(results, `[[`, "results")

# find the sigmasq for each model
ssq <- lapply(lapply(summs, `[[`, "hyperparameters"), `[[`, "sigmasq")

# unfold the sigmasq matrices into a n_trials x p nested list of n-vectors
ssq_nest <- lapply(1:n_trials, function(trial_ind) lapply(1:p, function(var_ind) ssq[[trial_ind]][ , var_ind]))
names(ssq_nest) <- paste0("trial", 1:n_trials)
for (j in 1:n_trials) names(ssq_nest[[j]]) <- paste0("variable", 1:p)

# find ssq that are blown: either NA or in excess of  1
blown_sig <- lapply(ssq_nest, lapply, function(sigma) (sigma > 1 | is.na(sigma)))

# find counts of individuals by variable and trial
blown_ct_var <- lapply(blown_sig, sapply, sum)

# find counts by trials of variables with blown sigmas
blown_ct_tr <- sapply(blown_ct_var, function(ssq) sum(ssq > 0))

# logical for trials with at least one blown variable
blown_trials <- blown_ct_tr > 0

# get the data for the blown trials
blown_res <- results[blown_trials]
blown_data <- lapply(blown_res, `[[`, "data")

library(covdepGE)
n_trials <- length(blown_data)
results <- vector("list", n_trials)
names(results) <- paste0("trial", 1:n_trials)

doParallel::registerDoParallel(13)

pb <- txtProgressBar(min = 0, max = n_trials, style = 3)

for (j in 1:n_trials){

  # get the data
  cont <- blown_data[[j]]
  X <- cont$data
  Z <- cont$covts

  # run the algorithm
  out <- tryCatch(covdepGE(X, Z, max_iter = 1e2, parallel = T, warnings = F,
                           stop_cluster = F, alpha_tol = 1e-12, grid_search = F),
                  error = function(msg) as.character(msg))

  # save the data and the results
  results[[j]] <- list(data = cont, results = out)

  setTxtProgressBar(pb, j)
}

close(pb)

doParallel::stopImplicitCluster()

Sys.time() - start

save(results, file = paste0("grid_sch_models_no_ELBO", n_trials, ".Rda"))

load("grid_sch_models27.Rda")
results_gr_sch <- results
load("analyses_demos_experiments/condition_number_analysis/grid_sch_models_no_ELBO27.Rda")
results_no_elbo <- results

for (j in 1:length(results_gr_sch)){
  res_j_grs <- results_gr_sch[[j]]$results
  time_grs <- as.numeric(res_j$model_details$elapsed, units = "secs")

  res_j_nelbo <- results_no_elbo[[j]]$results
  time_nelbo <- as.numeric(res_j_nelbo$model_details$elapsed, units = "secs")

  identical(res_j_grs$graphs, res_j_nelbo$graphs)
}

dt <- blown_data[[1]]

# out_gs <- covdepGE(dt$data, dt$covts, parallel = T, num_workers = 13, alpha_tol = 1e-12, elbo_tol = 1e-3)
# out_ngs <- covdepGE(dt$data, dt$covts, parallel = T, num_workers = 13, alpha_tol = 1e-12, grid_search = F, elbo_tol = 1e-3)
#
#
# mean(sapply(lapply(lapply(results0, `[[`, "results"), `[[`, "model_details"), `[[`, "ELBO"), na.rm = T)
# min(sapply(lapply(lapply(results_gr_sch, `[[`, "results"), `[[`, "model_details"), `[[`, "ELBO"), na.rm = T)
# mean(sapply(lapply(lapply(results0, `[[`, "results"), `[[`, "model_details"), `[[`, "num_unique"), na.rm = T)
# mean(sapply(lapply(lapply(results_gr_sch, `[[`, "results"), `[[`, "model_details"), `[[`, "num_unique"), na.rm = T)
