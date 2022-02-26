load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_demos_experiments/condition_number_analysis/cond_number_models.Rda")
results1000 <- results
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_demos_experiments/condition_number_analysis/cond_number_models100.Rda")
results100 <- results

# get dimensions of the data
dt_dim <- dim(results$trial1$data$data)
n <- dt_dim[1]; p <- dt_dim[2]
n_trials <- length(results)

# get all of the summaries
summs100 <- lapply(results100, `[[`, "results")

# get all of the CAVI details
cavi_dets100 <- lapply(summs100, `[[`, "CAVI_details")

# get all of the pi_stable_iter
pi_stable100 <- sapply(unlist(cavi_dets100, F), `[[`, "pi_stable_iter")

summary(pi_stable100)
sum(pi_stable100 > 1)
which(pi_stable100 > 1)

# get all of the summaries
summs1000 <- lapply(results1000, `[[`, "results")

# get all of the CAVI details
cavi_dets1000 <- lapply(summs1000, `[[`, "CAVI_details")

# get all of the pi_stable_iter
pi_stable1000 <- sapply(unlist(cavi_dets1000, F), `[[`, "pi_stable_iter")

summary(pi_stable1000)
sum(pi_stable1000 > 1)
which(pi_stable1000 > 1)
