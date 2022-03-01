setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_demos_experiments/hyperparameter_specification/condition_number_analysis")
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

set.seed(1)
n <- 180; p <- 24
n_trials <- 100

results <- vector("list", n_trials)

names(results) <- paste0("trial", 1:n_trials)

doParallel::registerDoParallel(14)

pb <- txtProgressBar(min = 0, max = n_trials, style = 3)

for (j in 1:n_trials){

  # generate the data
  cont <- generate_continuous(p = p)
  X <- cont$data
  Z <- cont$covts

  # run the algorithm
  out <- tryCatch(covdepGE(X, Z, max_iter = 1e2, parallel = T, warnings = F,
                           stop_cluster = F),
                  error = function(msg) as.character(msg))

  # save the data and the results
  results[[j]] <- list(data = cont, results = out)

  setTxtProgressBar(pb, j)
}

close(pb)

doParallel::stopImplicitCluster()

Sys.time() - start

#save(results, file = "cond_number_models100.Rda")

library(ggplot2)
library(latex2exp)
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/analyses_demos_experiments/condition_number_analysis/cond_number_models.Rda")

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

# find ssq that are blown: either NA or in excess of 15
blown_sig <- lapply(ssq_nest, lapply, function(sigma) (sigma > 15 | is.na(sigma)))

max(unlist(ssq_nest)[!unlist(blown_sig)])
min(unlist(ssq_nest)[unlist(blown_sig)], na.rm = T)
length(unlist(ssq_nest))
sum(unlist(blown_sig))
impute_missing <- rep(max(unlist(ssq_nest), na.rm = T), sum(is.na(unlist(ssq_nest))))
sum(is.na(unlist(ssq_nest)))

# find the values of the blown sigmasq
ggplot(data.frame(ssq = sort(c(na.omit(unlist(ssq_nest)), impute_missing))), aes(ssq)) +
  geom_histogram(color = "black", fill = "tomato3") +
  scale_x_continuous(trans = "log10", n.breaks = 15) + scale_y_continuous(trans = "log10", n.breaks = 8) +
  theme_bw() + ggtitle(TeX("Distribution of the MAPE fitted $\\sigma^2$")) +
  xlab(TeX("$\\sigma^2$")) + theme(plot.title = element_text(hjust = 0.5))

# find counts of individuals by variable and trial
blown_ct_var <- lapply(blown_sig, sapply, sum)

# find counts by trials of variables with blown sigmas
blown_ct_tr <- sapply(blown_ct_var, function(ssq) sum(ssq > 0))

# logical for trials with at least one blown variable
blown_trials <- blown_ct_tr > 0
sum(blown_trials)

# visualize the counts of the number of variables that blew up for each trial
blown_tr_df <- data.frame(sort(blown_ct_tr[blown_trials]))
blown_tr_df$trial <- factor(sapply(strsplit(row.names(blown_tr_df), "trial"), `[[`, 2))
blown_tr_df$trial <- factor(blown_tr_df$trial, levels = blown_tr_df$trial)
names(blown_tr_df)[1] <- "blowups"
ggplot(blown_tr_df, aes(trial, blowups, fill = trial)) + geom_bar(stat = "identity") +
  theme_bw() + ggtitle("Number of variables with blowups by trial") + theme(plot.title = element_text(hjust = 0.5)) +
  ggsci::scale_fill_igv() + guides(fill = "none")

# find blown individual indices for each variable
blown_inds <- lapply(blown_sig[blown_trials], sapply, which)

# filter out the zero length variables from blown_inds
blown_inds0 <- lapply(blown_inds, function(trial) trial[sapply(trial, function(variable) length(variable) > 0)])
tr_var_cts <- sapply(unlist(blown_inds0, F), length)
tr_var_names <- strsplit(names(tr_var_cts), split = "\\.")
tr_var_df <- cbind.data.frame(t(sapply(tr_var_names, c)), tr_var_cts)
row.names(tr_var_df) <- NULL
colnames(tr_var_df) <- c("Trial", "Variable", "Number of individuals that blew up")
tr_var_df$var_num <- factor(sapply(strsplit(tr_var_df$Variable, "variable"), `[[`, 2))
tr_var_df$tr_num <- factor(sapply(strsplit(tr_var_df$Trial, "trial"), `[[`, 2), levels = levels(blown_tr_df$trial))
tr_var_df$label <- paste0("Trial ", tr_var_df$tr_num, ", Variable ", tr_var_df$var_num)
tr_var_df$label <- factor(tr_var_df$label, levels = tr_var_df$label)
ggplot(tr_var_df, aes(label, `Number of individuals that blew up`, fill = tr_num)) +
  geom_bar(stat = "identity") + ggsci::scale_fill_igv() +
  guides(fill = "none") + xlab("Trial and variable number") +
  ggtitle("Number of Individuals that blew up by trial and variable number") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90))


# all of the data
datas <- lapply(results, `[[`, "data")

# all of the data mats
data_mats <- lapply(datas, `[[`, "data")

# all data matrices with/ without at least one blown ind
blown_mats <- data_mats[blown_trials]
unblown_mats <- data_mats[!blown_trials]

# find the condition number for both
blown_cond_nums <- sapply(blown_mats, kappa)
unblown_cond_nums <- sapply(unblown_mats, kappa)
ggplot(data.frame(condition_number = blown_cond_nums),aes(condition_number)) +
  geom_histogram(color = "black", fill = "tomato3", binwidth = 0.1) + theme_bw() +
  ggtitle("Distribution of condition numbers for data with at least one blowup (unweighted)") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(2.9, 4.3))
summary(blown_cond_nums)
ggplot(data.frame(condition_number = unblown_cond_nums),aes(condition_number)) +
  geom_histogram(color = "black", fill = "forestgreen", binwidth = 0.1) + theme_bw() +
  ggtitle("Distribution of condition numbers for data with no blowups (unweighted)") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(2.9, 4.3))
summary(unblown_cond_nums)

# find the weights
wts <- lapply(summs, `[[`, "weights")

# find the square root of the weights
wts_2 <- lapply(wts, sqrt)

# find the weighted data_mats for everyone
wted_mats <- lapply(1:n_trials, function(tr_ind)
  lapply(1:p, function(var_ind)
    lapply(1:n, function(indv_ind) wts_2[[tr_ind]][ , indv_ind] * data_mats[[tr_ind]][ , -var_ind])))

# length(wted_mats)
# length(wted_mats[[1]])
# length(wted_mats[[1]][[1]])
#
# length(blown_sig)
# length(blown_sig[[1]])
# length(blown_sig[[1]][[1]])

# find the weighted matrices corresponding to the individuals that did and did not blow up
blown_wt_mat <- lapply(1:n_trials, function(trial_ind)
  lapply(1:p, function(var_ind) wted_mats[[trial_ind]][[var_ind]][blown_sig[[trial_ind]][[var_ind]]]))
unblown_wt_mat <- lapply(1:n_trials, function(trial_ind)
  lapply(1:p, function(var_ind) wted_mats[[trial_ind]][[var_ind]][!blown_sig[[trial_ind]][[var_ind]]]))
for (trial_ind in 1:n_trials){
  for (var_ind in 1:p){
    if (length(blown_wt_mat[[trial_ind]][[var_ind]]) > 0){
      names(blown_wt_mat[[trial_ind]][[var_ind]]) <- paste0("individual", which(blown_sig[[trial_ind]][[var_ind]]))
    }
    names(unblown_wt_mat[[trial_ind]][[var_ind]]) <- paste0("individual", which(!blown_sig[[trial_ind]][[var_ind]]))
  }
}
# length(unlist(unlist(blown_wt_mat, F), F))
# sum(unlist(blown_sig))
# length(unlist(unlist(unblown_wt_mat, F), F))
# sum(!unlist(blown_sig))

# unlist the blown and unblown mats so that the list consists just of matrices
blown_wt_mat2 <- unlist(unlist(blown_wt_mat, F), F)
unblown_wt_mat2 <- unlist(unlist(unblown_wt_mat, F), F)

# find the condition number for each of the blown up matrices
blown_cond <- sapply(blown_wt_mat2, kappa)
unblown_cond <- sapply(unblown_wt_mat2, kappa)
ggplot(data.frame(condition_number = blown_cond), aes(condition_number)) +
  geom_histogram(color = "black", fill = "tomato3", binwidth = 0.25) + theme_bw() +
  ggtitle("Distribution of condition numbers for weighted matrices that blew up") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(2.8, 9))
summary(blown_cond)
ggplot(data.frame(condition_number = unblown_cond), aes(condition_number)) +
  geom_histogram(color = "black", fill = "forestgreen", binwidth = 0.25) + theme_bw() +
  ggtitle("Distribution of condition numbers for weighted matrices that did not blow up") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(2.8, 9))
summary(unblown_cond)

# collect weighted matrices by individual
wt_mats2500 <- unlist(wted_mats, F)
wt_mat_indv <- lapply(1:n, function(indv_ind) lapply(1:(n_trials * p), function (ind25) wt_mats2500[[ind25]][[indv_ind]]))

# calculate the condition number for each of these matrices
cond_num_indv_180 <- lapply(wt_mat_indv, sapply, kappa)

# calculate the median condition number for each of the individuals
med_cond_num_180 <- sapply(cond_num_indv_180, median)

# visualize these condition numbers
ggplot(data.frame(Individual = 1:n, median_condition_number = med_cond_num_180), aes(
  Individual, median_condition_number)) + geom_bar(stat = "identity", fill = "tomato3") + theme_bw() +
  ggtitle("Median Condition Number By Individual") + theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(4.5, 6))

# find all of the unique individuals that had at least one blown up sigmasq
blown_indvs <- sort(unique(unlist(blown_inds)))
# table(unlist(blown_inds))

# vizualize the counts of these individuals
blown_ct_df <- data.frame(table(unlist(blown_inds)))
blown_ct_df$Var1 <- as.numeric(as.character(blown_ct_df$Var1))
new_rows <- cbind(setdiff(1:n, blown_ct_df$Var1), 0)
colnames(new_rows) <- colnames(blown_ct_df)
blown_ct_df <- rbind.data.frame(blown_ct_df, new_rows)
names(blown_ct_df) <- c("Individual", "Number of Blow Ups")
ggplot(blown_ct_df, aes(Individual, `Number of Blow Ups`)) +
  geom_bar(stat = "identity", fill = "tomato3") +
  ggtitle("Number of Blowups by Individual") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# var 14, trial 86
# tr86_wt_mat <- lapply(1:n, function(var_ind) data_mats[[86]][ , -14] * wts_2[[86]][ , var_ind])
# tr86_conds <- sapply(tr86_wt_mat, kappa)
# ggplot(data.frame(Individual = 1:n, condition_number = tr86_conds), aes(Individual, condition_number)) +
#   geom_bar(stat = "identity", fill = "tomato3") + coord_cartesian(ylim = c(4.2, 6.1)) +
#   theme_bw() + ggtitle("Condition Number by Individual, Trial 86, Variable 14") +
#   theme(plot.title = element_text(hjust = 0.5))

# for each of the individuals who has blown up, find the matrices corresponding
# to the variables where they blew up and where they did not blow up
blown_indvs_blow_mats <- lapply(blown_indvs, function(indv_idx)
  lapply(1:n_trials, function(trial_ind)
    lapply(1:p, function(var_ind) blown_wt_mat[[trial_ind]][[var_ind]][paste0("individual", indv_idx)])))
names(blown_indvs_blow_mats) <- paste0("individual", blown_indvs)
blown_indvs_unblow_mats <- lapply(blown_indvs, function(indv_idx)
  lapply(1:n_trials, function(trial_ind)
    lapply(1:p, function(var_ind) unblown_wt_mat[[trial_ind]][[var_ind]][paste0("individual", indv_idx)])))
names(blown_indvs_unblow_mats) <- paste0("individual", blown_indvs)

indv_171_blow_mats <- unlist(unlist(blown_indvs_blow_mats[["individual171"]], F), F)
length(indv_171_blow_mats[!sapply(indv_171_blow_mats, is.null)])
indv_1_blow_mats <- unlist(unlist(blown_indvs_blow_mats[[1]], F), F)
# length(indv_1_blow_mats[!sapply(indv_1_blow_mats, is.null)])
# sum(unlist(blown_inds) == 1)

# find condition numbers for each of the matrices
blown_indvs_blow_conds <- lapply(blown_indvs_blow_mats, lapply, lapply, lapply,
                                 function(mat) if (!is.null(mat)) kappa(mat))
blown_indvs_unblow_conds <- lapply(blown_indvs_unblow_mats, lapply, lapply, lapply,
                                   function(mat) if (!is.null(mat)) kappa(mat))


# visualize these
ggplot(data.frame(condition_numbers = unlist(blown_indvs_blow_conds)), aes(condition_numbers)) +
  geom_histogram(binwidth = 0.25, color = "black", fill = "tomato3") +
  ggtitle("Condition numbers for matrices that blew up for individuals that blew up") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(3, 9))
summary(unlist(blown_indvs_blow_conds))
ggplot(data.frame(condition_numbers = unlist(blown_indvs_unblow_conds)), aes(condition_numbers)) +
  geom_histogram(binwidth = 0.25, color = "black", fill = "forestgreen") +
  ggtitle("Condition numbers for matrices that did not blow up for individuals that blew up") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(3, 9))
summary(unlist(blown_indvs_unblow_conds))

# summary(unlist(blown_indvs_blow_conds[["individual171"]]))
# summary(unlist(blown_indvs_unblow_conds[["individual171"]]))

# find condition numbers for each of the individuals with nulls removed
blown_indvs_blow_conds2 <- lapply(blown_indvs_blow_conds, unlist)
blown_indvs_unblow_conds2 <- lapply(blown_indvs_unblow_conds, unlist)

# find medians for each
# med_blow_conds <- sapply(blown_indvs_blow_conds2, median)
# med_unblow_conds <- sapply(blown_indvs_unblow_conds2, median)
# round(med_blow_conds - med_unblow_conds, 2)

# find means for each
mean_blow_conds <- sapply(blown_indvs_blow_conds2, mean)
mean_unblow_conds <- sapply(blown_indvs_unblow_conds2, mean)
ggplot(data.frame(X = mean_blow_conds - mean_unblow_conds), aes(X)) +
  geom_histogram(binwidth = .05, color = "black", fill = "forestgreen") + theme_bw() +
  ggtitle("Condition number for blown up matrices minus non-blown up matrices by Individual") +
  theme(plot.title = element_text(hjust = 0.5))
summary(mean_blow_conds - mean_unblow_conds)

# visualize the condition numbers for individual 171
ggplot(data.frame(blown_indvs_blow_conds2["individual171"]), aes(individual171)) +
  geom_histogram(binwidth = 0.25, color = "black", fill = "tomato3") + theme_bw() +
  ggtitle("Condition number for matrices that blew up, individual 171") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(xlim = c(4, 9))
ggplot(data.frame(blown_indvs_unblow_conds2["individual171"]), aes(individual171)) +
  geom_histogram(binwidth = 0.25, color = "black", fill = "forestgreen") + theme_bw() +
  ggtitle("Condition number for matrices that did not blow up, individual 171") +
  theme(plot.title = element_text(hjust = 0.5))
