# ------------------------------------------------------------------------------
# discrete covariate independent analysis
rm(list = ls())
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study_graph_learning")
library(kableExtra)

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each
exper <- "disc_cov_indep_"
trial_str <- "ntrials50_"
p <- 11
lambdas <- c(3, 15)
n1s <- c(50, 80)
path <- "./experiments/discrete/covariate independent"
files <- list.files(path)
prec <- 4
df <- list()
res <- list()
methods <- list(covdepGE = "W-PL", mgm = "mgm", varbvs = "CS")
for (lambda in lambdas){
  for (n1 in n1s){

    # if (lambda == 3 & n1 == 80) next

    # format file name
    n2 <- 100 - n1
    exp_name <- paste0(exper, trial_str, "p", p, "_n1_", n1, "_n2_", n2, "_lambda", lambda)

    # load results
    file_name <- files[startsWith(files, exp_name)]
    file_path <- file.path(path, file_name)
    load(file_path)
    results <- results[setdiff(names(results), "sample_data")]

    # process sensitivity results
    sens <- sapply(results, sapply, `[[`, "sens")
    res[[exp_name]]$sens <- sens
    sens_mean <- rowMeans(sens)
    max_sens_ind <- which.max(sens_mean)
    sens_mean <- sprintf(paste0("%.", prec, "f"), sens_mean)
    sens_mean[max_sens_ind] <- paste0("\\mathbf{", sens_mean[max_sens_ind], "}")
    sens_sd  <- sprintf(paste0("%.", prec, "f"), apply(sens, 1, sd))
    sens_str <- paste0(sens_mean, " (", sens_sd, ")")

    # process specificity results
    spec <- sapply(results, sapply, `[[`, "spec")
    res[[exp_name]]$spec <- spec
    spec_mean <- rowMeans(spec)
    max_spec_ind <- which.max(spec_mean)
    spec_mean <- sprintf(paste0("%.", prec, "f"), spec_mean)
    spec_mean[max_spec_ind] <- paste0("\\mathbf{", spec_mean[max_spec_ind], "}")
    spec_sd  <- sprintf(paste0("%.", prec, "f"), apply(spec, 1, sd))
    spec_str <- paste0(spec_mean, " (", spec_sd, ")")

    # combine summary strings
    perf_str <- cbind(sens_str, spec_str)
    perf_str <- matrix(paste0("$", perf_str, "$"), dim(perf_str))
    row.names(perf_str) <- row.names(spec)

    # create storage
    c_str <- paste0(c(lambda, rep("!", length(df))), collapse = "")
    n1_str <- paste0(c(n1, rep("!", length(df))), collapse = "")
    n2_str <- paste0(c(n2, rep("!", length(df))), collapse = "")
    df_exp <- data.frame(c = c_str, n1 = n1_str, n2 = n2_str,
                         method = c("covdepGE" ,"mgm", "varbvs"), sens = NA,
                         spec = NA)
    df_exp[ , c("sens", "spec")] <- perf_str[df_exp$method, ]
    df_exp$method <- unlist(methods[df_exp$method])

    df[[length(df) + 1]] <- df_exp
    rm("results")
  }
}

df <- Reduce(rbind, df)
colnames(df) <- c("$c$", "$n_1$", "$n_2$", "Method", "Sensitivity$(\\uparrow)$", "Specificity$(\\uparrow)$")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")

# ------------------------------------------------------------------------------
# discrete covariate free analysis
rm(list = ls())
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study_graph_learning")
library(kableExtra)

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each
exper <- "disc_cov_free_"
trial_str <- "ntrials50_"
p <- 11
lambdas <- c(3, 15)
samp_str <- "_n1_50_n2_50"
path <- "./experiments/discrete/covariate free"
files <- list.files(path)
prec <- 4
df <- list()
res <- list()
methods <- list(covdepGE = "W-PL", varbvs = "CS")
for (lambda in lambdas){

  # format file name
  exp_name <- paste0(exper, trial_str, "p", p, samp_str, "_lambda", lambda)

  # load results
  file_name <- files[startsWith(files, exp_name)]
  file_path <- file.path(path, file_name)
  load(file_path)
  results <- results[setdiff(names(results), "sample_data")]

  # process sensitivity results
  sens <- sapply(results, sapply, `[[`, "sens")
  res[[exp_name]]$sens <- sens
  sens_mean <- rowMeans(sens)
  max_sens_ind <- which.max(sens_mean)
  sens_mean <- sprintf(paste0("%.", prec, "f"), sens_mean)
  sens_mean[max_sens_ind] <- paste0("\\mathbf{", sens_mean[max_sens_ind], "}")
  sens_sd  <- sprintf(paste0("%.", prec, "f"), apply(sens, 1, sd))
  sens_str <- paste0(sens_mean, " (", sens_sd, ")")

  # process specificity results
  spec <- sapply(results, sapply, `[[`, "spec")
  res[[exp_name]]$spec <- spec
  spec_mean <- rowMeans(spec)
  max_spec_ind <- which.max(spec_mean)
  spec_mean <- sprintf(paste0("%.", prec, "f"), spec_mean)
  spec_mean[max_spec_ind] <- paste0("\\mathbf{", spec_mean[max_spec_ind], "}")
  spec_sd  <- sprintf(paste0("%.", prec, "f"), apply(spec, 1, sd))
  spec_str <- paste0(spec_mean, " (", spec_sd, ")")

  # combine summary strings
  perf_str <- cbind(sens_str, spec_str)
  perf_str <- matrix(paste0("$", perf_str, "$"), dim(perf_str))
  row.names(perf_str) <- row.names(spec)

  # create storage
  c_str <- paste0(c(lambda, rep("!", length(df))), collapse = "")
  df_exp <- data.frame(c = c_str, method = c("covdepGE", "varbvs"), sens = NA,
                       spec = NA)
  df_exp[ , c("sens", "spec")] <- perf_str[df_exp$method, ]
  df_exp$method <- unlist(methods[df_exp$method])

  df[[length(df) + 1]] <- df_exp
  rm("results")
}

df <- Reduce(rbind, df)
colnames(df) <- c("$c$", "Method", "Sensitivity$(\\uparrow)$", "Specificity$(\\uparrow)$")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")
# ------------------------------------------------------------------------------
# discrete covariate dependent analysis
rm(list = ls())
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study_graph_learning")
library(kableExtra)

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each
exper <- "disc_cov_dep_"
trial_str <- "ntrials50_"
dims <- c(11, 31, 51)
lambdas <- c(3, 15)
n1s <- c(50, 80)
path <- "./experiments/discrete/covariate dependent"
files <- list.files(path)
prec <- 4
df <- list()
methods <- list(covdepGE = "W-PL", mgm = "mgm", varbvs = "CS")
for (p in dims){
  for (lambda in lambdas){
    for (n1 in n1s){

      if (p != 11 & (lambda == 3 | n1 == 80)) next
      if (lambda == 3 & n1 == 80) next

      # format file name
      n2 <- 100 - n1
      exp_name <- paste0(exper, trial_str, "p", p, "_n1_", n1, "_n2_", n2, "_lambda", lambda)

      # load results
      file_name <- files[startsWith(files, exp_name)]
      file_path <- file.path(path, file_name)
      load(file_path)
      results <- results[setdiff(names(results), "sample_data")]

      # process sensitivity results
      sens <- sapply(results, sapply, `[[`, "sens")
      sens_mean <- rowMeans(sens)
      max_sens_ind <- which.max(sens_mean)
      sens_mean <- sprintf(paste0("%.", prec, "f"), sens_mean)
      sens_mean[max_sens_ind] <- paste0("\\mathbf{", sens_mean[max_sens_ind], "}")
      sens_sd  <- sprintf(paste0("%.", prec, "f"), apply(sens, 1, sd))
      sens_str <- paste0(sens_mean, " (", sens_sd, ")")

      # process specificity results
      spec <- sapply(results, sapply, `[[`, "spec")
      spec_mean <- rowMeans(spec)
      max_spec_ind <- which.max(spec_mean)
      spec_mean <- sprintf(paste0("%.", prec, "f"), spec_mean)
      spec_mean[max_spec_ind] <- paste0("\\mathbf{", spec_mean[max_spec_ind], "}")
      spec_sd  <- sprintf(paste0("%.", prec, "f"), apply(spec, 1, sd))
      spec_str <- paste0(spec_mean, " (", spec_sd, ")")

      # combine summary strings
      perf_str <- cbind(sens_str, spec_str)
      perf_str <- matrix(paste0("$", perf_str, "$"), dim(perf_str))
      row.names(perf_str) <- row.names(spec)

      # create storage
      p_str <- paste0(c(p - 1, rep("!", length(df))), collapse = "")
      c_str <- paste0(c(lambda, rep("!", length(df))), collapse = "")
      n1_str <- paste0(c(n1, rep("!", length(df))), collapse = "")
      n2_str <- paste0(c(n2, rep("!", length(df))), collapse = "")
      df_exp <- data.frame(p = p_str, c = c_str, n1 = n1_str, n2 = n2_str,
                           method = c("covdepGE" ,"mgm", "varbvs"), sens = NA,
                           spec = NA)
      df_exp[ , c("sens", "spec")] <- perf_str[df_exp$method, ]
      df_exp$method <- unlist(methods[df_exp$method])

      df[[length(df) + 1]] <- df_exp
      rm("results")
    }
  }
}

df <- Reduce(rbind, df)
colnames(df) <- c("$p$", "$c$", "$n_1$", "$n_2$", "Method", "Sensitivity$(\\uparrow)$", "Specificity$(\\uparrow)$")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")

# ------------------------------------------------------------------------------
# discrete covariate dependent high dimensional analysis
rm(list = ls())
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study_graph_learning")
library(kableExtra)

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each
exper <- "disc_cov_dep_"
trial_str <- "ntrials50_"
dims <- c(51, 101)
lambdas <- list()
lambdas[[51]] <- c(15, 3, 1, 0.6, 0.3)
lambdas[[101]] <- c(9, 6, 3, 1)
samp_str <- list()
samp_str[[51]] <- "_n1_40_n2_10_lambda"
samp_str[[101]] <- "_n1_30_n2_20_lambda"
path <- "./experiments/discrete/high dimensional"
files <- list.files(path)
prec <- 4
df <- list()
methods <- list(covdepGE = "W-PL", mgm = "mgm", varbvs = "CS")
for (p in dims){
  lambda <- lambdas[[p]]
  for (lamb in lambda){

    # format file name
    exp_name <- paste0(exper, trial_str, "p", p, samp_str[[p]], lamb, "_")

    # load results
    file_name <- files[startsWith(files, exp_name)]
    file_path <- file.path(path, file_name)
    load(file_path)
    results <- results[setdiff(names(results), "sample_data")]

    # process sensitivity results
    sens <- sapply(results, sapply, `[[`, "sens")
    sens_mean <- rowMeans(sens)
    max_sens_ind <- which.max(sens_mean)
    sens_mean <- sprintf(paste0("%.", prec, "f"), sens_mean)
    sens_mean[max_sens_ind] <- paste0("\\mathbf{", sens_mean[max_sens_ind], "}")
    sens_sd  <- sprintf(paste0("%.", prec, "f"), apply(sens, 1, sd))
    sens_str <- paste0(sens_mean, " (", sens_sd, ")")

    # process specificity results
    spec <- sapply(results, sapply, `[[`, "spec")
    spec_mean <- rowMeans(spec)
    max_spec_ind <- which.max(spec_mean)
    spec_mean <- sprintf(paste0("%.", prec, "f"), spec_mean)
    spec_mean[max_spec_ind] <- paste0("\\mathbf{", spec_mean[max_spec_ind], "}")
    spec_sd  <- sprintf(paste0("%.", prec, "f"), apply(spec, 1, sd))
    spec_str <- paste0(spec_mean, " (", spec_sd, ")")

    # combine summary strings
    perf_str <- cbind(sens_str, spec_str)
    perf_str <- matrix(paste0("$", perf_str, "$"), dim(perf_str))
    row.names(perf_str) <- row.names(spec)

    # create storage
    p_str <- p #paste0(c(p - 1, rep("!", length(df))), collapse = "")
    c_str <- lamb #paste0(c(lambda, rep("!", length(df))), collapse = "")
    df_exp <- data.frame(p = p_str, c = c_str,
                         method = c("covdepGE" ,"mgm", "varbvs"), sens = NA,
                         spec = NA)
    df_exp[ , c("sens", "spec")] <- perf_str[df_exp$method, ]
    df_exp$method <- unlist(methods[df_exp$method])

    df[[length(df) + 1]] <- df_exp
    rm("results")
  }
}

df <- Reduce(rbind, df)
colnames(df) <- c("$p$", "$c$", "$n_1$", "$n_2$", "Method", "Sensitivity$(\\uparrow)$", "Specificity$(\\uparrow)$")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")

# ------------------------------------------------------------------------------
# continuous covariate dependent analysis
rm(list = ls())
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study_graph_learning")
library(kableExtra)

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each
exper <- "cont_cov_dep_"
trial_str <- "ntrials70_"
samp_str <- "_n1_50_n2_50_n3_50"
dims <- c(11, 31, 51)
path <- "./experiments/continuous/covariate dependent"
files <- list.files(path)
prec <- 4
df <- list()
methods <- list(covdepGE = "W-PL", loggle = "loggle", mgm = "mgm")
for (p in dims){

  # format file name
  exp_name <- paste0(exper, trial_str, "p", p, samp_str)

  # load results
  file_name <- files[startsWith(files, exp_name)]
  file_path <- file.path(path, file_name)
  load(file_path)
  results <- results[1:50]
  results <- results[setdiff(names(results), "sample_data")]

  # process sensitivity results
  sens <- sapply(results, sapply, `[[`, "sens")
  sens_mean <- rowMeans(sens)
  max_sens_ind <- which.max(sens_mean)
  sens_mean <- sprintf(paste0("%.", prec, "f"), sens_mean)
  sens_mean[max_sens_ind] <- paste0("\\mathbf{", sens_mean[max_sens_ind], "}")
  sens_sd  <- sprintf(paste0("%.", prec, "f"), apply(sens, 1, sd))
  sens_str <- paste0(sens_mean, " (", sens_sd, ")")

  # process specificity results
  spec <- sapply(results, sapply, `[[`, "spec")
  spec_mean <- rowMeans(spec)
  max_spec_ind <- which.max(spec_mean)
  spec_mean <- sprintf(paste0("%.", prec, "f"), spec_mean)
  spec_mean[max_spec_ind] <- paste0("\\mathbf{", spec_mean[max_spec_ind], "}")
  spec_sd  <- sprintf(paste0("%.", prec, "f"), apply(spec, 1, sd))
  spec_str <- paste0(spec_mean, " (", spec_sd, ")")

  # combine summary strings
  perf_str <- cbind(sens_str, spec_str)
  perf_str <- matrix(paste0("$", perf_str, "$"), dim(perf_str))
  row.names(perf_str) <- row.names(spec)

  # create storage
  p_str <- paste0(c(p - 1, rep("!", length(df))), collapse = "")
  df_exp <- data.frame(p = p_str, method = c("covdepGE" , "loggle", "mgm"),
                       sens = NA, spec = NA)
  df_exp[ , c("sens", "spec")] <- perf_str[df_exp$method, ]
  df_exp$method <- unlist(methods[df_exp$method])

  df[[length(df) + 1]] <- df_exp
  rm("results")
}


df <- Reduce(rbind, df)
colnames(df) <- c("$p$", "Method", "Sensitivity$(\\uparrow)$", "Specificity$(\\uparrow)$")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")

# ------------------------------------------------------------------------------
# continuous multivariate covariate dependent analysis
rm(list = ls())
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study_graph_learning")
library(kableExtra)

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each
exper <- "cont_multi_cov_dep_"
trial_str <- "ntrials70_"
samp_str <- "_n_25"
dims <- c(11, 31)#, 51)
path <- "./experiments/continuous/multi covariate"
files <- list.files(path)
prec <- 4
df <- list()
methods <- list(covdepGE = "W-PL", covdepGE_sortZ = "tv W-PL", loggle = "loggle", mgm = "mgm")
for (p in dims){

  # format file name
  exp_name <- paste0(exper, trial_str, "p", p, samp_str)

  # load results
  file_name <- files[startsWith(files, exp_name)]
  file_path <- file.path(path, file_name)
  load(file_path)
  results <- results[1:50]
  results <- results[setdiff(names(results), "sample_data")]

  # process sensitivity results
  sens <- sapply(results, sapply, `[[`, "sens")
  sens_mean <- rowMeans(sens)
  max_sens_ind <- which.max(sens_mean)
  sens_mean <- sprintf(paste0("%.", prec, "f"), sens_mean)
  sens_mean[max_sens_ind] <- paste0("\\mathbf{", sens_mean[max_sens_ind], "}")
  sens_sd  <- sprintf(paste0("%.", prec, "f"), apply(sens, 1, sd))
  sens_str <- paste0(sens_mean, " (", sens_sd, ")")

  # process specificity results
  spec <- sapply(results, sapply, `[[`, "spec")
  spec_mean <- rowMeans(spec)
  max_spec_ind <- which.max(spec_mean)
  spec_mean <- sprintf(paste0("%.", prec, "f"), spec_mean)
  spec_mean[max_spec_ind] <- paste0("\\mathbf{", spec_mean[max_spec_ind], "}")
  spec_sd  <- sprintf(paste0("%.", prec, "f"), apply(spec, 1, sd))
  spec_str <- paste0(spec_mean, " (", spec_sd, ")")

  # combine summary strings
  perf_str <- cbind(sens_str, spec_str)
  perf_str <- matrix(paste0("$", perf_str, "$"), dim(perf_str))
  row.names(perf_str) <- row.names(spec)

  # create storage
  p_str <- paste0(c(p - 1, rep("!", length(df))), collapse = "")
  df_exp <- data.frame(p = p_str, method = c("covdepGE" , "covdepGE_sortZ", "loggle", "mgm"),
                       sens = NA, spec = NA)
  df_exp[ , c("sens", "spec")] <- perf_str[df_exp$method, ]
  df_exp$method <- unlist(methods[df_exp$method])

  df[[length(df) + 1]] <- df_exp
  rm("results")
}


df <- Reduce(rbind, df)
colnames(df) <- c("$p$", "Method", "Sensitivity$(\\uparrow)$", "Specificity$(\\uparrow)$")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")

# ------------------------------------------------------------------------------
# continuous covariate dependent visualization
rm(list = ls())
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study_graph_learning")
library(ggplot2)
library(latex2exp)

# load the covdep results for the p=11 case
exper <- "cont_cov_dep_"
trial_str <- "ntrials70_"
samp_str <- "_n1_50_n2_50_n3_50"
exp_name <- paste0(exper, trial_str, "p11", samp_str)
path <- "./experiments/continuous/covariate dependent"
files <- list.files(path)
file_name <- files[startsWith(files, exp_name)]
file_path <- file.path(path, file_name)
load(file_path)

# extract the pip for covdepGE
true <- results$sample_data$true_precision
results <- results[1:50]
results <- results[setdiff(names(results), "sample_data")]
results <- lapply(results, `[[`, "covdepGE")
results <- lapply(results, `[[`, "pip")

# extract pip for 1, 2 and 1, 3 entries
n_trial <- length(results)
pip12 <- sapply(results, sapply, `[`, 1, 2)
pip13 <- sapply(results, sapply, `[`, 1, 3)

# add all lines to plots
plot12 <- ggplot(data.frame(Individual = 1:nrow(pip12)), aes(Individual)) + theme_classic() + ggtitle(TeX("$\\alpha_{12}$ across 50 trials")) + theme(plot.title = element_text(hjust = 0.5))
mean12 <- rowMeans(pip12)
upper12 <- apply(pip12, 1, quantile, 0.95)
lower12 <- apply(pip12, 1, quantile, 0.05)
for (j in 1:ncol(pip12)){
  plot12 <- plot12 + geom_line(data = data.frame(Individual = 1:nrow(pip12), PIP = pip12[ , j]), aes(Individual, PIP), alpha = 0.05)
}
plot12 <- plot12 +
  geom_line(aes(y = mean12), color = "tomato2", size = 1) +
  geom_line(aes(y = lower12), linetype = "dashed", color = "navy") +
  geom_line(aes(y = upper12), linetype = "dashed", color = "navy")

plot13 <- ggplot(data.frame(Individual = 1:nrow(pip13)), aes(Individual)) + theme_classic() + ggtitle(TeX("$\\alpha_{13}$ across 50 trials")) + theme(plot.title = element_text(hjust = 0.5))
mean13 <- rowMeans(pip13)
upper13 <- apply(pip13, 1, quantile, 0.95)
lower13 <- apply(pip13, 1, quantile, 0.05)
for (j in 1:ncol(pip13)){
  plot13 <- plot13 + geom_line(data = data.frame(Individual = 1:nrow(pip13), PIP = pip13[ , j]), aes(Individual, PIP), alpha = 0.05)
}
plot13 <- plot13 +
  geom_line(aes(y = mean13), color = "tomato2", size = 1) +
  geom_line(aes(y = lower13), linetype = "dashed", color = "navy") +
  geom_line(aes(y = upper13), linetype = "dashed", color = "navy")

plot12
plot13

# ------------------------------------------------------------------------------
# cdge analysis
rm(list = ls())
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study_graph_learning")
library(kableExtra)

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each
exper <- "covdepGE_disc_cov_dep_"
trial_str <- "ntrials50_"
dims <- c(11, 31, 51)
lambdas <- c(3, 15)
n1s <- c(50, 80)
path <- "./experiments/cdge2" # cdge is using the hp from varbvs, cdge2 does not use
files <- list.files(path)
prec <- 4
df <- list()
for (p in dims){
  for (lambda in lambdas){
    for (n1 in n1s){

      # format file name
      n2 <- 100 - n1
      exp_name <- paste0(exper, trial_str, "p", p, "_n1_", n1, "_n2_", n2, "_lambda", lambda)

      # load results
      file_name <- files[startsWith(files, exp_name)]
      if (length(file_name) == 0) next
      file_path <- file.path(path, file_name)
      load(file_path)
      results <- results[setdiff(names(results), "sample_data")]
      results <- lapply(results, `[[`, "covdepGE")
      if (any(is.na(results))) next

      # process sensitivity results
      sens <- sapply(results, `[[`, "sens")
      sens_mean <- mean(sens)
      sens_mean <- sprintf(paste0("%.", prec, "f"), sens_mean)
      sens_sd  <- sprintf(paste0("%.", prec, "f"), sd(sens))
      sens_str <- paste0(sens_mean, " (", sens_sd, ")")

      # process specificity results
      spec <- sapply(results, `[[`, "spec")
      spec_mean <- mean(spec)
      spec_mean <- sprintf(paste0("%.", prec, "f"), spec_mean)
      spec_sd  <- sprintf(paste0("%.", prec, "f"), sd(spec))
      spec_str <- paste0(spec_mean, " (", spec_sd, ")")

      # combine summary strings
      perf_str <- c(sens_str, spec_str)

      # create storage
      df_exp <- data.frame(p = p, c = lambda, n1 = n1, n2 = n2,
                           method = "covdepGE", sens = NA, spec = NA)
      df_exp[ , c("sens", "spec")] <- perf_str

      df[[length(df) + 1]] <- df_exp
      rm("results")
    }
  }
}

df <- Reduce(rbind, df)
colnames(df) <- c("$p$", "$c$", "$n_1$", "$n_2$", "Method", "Sensitivity$(\\uparrow)$", "Specificity$(\\uparrow)$")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")
