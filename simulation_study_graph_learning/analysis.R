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
n1s <- c(50, 80)
path <- "./experiments/discrete/covariate free"
files <- list.files(path)
prec <- 4
df <- list()
res <- list()
methods <- list(covdepGE = "W-PL", varbvs = "CS")
for (lambda in lambdas){
  for (n1 in n1s){

    if (n1 == 80) next

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
                         method = c("covdepGE", "varbvs"), sens = NA,
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

load("disc_cov_dep_ntrials3_p11_n1_50_n2_50_lambda15_20221123_155056.Rda")
#load("disc_cov_dep_ntrials3_p11_n1_50_n2_50_lambda15_20221123_161228.Rda")

res <- results[setdiff(names(results), "sample_data")]
res <- lapply(res, lapply, `[`, c("sens", "spec"))
sens <- sapply(res, sapply, `[[`, "sens")
spec <- sapply(res, sapply, `[[`, "spec")
rowMeans(sens)
rowMeans(spec)
apply(sens, 1, sd)
apply(spec, 1, sd)


