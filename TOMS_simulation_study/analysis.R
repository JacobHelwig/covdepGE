# ------------------------------------------------------------------------------
# continuous covariate dependent analysis
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(kableExtra)

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each

# univariate extraneous covariate
exper <- "cont_cov_dep_"
n <- 75
samp_str <- paste0("_n1_", n, "_n2_", n, "_n3_", n)
path <- "./experiments/z1"

# multivariate extraneous covariate
exper <- "cont_multi_cov_dep_"
n <- 25
samp_str <- paste0("_n_", n)
path <- "./experiments/z2"

ntrials <- 50
trial_str <- paste0("ntrials", ntrials, "_")
dims <- c(10, 25, 50)
files <- list.files(path)
prec <- 4
df <- vector("list", length(dims))
names(df) <- dims
for (j in 1:length(dims)){

  p <- dims[j]

  # format file name
  exp_name <- paste0(exper, trial_str, "p", p, samp_str)

  # load results
  file_name <- files[grepl(exp_name, files) & endsWith(files, ".Rda")]
  if (length(file_name) != 1) stop(paste(length(file_name), "files found"))
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

  # process time results
  time <- sapply(results, sapply, `[[`, "time")
  time_mean <- rowMeans(time)
  min_time_ind <- which.min(time_mean)
  time_mean <- sprintf(paste0("%.", prec, "f"), time_mean)
  time_mean[min_time_ind] <- paste0("\\mathbf{", time_mean[min_time_ind], "}")
  time_sd  <- sprintf(paste0("%.", prec, "f"), apply(time, 1, sd))
  time_str <- paste0(time_mean, " (", time_sd, ")")

  # combine summary strings
  perf_str <- cbind(sens_str, spec_str, time_str)
  perf_str <- matrix(paste0("$", perf_str, "$"), dim(perf_str))
  row.names(perf_str) <- row.names(spec)

  # create storage
  p_str <- paste0(c(p, rep("!", which(p == dims))), collapse = "")
  df_exp <- data.frame(p = p_str, method = c("covdepGE" , "JGL", "mgm"),
                       sens = NA, spec = NA, time = NA)
  df_exp[ , c("sens", "spec", "time")] <- perf_str[df_exp$method, ]
  df[[p]] <- df_exp

  rm("results")
}

df <- Reduce(rbind, df)
colnames(df) <- c("$p$", "Method", "Sensitivity$(\\uparrow)$", "Specificity$(\\uparrow)$", "Time (seconds)")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")

# ------------------------------------------------------------------------------
# continuous covariate dependent analysis
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(kableExtra)

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each

# univariate extraneous covariate
exper <- "cont_cov_dep_"
n <- 75
samp_str <- paste0("_n1_", n, "_n2_", n, "_n3_", n)
path <- "./experiments/z1_hp"

ntrials <- 50
trial_str <- paste0("ntrials", ntrials, "_")
dims <- c(10, 25, 50)
files <- c(list.files(path), list.files(substr(path, 1, nchar(path) - 3)))
prec <- 4

df <- list(grid_search = NA, hybrid = NA, model_average = NA)
df <- replicate(length(dims), df, simplify = F)
names(df) <- dims

for (p in as.character(dims)){

  best_sens <- best_spec <- 0
  best_time <- Inf
  sens_method <- spec_method <- time_method <- NA

  for (hp_method in names(df[[p]])){

    # format file name
    exp_name <- paste0(exper, trial_str, "p", p, samp_str)

    # load results

    file_name <- files[grepl(exp_name, files) & grepl(hp_method, files) & endsWith(files, ".Rda")]

    if (length(file_name) != 1) stop(paste(length(file_name), "files found"))
    file_path <- file.path(path, file_name)
    if (!file.exists(file_path)){
      file_path <- file.path(substr(path, 1, nchar(path) - 3), file_name)
    }
    load(file_path)
    results <- results[setdiff(names(results), "sample_data")]

    results <- lapply(results, `[`, "covdepGE")

    # process sensitivity results
    sens <- sapply(results, sapply, `[[`, "sens")
    sens_mean <- mean(sens)
    if (sens_mean > best_sens){
      best_sens <- sens_mean
      sens_method <- hp_method
    }
    sens_mean <- sprintf(paste0("%.", prec, "f"), sens_mean)
    sens_sd  <- sprintf(paste0("%.", prec, "f"), sd(sens))
    sens_str <- paste0(sens_mean, " (", sens_sd, ")")

    # process specificity results
    spec <- sapply(results, sapply, `[[`, "spec")
    spec_mean <- mean(spec)
    if (spec_mean > best_spec){
      best_spec <- spec_mean
      spec_method <- hp_method
    }
    spec_mean <- sprintf(paste0("%.", prec, "f"), spec_mean)
    spec_sd  <- sprintf(paste0("%.", prec, "f"), sd(spec))
    spec_str <- paste0(spec_mean, " (", spec_sd, ")")

    # process time results
    time <- sapply(results, sapply, `[[`, "time")
    time_mean <- mean(time)
    if (time_mean < best_time){
      best_time <- time_mean
      time_method <- hp_method
    }
    time_mean <- sprintf(paste0("%.", prec, "f"), time_mean)
    time_sd  <- sprintf(paste0("%.", prec, "f"), sd(time))
    time_str <- paste0(time_mean, " (", time_sd, ")")

    # combine summary strings
    perf_str <- cbind(sens_str, spec_str, time_str)
    perf_str <- matrix(paste0("$", perf_str, "$"), dim(perf_str))
    row.names(perf_str) <- row.names(spec)

    # create storage
    p_str <- paste0(c(p, rep("!", which(p == dims))), collapse = "")

    df_exp <- data.frame(p = p_str, method = hp_method,
                         sens = NA, spec = NA, time = NA)

    df_exp[ , c("sens", "spec", "time")] <- perf_str

    df[[p]][[hp_method]] <- df_exp
    rm("results")
  }
  best_sens <- sprintf(paste0("%.", prec, "f"), best_sens)
  best_spec <- sprintf(paste0("%.", prec, "f"), best_spec)
  best_time <- sprintf(paste0("%.", prec, "f"), best_time)
  df[[p]][[sens_method]]$sens <- gsub(best_sens, paste0("\\\\mathbf{", best_sens, "}"), df[[p]][[sens_method]]$sens)
  df[[p]][[spec_method]]$spec <- gsub(best_spec, paste0("\\\\mathbf{", best_spec, "}"), df[[p]][[spec_method]]$spec)
  df[[p]][[time_method]]$time <- gsub(best_time, paste0("\\\\mathbf{", best_time, "}"), df[[p]][[time_method]]$time)
}

df <- lapply(df, function(x) Reduce(rbind, x))
df <- Reduce(rbind, df)
colnames(df) <- c("$p$", "Method", "Sensitivity$(\\uparrow)$", "Specificity$(\\uparrow)$")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")
