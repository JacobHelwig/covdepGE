# ------------------------------------------------------------------------------
# Baseline comparisons
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(kableExtra)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(latex2exp)

subgroups_list <- list()

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each

# results config; med outputs median in place of mean and univariate is for
# switching setting from q=1->q=2
med <- F
univariate <- !F
sine <- F
four <- !T
seq <- F
if (univariate){

  # univariate extraneous covariate
  exper <- "cont_cov_dep_"
  methods <- c("covdepGE" , "JGL", "mgm")
  n <- 75
  samp_str <- paste0("_n1_", n, "_n2_", n, "_n3_", n)
  path <- "./experiments/z1"
  if (sine){
    path <- paste0(path, "_sine")
    exper <- paste0(exper, "sine_")
  }else if (seq){
    path <- paste0(path, "_seq")
    methods <- c("covdepGE", "covdepGE_seq")
  }
}else{

  # multivariate extraneous covariate
  exper <- ifelse(four, "cont_4_cov_dep_", "cont_multi_cov_dep_")
  methods <- c("covdepGE" , "JGL", "mgm", "covdepGE_sortZ")
  n <- 25 + 200 * four
  samp_str <- paste0("_n_", n)
  path <- ifelse(four, "./experiments/z4", "./experiments/z2")
}

ntrials <- 50
trial_str <- paste0("ntrials", ntrials, "_")
dims <- c(10, 25, 50, 100)
files <- list.files(path)
prec <- 2
df <- subgroups <- htest <- ngraphs <- vector("list", length(dims))
names(df) <- names(subgroups) <- names(htest) <- names(ngraphs) <- dims

for (p in as.character(dims)){

  # format file name
  exp_name <- paste0(exper, trial_str, "p", p, samp_str)

  # load results
  # mean(sapply(results, function(trial)length(trial$covdepGE$graphs$unique_graphs)))
  file_name <- files[grepl(exp_name, files) & endsWith(files, ".Rda")]
  if (length(file_name) != 1) stop(paste(length(file_name), "files found for", exp_name))
  file_path <- file.path(path, file_name)
  load(file_path)
  print(paste0("n: ", results$sample_data[1], ", p: ", results$sample_data[2]))
  results <- results[setdiff(names(results), "sample_data")]

  # extract number of graphs
  ngraphs[[p]] <- sapply(results, function(trial)length(trial$covdepGE$graphs$unique_graphs))

  # out <- results$trial2$covdepGE
  # for (j in 1:length(out$graphs$unique_graphs)){
  #   out$graphs$unique_graphs[[j]]$graph <- as.matrix(out$graphs$unique_graphs[[j]]$graph)
  # }
  # plot(out)

  # process sensitivity results
  sens <- sapply(results, sapply, `[[`, "sens")
  sens_mean <- rowMeans(sens)
  if (med) sens_mean <- apply(sens, 1, median)
  max_sens_ind <- which.max(sens_mean)
  sens_mean <- sprintf(paste0("%.", prec, "f"), sens_mean * 100)
  sens_mean[max_sens_ind] <- paste0("\\mathbf{", sens_mean[max_sens_ind], "}")
  sens_sd  <- sprintf(paste0("%.", prec, "f"), apply(sens * 100, 1, sd))
  if (med) sens_sd <- sprintf(paste0("%.", prec, "f"), apply(sens * 100, 1, IQR))
  sens_str <- paste0(sens_mean, " (", sens_sd, ")")

  # hypothesis testing
  if (!seq){
    t_test <- t.test(sens['covdepGE',], sens['JGL',], alternative="greater")
    pval <- sprintf(paste0("%.", prec + 1, "f"), t_test$p.value)
    pval <- ifelse(round(t_test$p.value,prec+1) == 0, paste0(c("<0.",rep(0, prec), 1), collapse=""), pval)
    # pval <- ifelse(t_test$p.value < 0.05, paste0("\\mathbf{", pval, "}"), pval)
    pval <- paste0("$", pval, "$")
    stat <- paste0("$", sprintf(paste0("%.", prec, "f"), t_test$statistic), "$")
    degf <- paste0("$", sprintf(paste0("%.", prec, "f"), t_test$parameter), "$")
    mean_covdepGE <- paste0("$", sprintf(paste0("%.", prec, "f"), t_test$estimate[1] * 100), "$")
    # mean_JGL <- paste0("$", sprintf(paste0("%.", prec, "f"), t_test$estimate[2] * 100), "$")
    se <- paste0("$", sprintf(paste0("%.", prec, "f"), t_test$stderr * 100), "$")
    htest[[p]] <- data.frame(covdepGE=paste0("$", sens_str[row.names(sens) == "covdepGE"], "$"),
                             JGL=paste0("$", sens_str[row.names(sens) == "JGL"], "$"),
                             stderr=se, stat=stat, df=degf, p=pval)
  }

  # process specificity results
  spec <- sapply(results, sapply, `[[`, "spec")
  spec_mean <- rowMeans(spec)
  if (med) spec_mean <- apply(spec, 1, median)
  max_spec_ind <- which.max(spec_mean)
  spec_mean <- sprintf(paste0("%.", prec, "f"), spec_mean * 100)
  spec_mean[max_spec_ind] <- paste0("\\mathbf{", spec_mean[max_spec_ind], "}")
  spec_sd  <- sprintf(paste0("%.", prec, "f"), apply(spec * 100, 1, sd))
  if (med) spec_sd <- sprintf(paste0("%.", prec, "f"), apply(spec * 100, 1, IQR))
  spec_str <- paste0(spec_mean, " (", spec_sd, ")")

  # process time results
  time <- sapply(results, sapply, `[[`, "time")
  time_mean <- rowMeans(time)
  if (seq){
    time_ratios <- time['covdepGE_seq',] / time['covdepGE',]
    ratio_mean <- sprintf(paste0("%.", prec, "f"), mean(time_ratios))
    ratio_sd <- sprintf(paste0("%.", prec, "f"), sd(time_ratios))
    ratio_str <- paste0(ratio_mean, " (", ratio_sd, ")")
  }
  if (med) time_mean <- apply(time, 1, median)
  min_time_ind <- which.min(time_mean)
  time_mean <- sprintf(paste0("%.", prec, "f"), time_mean)
  time_mean[min_time_ind] <- paste0("\\mathbf{", time_mean[min_time_ind], "}")
  time_sd  <- sprintf(paste0("%.", prec, "f"), apply(time, 1, sd))
  if (med) time_sd <- sprintf(paste0("%.", prec, "f"), apply(time, 1, IQR))
  time_str <- paste0(time_mean, " (", time_sd, ")")

  # combine summary strings
  perf_str <- cbind(sens_str, spec_str, time_str)
  perf_str <- matrix(paste0("$", perf_str, "$"), dim(perf_str))
  row.names(perf_str) <- row.names(spec)

  # create storage
  df_exp <- data.frame(p = p, method = methods,
                       sens = NA, spec = NA, time = NA)
  df_exp[ , c("sens", "spec", "time")] <- perf_str[df_exp$method, ]
  if (seq){
    df_exp$time_ratio <- ratio_str
  }
  df[[p]] <- df_exp

  subgroups[[p]] <- sapply(results, function(trial) length(unique(trial$JGL$classification)))

  rm("results")
}

if (seq){
  df <- Reduce(rbind, df)
  df <- df[,c("p", "method", "time", "time_ratio")]
  df <- data.frame(p = df[df$method == "covdepGE", "p"],
                   ptime = df[df$method == "covdepGE","time"],
                   stime = df[df$method == "covdepGE_seq","time"],
                   ratio = df[df$method == "covdepGE_seq","time_ratio"])
  colnames(df) <- colnames(df) <- c("$p$", "Parallel Time (seconds)", "Sequential Time (seconds)", "Ratio")
  kbl(df, format = "latex", booktabs = T, escape = FALSE) # %>%
    collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")
}

htest_df <- Reduce(rbind, htest)
colnames(htest_df) <- c("\\texttt{covdepGE}", "\\texttt{JGL}", "Standard Error", "T-statistic", "DF", "p-value")
kbl(htest_df, format = "latex", booktabs = T, escape = FALSE)

df <- Reduce(rbind, df)
df$method <- gsub("_sortZ", "\\\\_time", df$method)
df$method <- paste0("\\texttt{", df$method, "}")
colnames(df) <- c("$p$", "Method", "Sensitivity$(\\%)$", "Specificity$(\\%)$", "Time (seconds)")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")

subgroups_list[[exper]] <- factor(Reduce(c, subgroups), levels = 2:6)

windowsFonts("Times" = windowsFont("Times"))
colors <- c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF")
plots[[exper]] <- list(x = factor(subgroups, levels = 2:6), plot = NULL)
plots <- list(ggplot() +
                geom_bar(aes(x = subgroups_list[["cont_cov_dep_"]]),
                         color = "black", fill = colors[1]) +
                theme_pubclean() +
                theme(text = element_text(family = "Times", size = 18),
                      plot.title = element_text(hjust = 0.5)) +
                ggtitle(TeX(paste0("Optimal $\\textit{K}, \\textit{q}=1$"))) +
                scale_y_continuous(breaks = seq(0, 180, 30), limits = c(0, 180)) +
                scale_x_discrete(drop=FALSE) +
                labs(x = TeX("$\\textit{K}$")),
              ggplot() +
                geom_bar(aes(x = subgroups_list[["cont_multi_cov_dep_"]]),
                         color = "black", fill = colors[2]) +
                theme_pubclean() +
                theme(text = element_text(family = "Times", size = 18),
                      plot.title = element_text(hjust = 0.5)) +
                ggtitle(TeX(paste0("Optimal $\\textit{K}, \\textit{q}=2$"))) +
                scale_y_continuous(breaks = seq(0, 180, 30), limits = c(0, 180)) +
                scale_x_discrete(drop=FALSE) +
                labs(x = TeX("$\\textit{K}$")))
fig <- ggarrange(plotlist = plots, nrow = 2)
ggsave("plots/subgroup_bar_newcol.pdf", height = 6, width = 11)

# ------------------------------------------------------------------------------
# HP comparisons
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(kableExtra)

# load and store results from each experiment; get sensitivity and specificity
# and calculate mean/ sd for each

# results config; switches from comparisons with max_iter_grid=max_iter=100 to
# max_iter_grid=10
max_iter_grid <- F

if (max_iter_grid){
  # comparisons with max_iter_grid = max_iter (set max_iter_grid = T)
  exper <- "cont_cov_dep_"
  n <- 75
  samp_str <- paste0("_n1_", n, "_n2_", n, "_n3_", n)
  path <- "./experiments/z1_hp"
  files <- list.files(path)
  files <- files[grepl("model_average", files)]
  files <- c(files, list.files(paste0(path, "_full")))
}else{
  # default comparisons
  exper <- "cont_cov_dep_"
  n <- 75
  samp_str <- paste0("_n1_", n, "_n2_", n, "_n3_", n)
  path <- "./experiments/z1_hp"
  files <- c(list.files(path), list.files(substr(path, 1, nchar(path) - 3)))
}
ntrials <- 50
trial_str <- paste0("ntrials", ntrials, "_")
dims <- c(10, 25, 50, 100)
prec <- 2

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
      if (max_iter_grid){
        file_path <- file.path(paste0(path, "_full"), file_name)
      }else{
        file_path <- file.path(substr(path, 1, nchar(path) - 3), file_name)
      }
    }
    load(file_path)
    print(paste0("n: ", results$sample_data[1], ", p: ", results$sample_data[2]))
    results <- results[setdiff(names(results), "sample_data")]

    results <- lapply(results, `[`, "covdepGE")

    # process sensitivity results
    sens <- sapply(results, sapply, `[[`, "sens")
    sens_mean <- mean(sens)
    if (sens_mean > best_sens){
      best_sens <- sens_mean
      sens_method <- hp_method
    }
    sens_mean <- sprintf(paste0("%.", prec, "f"), sens_mean * 100)
    sens_sd  <- sprintf(paste0("%.", prec, "f"), sd(sens * 100))
    sens_str <- paste0(sens_mean, " (", sens_sd, ")")

    # process specificity results
    spec <- sapply(results, sapply, `[[`, "spec")
    spec_mean <- mean(spec)
    if (spec_mean > best_spec){
      best_spec <- spec_mean
      spec_method <- hp_method
    }
    spec_mean <- sprintf(paste0("%.", prec, "f"), spec_mean * 100)
    spec_sd  <- sprintf(paste0("%.", prec, "f"), sd(spec * 100))
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
    df_exp <- data.frame(p = p, method = hp_method,
                         sens = NA, spec = NA, time = NA)

    df_exp[ , c("sens", "spec", "time")] <- perf_str

    df[[p]][[hp_method]] <- df_exp
    rm("results")
  }
  best_sens <- sprintf(paste0("%.", prec, "f"), best_sens * 100)
  best_spec <- sprintf(paste0("%.", prec, "f"), best_spec * 100)
  best_time <- sprintf(paste0("%.", prec, "f"), best_time)
  df[[p]][[sens_method]]$sens <- gsub(best_sens, paste0("\\\\mathbf{", best_sens, "}"), df[[p]][[sens_method]]$sens)
  df[[p]][[spec_method]]$spec <- gsub(best_spec, paste0("\\\\mathbf{", best_spec, "}"), df[[p]][[spec_method]]$spec)
  df[[p]][[time_method]]$time <- gsub(best_time, paste0("\\\\mathbf{", best_time, "}"), df[[p]][[time_method]]$time)
}

df <- lapply(df, function(x) Reduce(rbind, x))
df <- Reduce(rbind, df)
df$method <- gsub("_", "\\\\_", df$method)
df$method <- paste0("\\texttt{", df$method, "}")
colnames(df) <- c("$p$", "HP Method", "Sensitivity$(\\%)$", "Specificity$(\\%)$", "Time (seconds)")
kbl(df, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2, 3, 4), latex_hline = "major", valign = "middle")
