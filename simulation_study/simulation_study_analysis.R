rm(list = ls())
library(covdepGE)
library(JGL)
library(mgm)
library(kableExtra)

# p = 25, n = 150
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/res_p25_n150_20220820_205707.Rda")
results25_150 <- results

# p = 50, n = 150
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/res_p50_n150_20220820_205746.Rda")
results50_150 <- results

# p = 100, n = 300
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/res_p100_n300_20220820_210016.Rda")
results100_300 <- results

# p = 100, n = 600
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/res_p100_n600_20220820_205829.Rda")
results100_600 <- results

# delete this
results100_600 <- results100_300 <- results50_150

# remove NULLS
results25_150 <- results25_150[!sapply(results25_150, is.null)]
results50_150 <- results50_150[!sapply(results50_150, is.null)]
results100_300 <- results100_300[!sapply(results100_300, is.null)]
results100_600 <- results100_600[!sapply(results100_600, is.null)]

# put all results together
results <- list(
  p25_n150 = results25_150,
  p50_n150 = results50_150,
  p100_n300 = results100_300,
  p100_n600 = results100_600
)

# aggregate results by model and put into a models list
covdepGE_mods <- lapply(results, lapply, `[[`, "covdepGE")
jgl_mods <- lapply(results, lapply, `[[`, "JGL")
mgm_mods <- lapply(results, lapply, `[[`, "mgm")
mods <- list(covdepGE = covdepGE_mods,
             JGL = jgl_mods,
             mgm = mgm_mods)

# extract times, sensitivity, specificity, ect.
times <- lapply(mods, lapply, sapply, `[[`, "time")
sens <- lapply(mods, lapply, sapply, `[[`, "sens")
spec <- lapply(mods, lapply, sapply, `[[`, "spec")
TP <- lapply(mods, lapply, sapply, `[[`, "TP_n")
# TN <- lapply(mods, lapply, sapply, `[[`, "TN_n")
FP <- lapply(mods, lapply, sapply, `[[`, "FP_n")
# FN <- lapply(mods, lapply, sapply, `[[`, "FN_n")

# function to get the mean and standard deviation at a specified precision
mean_sd <- function(x, prec = 2, mean_format = "f", sd_format = "f") paste0(
  formatC(mean(x), prec, format = mean_format), "(",
  formatC(sd(x), prec, format = sd_format), ")")

# get summary stats for each
times_sum <- lapply(times, lapply, mean_sd)
sens_sum <- lapply(sens, lapply, mean_sd)
spec_sum <- lapply(spec, lapply, mean_sd, sd_format = "g")
TP_sum <- lapply(TP, lapply, mean_sd)
FP_sum <- lapply(FP, lapply, mean_sd)

# re-aggregate from by model to by experiment
times_exp <- lapply(names(results), function(exp_name) sapply(
  times_sum, `[[`, exp_name))
sens_exp <- lapply(names(results), function(exp_name) sapply(
  sens_sum, `[[`, exp_name))
spec_exp <- lapply(names(results), function(exp_name) sapply(
  spec_sum, `[[`, exp_name))
TP_exp <- lapply(names(results), function(exp_name) sapply(
  TP_sum, `[[`, exp_name))
FP_exp <- lapply(names(results), function(exp_name) sapply(
  FP_sum, `[[`, exp_name))
names(times_exp) <- names(sens_exp) <- names(spec_exp) <- names(TP_exp) <-
  names(FP_exp) <- names(results)

# create a matrix for each experiment
exp_sum <- list(Sensitivity = sens_exp, Specificity = spec_exp,
                `TP/graph` = TP_exp, `FP/graph` = FP_exp,
                `Time(s)` = times_exp)
p25_n150df <- as.matrix(data.frame(lapply(exp_sum, `[[`, "p25_n150")))
p50_n150df <- as.matrix(data.frame(lapply(exp_sum, `[[`, "p50_n150")))
p100_n300df <- as.matrix(data.frame(lapply(exp_sum, `[[`, "p100_n300")))
p100_n600df <- as.matrix(data.frame(lapply(exp_sum, `[[`, "p100_n600")))

# combine all of the matrices
res_mat <- rbind(p25_n150df, p50_n150df, p100_n300df, p100_n600df)
row.names(res_mat) <- rep(row.names(p25_n150df), 4)
colnames(res_mat) <- names(exp_sum)

kbl(res_mat, format = "latex", booktabs = T)

times_df <- data.frame(lapply(times_sum, unlist))
sens_df <- data.frame(lapply(sens_sum, unlist))

list(times = times_df, sens = sens_df)


unlist(times_sum$covdepGE)

perf_df <- cbind.data.frame(
  Sensitivity = sapply(sens, mean_sd),
  Specificity = sapply(spec, function(x)
    paste0(round(mean(x), 4), "(", formatC(sd(x), 2, format = "e"), ")")),
  `FP per graph` = sapply(FP, mean_sd),
  `FN per graph` = sapply(FN, mean_sd))
library(kableExtra)
kbl(perf_df, format = "latex", booktabs = T)

