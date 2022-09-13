rm(list = ls())
library(covdepGE)
library(JGL)
library(kableExtra)
library(mgm)

# function to add covdepGE results to reuslts from JGM and mgm
add_res <- function(covdepGE_res_path, JGM_mgm_res_path){

  # load the JGM and mgm results and save; do the same for covdepGE
  load(JGM_mgm_res_path)
  results_final <- results
  load(covdepGE_res_path)
  names(results_final) <- names(results)

  # add covdepGE results to JGM and mgm results and return
  for(j in 1:length(results)){
    results_final[[j]]$covdepGE <- results[[j]]
  }
  results_final
}

# p = 5, n = 90
results5_90 <- add_res(
  "~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/p5_n90/res_p5_n90_covdepGE_20220908_215120.Rda",
  "~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/p5_n90/res_p5_n90_JGL_mgm_20220908_215316.Rda")

# p = 15, n = 90
results15_90 <- add_res(
  "~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/p15_n90/res_p15_n90_covdepGE_20220908_215229.Rda",
  "~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/p15_n90/res_p15_n90_JGL_mgm_20220908_215604.Rda")

# p = 25, n = 150
results25_150 <- add_res(
  "~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/p25_n150/res_p25_n150_covdepGE_20220825_121750.Rda",
  "~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/p25_n150/res_p25_n150_JGL_mgm_20220825_090115.Rda")

# p = 50, n = 150
results50_150 <- add_res(
  "~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/p50_n150/res_p50_n150_covdepGE_20220825_090326.Rda",
  "~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/p50_n150/res_p50_n150_JGL_mgm_20220825_090205.Rda"
)

# p = 100, n = 300
results100_300 <- add_res(
  "~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/p100_n300/res_p100_n300_covdepGE_20220824_084919.Rda",
  "~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/p100_n300/res_p100_n300_JGL_mgm_20220823_102037.Rda"
)

# put all results together
results <- list(
  p5_n90 = results5_90,
  p15_n90 = results15_90,
  p25_n150 = results25_150,
  p50_n150 = results50_150,
  p100_n300 = results100_300
)
sapply(results, length)

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
mean_sd <- function(x, prec = 4, mean_format = "f", sd_format = "f") {
  paste0("$", formatC(mean(x), prec, format = mean_format), "(",
         formatC(sd(x), prec, format = sd_format), ")$")
}

# get summary stats for each
times_sum <- lapply(times, lapply, mean_sd)
sens_sum <- lapply(sens, lapply, mean_sd)
spec_sum <- lapply(spec, lapply, mean_sd)
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
exp_sum <- list("Sensitivity$(\\uparrow)$" = sens_exp,
                "Specificity$(\\uparrow)$" = spec_exp,
                # "TP/graph$(\\uparrow)$" = TP_exp,
                # "FP/graph$(\\downarrow)$" = FP_exp,
                "Time(s)$(\\downarrow)$" = times_exp)
p5_n90df <- as.matrix(data.frame(lapply(exp_sum, `[[`, "p5_n90")))
p15_n90df <- as.matrix(data.frame(lapply(exp_sum, `[[`, "p15_n90")))
p25_n150df <- as.matrix(data.frame(lapply(exp_sum, `[[`, "p25_n150")))
p50_n150df <- as.matrix(data.frame(lapply(exp_sum, `[[`, "p50_n150")))
p100_n300df <- as.matrix(data.frame(lapply(exp_sum, `[[`, "p100_n300")))

# combine all of the matrices
n <- rep(c(90, "90_", 150, "150_", 300), each = 3)
p <- rep(c(5, 15, 25, 50, 100), each = 3)
pkg_names <- rep(paste0("\\texttt{", row.names(p25_n150df), "}"), 5)
res_mat <- rbind(p5_n90df, p15_n90df, p25_n150df, p50_n150df, p100_n300df)
res_mat <- cbind(p = p, n = n, pkg = pkg_names, res_mat)
colnames(res_mat) <- c("$p$", "$n$", "Package", names(exp_sum))
rownames(res_mat) <- NULL

kbl(res_mat, format = "latex", booktabs = T, escape = FALSE) %>%
  collapse_rows(columns = c(1, 2), latex_hline = "major", valign = "middle")
