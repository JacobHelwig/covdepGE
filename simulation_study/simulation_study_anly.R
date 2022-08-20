# n = 75
# /home/jacob.a.helwig/covdepGE/simulation_study/results75_20220819_124138.Rda
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/results75_20220819_124138.Rda")
dim(results$trial1$data$X)

# n = 150
# /home/jacob.a.helwig/covdepGE/simulation_study/results150_20220819_124210.Rda
# /home/jacob.a.helwig/covdepGE/simulation_study/results150_20220819_205227.Rda
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/results150_20220819_124210.Rda")

# n = 300
# /home/jacob.a.helwig/covdepGE/simulation_study/results300_20220819_124221.Rda
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/results300_20220819_124221.Rda")

# n = 600
# /home/jacob.a.helwig/covdepGE/simulation_study/results600_20220819_124539.Rda
# /home/jacob.a.helwig/covdepGE/simulation_study/results600_20220819_203607.Rda
load("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/simulation_study/results600_20220819_124539.Rda")

# extract models
covdepGE_mods <- lapply(results, `[[`, "covdepGE")
jgl_mods <- lapply(results, `[[`, "JGL")
mgm_mods <- lapply(results, `[[`, "mgm")
mods <- list(covdepGE = covdepGE_mods,
             JGL = jgl_mods,
             mgm = mgm_mods)

# function for extracting
extractor <- function(lst, name) unlist(sapply(lst, `[[`, name))

# function for evaluating estimated graphs compared to ground truth
eval_est <- function(est, true){

  # get true number of edges and non-edges
  num_edge <- sum(true, na.rm = T)
  num_non <- sum(true == 0, na.rm = T)

  # calculate sensitivity and specificity
  true_edge <- sum(est == 1 & true == 1, na.rm = T)
  false_edge <- sum(est == 1 & true == 0, na.rm = T)
  true_non <- sum(est == 0 & true == 0, na.rm = T)
  false_non <- sum(est == 0 & true == 1, na.rm = T)
  sens <- true_edge / num_edge
  spec <- true_non / num_non

  list(sens = sens, spec = spec, TP = true_edge, FP = false_edge,
       TN = true_non, FN = false_non)
}

# extract times, sensitivity, and specificity
times <- lapply(mods, extractor, "time")
sens <- lapply(mods, extractor, "sens")
spec <- lapply(mods, extractor, "spec")

# analyze
lapply(times, summary)
lapply(sens, summary)
lapply(spec, summary)
length(spec$covdepGE)

# filter out the nulls
covdepGE_mods <- covdepGE_mods[!sapply(covdepGE_mods, is.null)]
jgl_mods <- jgl_mods[!sapply(jgl_mods, is.null)]
mgm_mods <- mgm_mods[!sapply(mgm_mods, is.null)]

# fit the sens and spec metrics
ddim <- dim(results$trial1$data$X)
n <- ddim[1]
p <- ddim[2]
true_graphs <- array(unlist(results$trial1$data$true_precision), c(p, p, n))
true_graphs <- (true_graphs != 0) * 1 + replicate(n, diag(rep(NA, p)) * 1)
res_covdepGE <- lapply(lapply(covdepGE_mods, `[[`, "str"), eval_est, true_graphs)
res_jgl <- lapply(lapply(jgl_mods, `[[`, "str"), eval_est, true_graphs)
res_mgm <- lapply(lapply(mgm_mods, `[[`, "str"), eval_est, true_graphs)
mod_res <- list(covdepGE = res_covdepGE,
                JGL = res_jgl,
                mgm = res_mgm)

# extract sensitivity and specificity
sens <- lapply(mod_res, extractor, "sens")
spec <- lapply(mod_res, extractor, "spec")

# extract TP, TN, FP, FN
TP <- lapply(lapply(mod_res, extractor, "TP"), `/`, n)
TN <- lapply(lapply(mod_res, extractor, "TN"), `/`, n)
FP <- lapply(lapply(mod_res, extractor, "FP"), `/`, n)
FN <- lapply(lapply(mod_res, extractor, "FN"), `/`, n)

# group together the summary statistics
mean_sd <- function(x) paste0(round(mean(x), 2), "(", round(sd(x), 2), ")")

perf_df <- cbind.data.frame(
  Sensitivity = sapply(sens, mean_sd),
  Specificity = sapply(spec, function(x)
    paste0(round(mean(x), 4), "(", formatC(sd(x), 2, format = "e"), ")")),
  `FP per graph` = sapply(FP, mean_sd),
  `FN per graph` = sapply(FN, mean_sd))
library(kableExtra)
kbl(perf_df, format = "latex", booktabs = T)

# visualizing the true graphs
library(ggplot2)
library(covdepGE)
library(ggpubr)
library(extrafont)
loadfonts(device = "win")
dat <- generateData(100, 1, 1, 1)
true_graphs <- lapply(lapply(lapply(dat$true_precision, `!=`, 0), `*`, 1), `-`, diag(100))
graph_inds <- c("1,...,zzz", "zzz,...,zzz", "zzz,...,zzz")
titles <- paste0("CDS ", 1:3, ", observations ", graph_inds)
true_list <- lapply(1:length(true_graphs), function(k)
  matViz(true_graphs[[k]]) +
    ggtitle(titles[k]) +
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)) +
    theme(plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          text = element_text(family = "Times")))
true_graph <- ggarrange(plotlist = true_list, nrow = 1, common.legend = T)
# true_graph <- ggarrange(true_list[[1]], NULL,
#                         true_list[[2]], NULL,
#                         true_list[[3]], nrow = 1, common.legend = T, widths = c(1, 0.5, 1, 0.5, 1))
ggsave("plots/true.pdf", true_graph, height = 4, width = 11)
