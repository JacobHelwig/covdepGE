# visualizing the true graphs
library(covdepGE)
library(extrafont)
library(ggplot2)
library(ggpubr)
library(latex2exp)
loadfonts(device = "win")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("data.R")

# z1
dat <- cont_cov_dep_data(p = 10, n1 = 1, n2 = 1, n3 = 1)
true_graphs <- lapply(lapply(lapply(dat$true_precision, `!=`, 0), `*`, 1), `-`, diag(10))
graph_inds <- c("$z_l\\in\\left[-3,-1\\right]$", "$z_l\\in\\left[-1,1\\right]$", "$z_l\\in\\left[-3,-1\\right]$")
titles <- lapply(graph_inds, TeX)
windowsFonts("Times" = windowsFont("Times"))
true_list <- lapply(1:length(true_graphs), function(k)
  matViz(true_graphs[[k]]) +
    ggtitle(titles[[k]]) +
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)) +
    theme(plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          text = element_text(family = "Times")
          ))
(true_graph <- ggarrange(plotlist = true_list, nrow = 1, common.legend = T))
ggsave("plots/ground_truth_CDS_z1.pdf", true_graph, height = 4, width = 11)

# z2
dat <- cont_multi_cov_dep_data(p = 10, n = 1)
true_graphs <- lapply(lapply(lapply(dat$true_precision, `!=`, 0), `*`, 1), `-`, diag(10))
graph_inds <- c("$z_l\\in\\left[-3,-1\\right]$", "$z_l\\in\\left[-1,1\\right]$", "$z_l\\in\\left[-3,-1\\right]$")
titles <- lapply(graph_inds, TeX)
windowsFonts("Times" = windowsFont("Times"))
true_list <- lapply(1:length(true_graphs), function(k)
  matViz(true_graphs[[k]]) +
    ggtitle(titles[[k]]) +
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)) +
    theme(plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          text = element_text(family = "Times")
    ))
(true_graph <- ggarrange(plotlist = true_list, nrow = 1, common.legend = T))
ggsave("plots/ground_truth_CDS_z1.pdf", true_graph, height = 4, width = 11)
