# visualizing the true graphs
rm(list=ls())
library(covdepGE)
library(extrafont)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(stringr)
loadfonts(device = "win")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("data.R")

# z1
title_font <- 36
dat <- cont_cov_dep_data(p = 10, n1 = 1, n2 = 1, n3 = 1)
true_graphs <- lapply(lapply(lapply(dat$true_precision, `!=`, 0), `*`, 1), `-`, diag(10))
windowsFonts("Times" = windowsFont("Times"))
true_list <- lapply(1:length(true_graphs), function(k)
  matViz(true_graphs[[k]]) +
    ggtitle(paste0("CDS ", k)) +
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)) +
    theme(plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          text = element_text(family = "Times")
          ) + theme(legend.position="none"))
true_list <- c(true_list[1], list(NULL), true_list[2], list(NULL), true_list[3])
space <- 0.5
(true_graph_z1 <- ggarrange(plotlist = true_list, nrow = 1, widths = c(1, space, 1, space, 1)))
ggsave("plots/ground_truth_CDS_z1.pdf", true_graph_z1, height = 4, width = 11)


# graph_inds <- c("$z_l\\in\\left[-3,-1\\right]$", "$z_l\\in\\left[-1,1\\right]$", "$z_l\\in\\left[-3,-1\\right]$")
# titles <- lapply(graph_inds, TeX)


# z2
dat <- cont_multi_cov_dep_data(10, 1)
true_graphs_z2 <- lapply(lapply(lapply(dat$true_precision, `!=`, 0), `*`, 1), `-`, diag(10))
ints <- apply(dat$Z, 2, cut, c(-3, -1, 1, 3))
ints <- apply(ints, 1, paste, collapse = "x")
unique_inds <- rep(NA, length(true_graphs_z2))
unique_graphs <- unique(true_graphs_z2)
unique_ints <- vector("list", length(unique_graphs))
for (i in 1:length(unique_inds)){
  unique_inds[i] <- which(sapply(unique_graphs, identical, true_graphs_z2[[i]]))
  unique_ints[[unique_inds[i]]] <- c(unique_ints[[unique_inds[i]]], ints[i])
}

true_list <- lapply(1:length(unique_graphs), function(k)
  matViz(unique_graphs[[k]], color2 = "steelblue") +
    ggtitle(paste0("CDS ", k)) +
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)) +
    theme(plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          text = element_text(family = "Times")
    ) + theme(legend.position="none"))
(true_graph_z2 <- ggarrange(plotlist = true_list, nrow = 1))
ggsave("plots/ground_truth_CDS_z2.pdf", true_graph_z2, height = 3, width = 11)

true_graph_z1 <- annotate_figure(true_graph_z1, top = text_grob(TeX("Ground Truth Conditional Dependence Structure, $\\textit{q}=1$"),
                                                                size = 18, family = "Times"))
true_graph_z2 <- annotate_figure(true_graph_z2, top = text_grob(TeX("Ground Truth Conditional Dependence Structure, $\\textit{q}=2$"),
                                                                size = 18, family = "Times"))
(true_graphs <- ggarrange(true_graph_z1, true_graph_z2, nrow = 2))
ggsave("plots/ground_truth_CDS.pdf", true_graphs, height = 6, width = 11)

true_list <- lapply(1:length(true_graphs_z2), function(k)
  matViz(true_graphs_z2[[k]]) +
    ggtitle(paste0("CDS ", k, " ", ints[k])) +
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)) +
    theme(plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),
          text = element_text(family = "Times")
    ))
for (j in 1:length(true_graphs)){
  ggsave(paste0("plots/ground_truth_CDS_z2/", ints[j], ".pdf"), true_list[[j]], height = 3, width = 11 / 3)
}



# ints <- gsub("(", "[", ints, fixed = T)
# ints <- gsub("[", "\\left[", ints, fixed = T)
# ints <- gsub("]", "\\right]", ints, fixed = T)
# unique_graphs <- unique(true_graphs)
# titles <- vector("list", length(unique_graphs))
# unique_inds <- rep(NA, length(true_graphs))
# for (i in 1:length(unique_inds)){
#   unique_inds[i] <- which(sapply(unique_graphs, identical, graphs[[i]]))
#
#   titles[[unique_inds[i]]] <- c(titles[[unique_inds[i]]], ints[i])
# }
# titles <- sapply(titles, paste0, collapse = "\\cup")
# titles <- unname(sapply(titles, function(title) substr(title, start = 5, stop = nchar(title))))
# titles <- paste0("$z_l\\in", titles, "$")
# titles <- sapply(titles, TeX)
# true_list <- lapply(1:length(unique_graphs), function(k)
#   matViz(unique_graphs[[k]]) +
#     ggtitle(titles[[k]]) +
#     coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)) +
#     theme(plot.title = element_text(size = 12),
#           axis.text = element_text(size = 18),
#           text = element_text(family = "Times")
#     ))
# (true_graph <- ggarrange(plotlist = true_list, nrow = 2, ncol = 3, common.legend = T))
# ggsave("plots/ground_truth_CDS_z2.pdf", true_graph, height = 8, width = 11)
#
#
# graphs <- lapply(true_graphs, function(g) matViz(g, color2 = "steelblue")+
#                    coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)))
# for (j in 1:length(true_graphs)){
#   ggsave(paste0("plots/ground_truth_CDS_z2/", ints[j], ".pdf"), graphs[[j]], height = 3, width = 11 / 3)
# }
