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
