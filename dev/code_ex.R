library(covdepGE)
set.seed(12)
data <- generateData()

# visualize true precision matrices
lapply(c(1, 90, 121), function(l) matViz(data$true_precision[[l]], incl_val = T) +
         ggplot2::ggtitle(paste("True precision matrix, observation", l)))

# fit the model and visualize the estimated graphs
(out <- covdepGE(data$X, data$Z, parallel = T, num_workers = 5))
plot(out)

# visualize the posterior inclusion probabilities for variables (1, 3) and (1, 2)
inclusionCurve(out, 1, 2)
inclusionCurve(out, 1, 3)


library(ggplot2)
library(svglite)
library(ggpubr)
library(covdepGE)

set.seed(12)
data <- generateData()

# visualize true precision matrices
tp <- lapply(c(1, 90, 121), function(l) matViz(data$true_precision[[l]], incl_val = T) +
               ggplot2::ggtitle(paste("Observation", l)))

# fit the model and visualize the estimated graphs
(out <- covdepGE(data$X, data$Z, parallel = T, num_workers = 5))
ot_pt <- plot(out)

# visualize the posterior inclusion probabilities for variables (1, 3) and (1, 2)
inc1 <- inclusionCurve(out, 1, 2)
inc2 <- inclusionCurve(out, 1, 3)
tp

path <- "C:/Users/jacob/OneDrive/Documents/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/dev/plots"
tp_ <- ggarrange(plotlist = tp, nrow = 1, legend = F)
ot_pt_ <- ggarrange(plotlist = ot_pt, nrow = 1, legend = F)
inc <- ggarrange(plotlist = list(inc1, inc2), nrow = 1)
ggsave(paste0(path, "/true_prec.pdf"), tp_, height = 4, width = 11)
ggsave(paste0(path, "/preds.pdf"), ot_pt_, height = 4, width = 11)
ggsave(paste0(path, "/pip.pdf"), inc, height = 4, width = 11)
