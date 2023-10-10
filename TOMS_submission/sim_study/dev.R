rm(list=ls())
setwd("~/TAMU/Research/An approximate Bayesian approach to covariate dependent/covdepGE/TOMS_submission/sim_study")
library(covdepGE)
source("data.R")
set.seed(1)
data <- cont_cov_dep_sine_data(p=5, n1=75, n2=75, n3=75)
data <- cont_4_cov_dep_data(p=11, n=225)
# data <- cont_4_cov_dep_data(p=11, Z=matrix(seq(-3,3,by=0.01), length(seq(-3,3,by=0.01)), 4))

# for (i in 100:200){
#   print(i)
#   print(matViz(data$true_precision[[i]], incl_val = T) + ggtitle(i))
# }

graphs <- lapply(data$true_precision, function(mat) (mat > 0) * 1)

# find the unique graphs
unique_graphs <- unique(graphs)

# create a list where the j-th element is the j-th unique graph and the
# indices of the observations corresponding to this graph
unique_sum <- vector("list", length(unique_graphs))
names(unique_sum) <- paste0("graph", 1:length(unique_graphs))

# iterate over each of the unique graphs
for (j in 1:length(unique_graphs)){

  # fix the unique graph
  graph <- unique_graphs[[j]]

  # find indices of the observations corresponding to this graph
  graph_inds <- which(sapply(graphs, identical, graph))

  # split up the contiguous subsequences of these indices
  cont_inds <- split(sort(graph_inds), cumsum(c(1, diff(sort(graph_inds))
                                                != 1)))

  # create a character summary for each of the contiguous sequences
  inds_sum <- sapply(cont_inds, function(idx_seq) ifelse(length(
    idx_seq) > 3, paste0(min(idx_seq), ",...,", max(idx_seq)),
    paste0(idx_seq, collapse = ",")))

  # combine the summary
  inds_sum <- paste0(inds_sum, collapse = ",")

  # add the graph, indices, and summary to the unique graphs summary list
  unique_sum[[j]] <- list(graph = graph, indices = graph_inds,
                          ind_sum = inds_sum)
}

# ng <- length(unique_sum)
# i <- 0
# for (graph in unique_sum){
#   print(paste0(i,"/",ng))
#   print(matViz(graph$graph) + ggtitle(graph$ind_sum))
#   i <- i + 1
#   if (i == 100){
#     break
#   }
# }


out <- covdepGE(data$X, data$Z, parallel = T)
plot(out)


rm(list=ls())


