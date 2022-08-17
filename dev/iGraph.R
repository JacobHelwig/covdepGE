library(igraph)
library(covdepGE)
set.seed(1)
data <- generateData()
graphs <- unique(lapply(lapply(lapply(data$true_precision, `!=`, 0), `*`, 1), `-`, diag(5)))
g1 <- graph_from_adjacency_matrix(graphs[[1]])
g2 <- graph_from_adjacency_matrix(graphs[[2]])
g3 <- graph_from_adjacency_matrix(graphs[[3]])
g1p <- plot(g1, edge.arrow.mode = "-", vertex.size = 30,
            vertex.color = "cornflowerblue", edge.width = 2, vertex.label.cex = 2)
g1p <- plot(g1, edge.arrow.mode = "-", vertex.size = 30,
            vertex.color = "cornflowerblue", edge.width = 2, vertex.label.cex = 2)

