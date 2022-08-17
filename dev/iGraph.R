library(igraph)
library(covdepGE)
set.seed(1)
data <- generateData()
graphs <- unique(lapply(lapply(lapply(data$true_precision, `!=`, 0), `*`, 1), `-`, diag(5)))
graphs <- lapply(graphs, `[`, 1:3, 1:3)
g1 <- graph_from_adjacency_matrix(graphs[[1]])
g2 <- graph_from_adjacency_matrix(graphs[[2]])
g3 <- graph_from_adjacency_matrix(graphs[[3]])
svg("g1.svg")
plot(g1, edge.arrow.mode = "-", vertex.size = 70, vertex.color = "cornflowerblue",
     label.font = 1, edge.width = 2, vertex.label.cex = 7,
     vertex.label.color = "black", edge.color = "black")
dev.off()
svg("g2.svg")
plot(g2, edge.arrow.mode = "-", vertex.size = 70, vertex.color = "coral1",
     label.font = 1, edge.width = 2, vertex.label.cex = 7,
     vertex.label.color = "black", edge.color = "black")
dev.off()
svg("g3.svg")
plot(g1, edge.arrow.mode = "-", vertex.size = 70, vertex.color = "darkseagreen2",
     label.font = 1, edge.width = 2, vertex.label.cex = 7,
     vertex.label.color = "black", edge.color = "black")
dev.off()

