library(igraph)
library(covdepGE)
set.seed(1)
data <- generateData(5, 1, 1, 1)
graphs <- lapply(lapply(lapply(data$true_precision, `!=`, 0), `*`, 1), `-`, diag(5))
#graphs <- lapply(graphs, `[`, 1:3, 1:3)
g1 <- graph_from_adjacency_matrix(graphs[[1]])
g2 <- graph_from_adjacency_matrix(graphs[[2]])
g3 <- graph_from_adjacency_matrix(graphs[[3]])
#svg("g1.svg")
pdf("g1.pdf")
plot(g1, edge.arrow.mode = "-", vertex.size = 50, vertex.color = "cornflowerblue",
     edge.width = 2, vertex.label.cex = 5,
     vertex.label.color = "black", edge.color = "black")
title(main = "Observations 1,...,50", cex.main = 3)
dev.off()
#svg("g2.svg")
pdf("g2.pdf")
plot(g2, edge.arrow.mode = "-", vertex.size = 50, vertex.color = "coral1",
     edge.width = 2, vertex.label.cex = 5,
     vertex.label.color = "black", edge.color = "black")
title(main = "Observations 51,...,100", cex.main = 3)
dev.off()
svg("g3.svg")
plot(g1, edge.arrow.mode = "-", vertex.size = 70, vertex.color = "darkseagreen2",
     label.font = 1, edge.width = 2, vertex.label.cex = 7,
     vertex.label.color = "black", edge.color = "black")
dev.off()

