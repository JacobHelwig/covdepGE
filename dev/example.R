library(covdepGE)
library(ggplot2)

# get the data
set.seed(1)
data <- generate_continuous()
X <- data$data
Z <- data$covts
interval <- data$interval
prec <- data$true_precision

# get overall and within interval sample sizes
n <- nrow(X)
n1 <- sum(interval == 1)
n2 <- sum(interval == 2)

# visualize the distribution of the extraneous covariate
ggplot(data.frame(Z = Z, interval = as.factor(interval))) +
  geom_histogram(aes(Z, fill = interval), color = "black", bins = n %/% 5)

# visualize the true precision matrices in each of the intervals

# interval 1
matViz(prec[[1]], incl_val = TRUE) + ggtitle("True precision matrix, interval 1")

# interval 2 (varies continuously with Z)
int2_mats <- prec[interval == 2]
int2_inds <- c(5, n2 %/% 2, n2 - 5)
lapply(int2_inds, function(j) matViz(int2_mats[[j]], incl_val = TRUE) +
         ggtitle(paste("True precision matrix, interval 2, individual", j)))

# interval 3
matViz(prec[[length(prec)]], incl_val = TRUE) +
  ggtitle("True precision matrix, interval 3")

# fit the model and visualize the estimated precision matrices
(out <- covdepGE(X, Z))
plot(out)

# visualize the posterior inclusion probabilities for variables (1, 3) and (1, 2)
inclusionCurve(out, 1, 2)
inclusionCurve(out, 1, 3)
