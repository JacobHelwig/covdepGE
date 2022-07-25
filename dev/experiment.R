library(covdepGE)
library(doParallel)
library(doRNG)
start <- Sys.time()
(cores <- detectCores() - 5)
doParallel::registerDoParallel()
n_trials <- 100
set.seed(1)
res <- foreach(j = 1:n_trials, .packages = "covdepGE")%dorng%{
  out <- vector("list", 7)
  names(out) <- c("data", "hybrid", "hybrid500", "grid_search", "grid_search500",
                                         "model_average", "model_average500")
  out$data <- generate_continuous()
  out$hybrid <- covdepGE(out$data$data, out$data$covts, nssq = 3, nsbsq = 3,
                         npip = 3, prog_bar = F, max_iter_grid = 100)
  out$hybrid500 <- covdepGE(out$data$data, out$data$covts, nssq = 20, nsbsq = 20,
                            npip = 20, prog_bar = F, max_iter_grid = 100)
  out$grid_search <- covdepGE(out$data$data, out$data$covts, "grid_search",
                              nssq = 3, nsbsq = 3, npip = 3, prog_bar = F,
                              max_iter_grid = 100)
  out$grid_search500 <- covdepGE(out$data$data, out$data$covts, "grid_search",
                                 nssq = 20, nsbsq = 20, npip = 20, prog_bar = F,
                                 max_iter_grid = 100)
  out$model_average <- covdepGE(out$data$data, out$data$covts, "model_average",
                                nssq = 3, nsbsq = 3, npip = 3, prog_bar = F)
  out$model_average500 <- covdepGE(out$data$data, out$data$covts, "model_average",
                                   nssq = 20, nsbsq = 20, npip = 20, prog_bar = F)
  out
}
save(res, file = "res3.Rda")
Sys.time() - start

# analysis
load("dev/res2.Rda")
names(res[[1]])
res[[1]]$hybrid
res[[1]]$hybrid500
res[[1]]$model_average500

# look at marginal grids
res[[1]]$hybrid$hyperparameters
hp_names <- c("ssq", "sbsq", "pip")
marg_grids <- lapply(hp_names, function(hp) lapply(lapply(lapply(
  res[[1]]$hybrid$hyperparameters, `[[`, "grid"), `[[`, hp), "unique"))
names(marg_grids) <- hp_names
marg_grids500 <- lapply(hp_names, function(hp) lapply(lapply(lapply(
  res[[1]]$hybrid500$hyperparameters, `[[`, "grid"), `[[`, hp), "unique"))
names(marg_grids500) <- hp_names

# get the true precision
true_prec <- res[[1]]$data$true_precision

# replace diagonals with NA and get structure
for (j in 1:length(true_prec)){
  diag(true_prec[[j]]) <- NA
  true_prec[[j]] <- (true_prec[[j]] != 0) * 1
}

# get number of true 1 and 0
true1 <- sum(sapply(true_prec, function(mat) sum(mat, na.rm = T)))
true0 <- sum(sapply(true_prec, function(mat) sum(-mat + 1, na.rm = T)))

# create storage for recording performance
perf <- vector("list", length(res[[1]]) - 1)
names(perf) <- setdiff(names(res[[1]]), "data")
perf <- rep(list(perf), length(res))

j <- 1
name <- names(perf[[j]])[1]

# loop over each of the trial and hyperparameter specifications
for (j in 1:length(res)){
  for (hp_spec in names(perf[[j]])){

    # get the correct ones and zeroes
    cor1 <- sum(unlist(res[[j]][[hp_spec]]$graphs$graphs) == unlist(true_prec) &
                  unlist(true_prec) == 1, na.rm = T)
    cor0 <- sum(unlist(res[[j]][[hp_spec]]$graphs$graphs) == unlist(true_prec) &
                  unlist(true_prec) == 0, na.rm = T)

    # calculate the sensitivity and specificity and add to the performance
    # storage
    sens <- cor1 / true1
    spec <- cor0 / true0
    perf[[j]][[hp_spec]] <- list(sens = sens, spec = spec)
  }
}

# extract sensivity and specificity for each of the hyperparameter specification
# methods
library(ggplot2)
hyb_sens <- sapply(lapply(perf, `[[`, "hybrid"), `[[`, "sens")
summary(hyb_sens)
ggplot() + geom_histogram(aes(hyb_sens), color = "black")
hyb_spec <- sapply(lapply(perf, `[[`, "hybrid"), `[[`, "spec")
summary(hyb_spec)
ggplot() + geom_histogram(aes(hyb_spec), color = "black")
hyb500_sens <- sapply(lapply(perf, `[[`, "hybrid500"), `[[`, "sens")
summary(hyb500_sens)
ggplot() + geom_histogram(aes(hyb500_sens), color = "black")
hyb500_spec <- sapply(lapply(perf, `[[`, "hybrid500"), `[[`, "spec")
summary(hyb500_spec)
ggplot() + geom_histogram(aes(hyb500_spec), color = "black")

grid_sens <- sapply(lapply(perf, `[[`, "grid_search"), `[[`, "sens")
summary(grid_sens)
ggplot() + geom_histogram(aes(grid_sens), color = "black")
grid_spec <- sapply(lapply(perf, `[[`, "grid_search"), `[[`, "spec")
summary(grid_spec)
ggplot() + geom_histogram(aes(grid_spec), color = "black")
grid500_sens <- sapply(lapply(perf, `[[`, "grid_search500"), `[[`, "sens")
summary(grid500_sens)
ggplot() + geom_histogram(aes(grid500_sens), color = "black")
grid500_spec <- sapply(lapply(perf, `[[`, "grid_search500"), `[[`, "spec")
summary(grid500_spec)
ggplot() + geom_histogram(aes(grid500_spec), color = "black")

avg_sens <- sapply(lapply(perf, `[[`, "model_average"), `[[`, "sens")
summary(avg_sens)
ggplot() + geom_histogram(aes(avg_sens), color = "black")
avg_spec <- sapply(lapply(perf, `[[`, "model_average"), `[[`, "spec")
summary(avg_spec)
ggplot() + geom_histogram(aes(avg_spec), color = "black")
avg500_sens <- sapply(lapply(perf, `[[`, "model_average500"), `[[`, "sens")
summary(avg500_sens)
ggplot() + geom_histogram(aes(avg500_sens), color = "black")
avg500_spec <- sapply(lapply(perf, `[[`, "model_average500"), `[[`, "spec")
summary(avg500_spec)
ggplot() + geom_histogram(aes(avg500_spec), color = "black")

# compare

# hybrid to grid
ggplot() + geom_histogram(aes(x = grid_sens - hyb_sens), color = "black")
summary(grid_sens - hyb_sens)
ggplot() + geom_histogram(aes(x = grid_spec - hyb_spec), color = "black")
summary(grid_spec - hyb_spec)

# model average to grid
ggplot() + geom_histogram(aes(x = grid_sens - avg_sens), color = "black")
summary(grid_sens - avg_sens)
ggplot() + geom_histogram(aes(x = grid_spec - avg_spec), color = "black")
summary(grid_spec - avg_spec)

# hybrid to model average
ggplot() + geom_histogram(aes(x = hyb_sens - avg_sens), color = "black")
summary(hyb_sens - avg_sens)
ggplot() + geom_histogram(aes(x = hyb_spec - avg_spec), color = "black")
summary(hyb_spec - avg_spec)

# hybrid to hybrid500
ggplot() + geom_histogram(aes(x = hyb500_sens - hyb_sens), color = "black")
summary(hyb500_sens - hyb_sens)
ggplot() + geom_histogram(aes(x = hyb500_spec - hyb_spec), color = "black")
summary(hyb500_spec - hyb_spec)

# grid to grid500
ggplot() + geom_histogram(aes(x = grid500_sens - grid_sens), color = "black")
summary(grid500_sens - grid_sens)
ggplot() + geom_histogram(aes(x = grid500_spec - grid_spec), color = "black")
summary(grid500_spec - grid_spec)

# avg to avg500
ggplot() + geom_histogram(aes(x = avg500_sens - avg_sens), color = "black")
summary(avg500_sens - avg_sens)
ggplot() + geom_histogram(aes(x = avg500_spec - avg_spec), color = "black")
summary(avg500_spec - avg_spec)

# hyb500 to grid500
ggplot() + geom_histogram(aes(x = grid500_sens - hyb500_sens), color = "black")
summary(grid500_sens - hyb500_sens)
ggplot() + geom_histogram(aes(x = grid500_spec - hyb500_spec), color = "black")
summary(grid500_spec - hyb500_spec)

# avg500 to grid500
ggplot() + geom_histogram(aes(x = grid500_sens - avg500_sens), color = "black")
summary(grid500_sens - avg500_sens)
ggplot() + geom_histogram(aes(x = grid500_spec - avg500_spec), color = "black")
summary(grid500_spec - avg500_spec)

# hyb500 to avg500
ggplot() + geom_histogram(aes(x = hyb500_sens - avg500_sens), color = "black")
summary(hyb500_sens - avg500_sens)
ggplot() + geom_histogram(aes(x = hyb500_spec - avg500_spec), color = "black")
summary(hyb500_spec - avg500_spec)
