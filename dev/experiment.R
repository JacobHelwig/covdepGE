library(covdepGE)
library(doParallel)
library(doRNG)
cores <- detectCores() - 5
doParallel::registerDoParallel()
n_trials <- 1
set.seed(1)
foreach(j = 1:n_trials)%dorng%{
  res <- vector("list", 7)
  names(res) <- c("data", "hybrid", "hybrid500", "grid_search", "grid_search500",
                  "model_average", "model_average500")
  res[["data"]] <- generate_continuous()
  res[["hybrid"]] <- covdepGE(data$data, data$covts)
  res[["hybrid500"]] <- covdepGE(data$data, data$covts,
                                 nssq = 1, nsbsq = 1, npip = 1)
  res[["grid_search"]] <- covdepGE(data$data, data$covts, "grid_search")
  res[["grid_search500"]] <- covdepGE(data$data, data$covts, "grid_search",
                                      nssq = 1, nsbsq = 1, npip = 1)
  res[["model_average"]] <- covdepGE(data$data, data$covts, "model_average")
  res[["model_average500"]] <- covdepGE(data$data, data$covts, "model_average",
                                        nssq = 1, nsbsq = 1, npip = 1)
}
save(res, file = "res.Rda")
