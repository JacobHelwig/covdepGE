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
