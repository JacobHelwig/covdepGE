# This function generates a square matrix of uniformly distributed random
# numbers, finds the corresponding (complex) eigenvalues and then selects the
# eigenvalue with the largest modulus. The dimensions of the matrix and the
# standard deviation of the random numbers are given as input parameters.
max.eig <- function(N, sigma) {
  d <- matrix(rnorm(N**2, sd = sigma), nrow = N)

  E <- eigen(d)$values

  abs(E)[[1]]
}

# The foreach() functionality can be applied to a cluster using the doSNOW
# library. We will start by using doSNOW to create a collection of R instances
# on a single machine using a SOCK cluster.
library(rbenchmark)
library(doSNOW)
getDoParRegistered()
cluster <- makeCluster(4, type = "SOCK")
registerDoSNOW(cluster)

benchmark(
  foreach(n = 1:50) %do% max.eig(n, 1),
  foreach(n = 1:50) %dopar% max.eig(n, 1)
)

stopCluster(cluster)




par_sleep <- function(stop_cluster = T){

  if (!foreach::getDoParRegistered() | nrow(showConnections()) == 0){
    warning("Active parallel backend not detected; registering parallel backend using 3 cores")
    doParallel::registerDoParallel(cores = 3)
  } else{
    message(paste("Detected", foreach::getDoParWorkers(), "registered workers on an active cluster"))
  }
  start <- Sys.time()
  foreach::`%dopar%`(foreach::foreach(j = 1:3), Sys.sleep(1))
  print(Sys.time() - start)
  if(stop_cluster) doParallel::stopImplicitCluster()
}

par_sleep()
par_sleep(F)

cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)
par_sleep()
parallel::stopCluster(cl)
par_sleep()


par_sleep()
par_sleep()


