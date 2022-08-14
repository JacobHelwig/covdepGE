data <- generateData()

test_that("Wrong size X and Z", {
  expect_error(covdepGE(data$X, data$Z[-1]))
})

test_that("With 1 HP candidate, all hp_methods give the same results", {
  out_hyb <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1,
                      max_iter = 0, max_iter_grid = 100)
  out_grid <- covdepGE(data$X, data$Z, hp_method = "grid_search",
                       ssq = 0.5, sbsq = 0.5, pip = 0.1,
                       max_iter = 0, max_iter_grid = 100)
  out_avg <- covdepGE(data$X, data$Z, hp_method = "model_average",
                       ssq = 0.5, sbsq = 0.5, pip = 0.1)
  expect_equal(out_hyb$variational_params$alpha, out_grid$variational_params$alpha)
  expect_equal(out_grid$variational_params$alpha, out_avg$variational_params$alpha)
})

test_that("With 8 HP candidates, all hp_methods give different results", {
  out_hyb <- covdepGE(data$X, data$Z, nssq = 2, nsbsq = 2, npip = 2,
                      max_iter = 0, max_iter_grid = 100)
  out_grid <- covdepGE(data$X, data$Z, hp_method = "grid_search",
                       nssq = 2, nsbsq = 2, npip = 2,
                       max_iter = 0, max_iter_grid = 100)
  out_avg <- covdepGE(data$X, data$Z, hp_method = "model_average",
                      nssq = 2, nsbsq = 2, npip = 2)
  expect_false(isTRUE(all.equal(
    out_hyb$variational_params$alpha, out_grid$variational_params$alpha)))
  expect_false(isTRUE(all.equal(
    out_hyb$variational_params$alpha, out_avg$variational_params$alpha)))
  expect_false(isTRUE(all.equal(
    out_grid$variational_params$alpha, out_avg$variational_params$alpha)))
})

test_that("HP grid may be specified directly", {
  ssq <- c(0.1, 0.5)
  sbsq <- c(0.1, 0.5)
  pip <- c(0.01, 0.1)
  out <- covdepGE(data$X, data$Z, ssq = ssq, sbsq = sbsq, pip = pip)
  expect_equal(0,
               sum(
                 out$hyperparameters$variable1$grid[, c("pip", "ssq", "sbsq")] -
                   expand.grid(pip, ssq, sbsq)))
})

test_that("HP grid may be specified by count", {
  nssq <- 2
  nsbsq <- 3
  npip <- 4
  out <- covdepGE(data$X, data$Z, nssq = nssq, nsbsq = nsbsq, npip = npip)
  expect_equal(out$model_details$grid_size,
               prod(nssq, nsbsq, npip))
})

test_that("Small pip results in sparser estimates", {
  out_small <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.01)
  out_large <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1)
  expect_lt(sum(unlist(out_small$graphs$graphs)),
            sum(unlist(out_large$graphs$graphs)))
})


test_that("Large tau produces 1 graph", {
  out <- covdepGE(data$X, data$Z, tau = 1e5, ssq = 0.5, sbsq = 0.5, pip = 0.1)
  expect_equal(1, out$model_details$num_unique)
})
