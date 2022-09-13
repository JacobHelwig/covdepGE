library(covdepGE)
library(testthat)
set.seed(1)
data <- generateData()

test_that("Runtime is reasonable", {
  skip_on_cran()
  out <- covdepGE(data$X, data$Z, prog_bar = F)
  print(out$model_details$elapsed)
  expect_lt(as.numeric(out$model_details$elapsed, units = "mins"), 5)
})

test_that("Wrong size X and Z", {
  expect_error(covdepGE(data$X, data$Z[-1], prog_bar = F))
})

test_that("Constant Z gives 2 warnings", {
  expect_equal(2,length(capture_warnings(covdepGE(
    data$X, ssq = 0.5, sbsq = 0.5, pip = 0.1, prog_bar = F))))
})

test_that("Constant Z gives 1 graph", {
  out <- covdepGE(data$X, ssq = 0.5, sbsq = 0.5, pip = 0.1, tau = 0.5,
                  scale_Z = F, prog_bar = F)
  expect_equal(out$model_details$num_unique, 1)
})

test_that("With 1 HP candidate, all hp_methods give the same results", {
  out_hyb <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1,
                      max_iter = 0, max_iter_grid = 100, prog_bar = F)
  out_grid <- covdepGE(data$X, data$Z, hp_method = "grid_search",
                       ssq = 0.5, sbsq = 0.5, pip = 0.1,
                       max_iter = 0, max_iter_grid = 100, prog_bar = F)
  out_avg <- covdepGE(data$X, data$Z, hp_method = "model_average",
                       ssq = 0.5, sbsq = 0.5, pip = 0.1, prog_bar = F)
  expect_equal(out_hyb$variational_params, out_grid$variational_params)
  expect_equal(out_grid$variational_params, out_avg$variational_params)
})

test_that("With 8 HP candidates, all hp_methods give different results", {
  out_hyb <- covdepGE(data$X, data$Z, nssq = 2, nsbsq = 2, npip = 2,
                      max_iter = 0, max_iter_grid = 100, prog_bar = F)
  out_grid <- covdepGE(data$X, data$Z, hp_method = "grid_search",
                       nssq = 2, nsbsq = 2, npip = 2,
                       max_iter = 0, max_iter_grid = 100, prog_bar = F)
  out_avg <- covdepGE(data$X, data$Z, hp_method = "model_average",
                      nssq = 2, nsbsq = 2, npip = 2, prog_bar = F)
  expect_false(isTRUE(all.equal(
    out_hyb$variational_params, out_grid$variational_params)))
  expect_false(isTRUE(all.equal(
    out_hyb$variational_params, out_avg$variational_params)))
  expect_false(isTRUE(all.equal(
    out_grid$variational_params, out_avg$variational_params)))
})

test_that("HP grid may be specified directly", {
  ssq <- c(0.1, 0.5)
  sbsq <- c(0.1, 0.5)
  pip <- c(0.01, 0.1)
  out <- covdepGE(data$X, data$Z, ssq = ssq, sbsq = sbsq, pip = pip,
                  prog_bar = F)
  expect_equal(0,
               sum(
                 out$hyperparameters$variable1$grid[, c("pip", "ssq", "sbsq")] -
                   expand.grid(pip, ssq, sbsq)))
})

test_that("HP grid may be specified by count", {
  nssq <- 2
  nsbsq <- 3
  npip <- 4
  out <- covdepGE(data$X, data$Z, nssq = nssq, nsbsq = nsbsq, npip = npip,
                  prog_bar = F)
  expect_equal(out$model_details$grid_size,
               prod(nssq, nsbsq, npip))
})

test_that("Small pip results in sparser estimates", {
  out_small <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.01,
                        prog_bar = F)
  out_large <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1,
                        prog_bar = F)
  expect_lt(sum(unlist(out_small$graphs$graphs)),
            sum(unlist(out_large$graphs$graphs)))
})

test_that("Larger ssq_mult gives larger ssq candidates", {
  out_small <- covdepGE(data$X, data$Z, nssq = 2, sbsq = 0.5, pip = 0.1,
                        prog_bar = F)
  out_large <- covdepGE(data$X, data$Z, ssq_mult = 2, nssq = 2, sbsq = 0.5,
                        pip = 0.1, prog_bar = F)
  expect_true(all(
    sort(unique(out_small$hyperparameters$variable1$grid$ssq)) <=
      sort(unique(out_large$hyperparameters$variable1$grid$ssq))))
})

test_that("Larger snr_upper gives larger sbsq candidates", {
  out_small <- covdepGE(data$X, data$Z, ssq = 0.5, nsbsq = 2, pip = 0.1,
                        prog_bar = F)
  out_large <- covdepGE(data$X, data$Z, ssq = 0.5, nsbsq = 2, pip = 0.1,
                        snr_upper = 35, prog_bar = F)
  expect_true(all(
    sort(unique(out_small$hyperparameters$variable1$grid$sbsq)) <=
      sort(unique(out_large$hyperparameters$variable1$grid$sbsq))))
})

test_that("Larger pip_upper gives larger pip candidates", {
  out_small <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, npip = 2,
                        prog_bar = F)
  out_large <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, npip = 2,
                        pip_upper = 0.9, prog_bar = F)
  expect_true(all(
    sort(unique(out_small$hyperparameters$variable1$grid$pip)) <=
      sort(unique(out_large$hyperparameters$variable1$grid$pip))))
})

test_that("The lower bound controls the lower endpoint of the HP candidates", {
  lower <- 0.1
  out <- covdepGE(data$X, data$Z, nssq = 2, nsbsq = 2, npip = 2,
                  ssq_lower = lower, sbsq_lower = lower, pip_lower = lower,
                  prog_bar = F)
  expect_true(all(
    rep(lower, 3) == apply(
      out$hyperparameters$variable1$grid[ , c("pip", "ssq", "sbsq")], 2, min)))
})

test_that("Large tau produces 1 graph", {
  out <- covdepGE(data$X, data$Z, tau = 1e5, ssq = 0.5, sbsq = 0.5, pip = 0.1,
                  prog_bar = F)
  expect_equal(1, out$model_details$num_unique)
})

test_that("Choice of norm does not change weights when Z is 1D", {
  out_inf <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1,
                      norm = Inf, prog_bar = F)
  out_2 <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1,
                    prog_bar = F)
  expect_equal(out_inf$weights$weights, out_2$weights$weights)
})

test_that("Choice of norm changes weights when Z is 2D", {
  Z2 <- cbind(data$Z, rnorm(length(data$Z)))
  out_3 <- covdepGE(data$X, Z2, ssq = 0.5, sbsq = 0.5, pip = 0.1, norm = 3,
                    prog_bar = F)
  out_2 <- covdepGE(data$X, Z2, ssq = 0.5, sbsq = 0.5, pip = 0.1, prog_bar = F)
  expect_false(isTRUE(all.equal(out_3$weights$weights, out_2$weights$weights)))
})

test_that("center_X is working", {
  out_cent <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, npip = 0.1,
                       prog_bar = F)
  out_nocent <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, npip = 0.1,
                       center_X = F, prog_bar = F)
  out_mancent <- covdepGE(scale(data$X, scale = F), data$Z, ssq = 0.5, sbsq = 0.5, npip = 0.1,
                         center_X = F, prog_bar = F)
  expect_false(isTRUE(all.equal(out_cent$variational_params,
                                out_nocent$variational_params)))
  expect_equal(out_cent$variational_params, out_mancent$variational_params)
})

test_that("scale_Z is working", {
  out_scale <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, npip = 0.1,
                        prog_bar = F)
  out_noscale <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, npip = 0.1,
                         scale_Z = F, prog_bar = F)
  out_manscale <- covdepGE(data$X, scale(data$Z), ssq = 0.5, sbsq = 0.5, npip = 0.1,
                            scale_Z = F, prog_bar = F)
  expect_false(isTRUE(all.equal(out_scale$variational_params,
                                out_noscale$variational_params)))
  expect_equal(out_scale$variational_params, out_manscale$variational_params)
})

test_that("Greater alpha_tol gives faster convergence", {
  skip_on_cran()
  out_slow <- covdepGE(data$X, data$Z, alpha_tol = 1e-12, prog_bar = F)
  out_fast <- covdepGE(data$X, data$Z, alpha_tol = 1e-1, prog_bar = F)
  expect_gt(out_slow$model_details$elapsed, out_fast$model_details$elapsed)
})

test_that("Greater max_iter and max_iter_grid gives slower convergence", {
  skip_on_cran()
  out_slow1 <- covdepGE(data$X, data$Z, nssq = 2, nsbsq = 2, npip = 2,
                        max_iter_grid = 1000, prog_bar = F)
  out_slow2 <- covdepGE(data$X, data$Z, nssq = 2, nsbsq = 2, npip = 2,
                        max_iter = 1000, prog_bar = F)
  out_fast <- covdepGE(data$X, data$Z, nssq = 2, nsbsq = 2, npip = 2,
                       prog_bar = F)
  expect_gt(out_slow1$model_details$elapsed, out_fast$model_details$elapsed)
  expect_gt(out_slow2$model_details$elapsed, out_fast$model_details$elapsed)
})

test_that("Greater edge_threshold gives more sparse estimates", {
  out_sparse <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1,
                         edge_threshold = 0.75, prog_bar = F)
  out <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1, prog_bar = F)
  expect_lt(sum(unlist(out_sparse$graphs$graphs)),
            sum(unlist(out$graphs$graphs)))
})

test_that("sym_method affects sparsity", {
  out_sparse <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1,
                         sym_method = "min", prog_bar = F)
  out <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1, prog_bar = F)
  out_full <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1,
                       sym_method = "max", prog_bar = F)
  expect_lt(sum(unlist(out_sparse$graphs$graphs)),
            sum(unlist(out$graphs$graphs)))
  expect_lt(sum(unlist(out$graphs$graphs)),
            sum(unlist(out_full$graphs$graphs)))
})

test_that("parallel gives same results as sequential", {
  out_par1 <- suppressMessages(covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5,
                                        pip = 0.1, parallel = T,
                                        num_workers = 1))
  doParallel::registerDoParallel(1)
  out_par2 <- suppressMessages(covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5,
                                       pip = 0.1, parallel = T))
  out_seq <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1,
                      prog_bar = F)
  expect_equal(out_par1$variational_params, out_seq$variational_params)
  expect_equal(out_par2$variational_params, out_seq$variational_params)
})

test_that("The progress bar can be turned off", {
  expect_false(0 == nchar(capture_output(covdepGE(
    data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1))))
  expect_equal(0, nchar(capture_output(covdepGE(
    data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1, prog_bar = F))))
})

test_that("print and summary give the same results", {
  out <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1, prog_bar = F)
  expect_equal(capture_output(print(out)),
               capture_output(summary(out)))
})

# generateData

test_that("p controls the size of the data", {
  p <- 9
  data <- generateData(p = p)
  expect_equal(p, ncol(data$X))
})

test_that("n1, n2, n3 control the sample size", {
  n1 <- 1
  n2 <- 2
  n3 <- 3
  data <- generateData(n1 = n1, n2 = n2, n3 = n3)
  expect_equal(data$interval, c(rep(1, n1), rep(2, n2), rep(3, n3)))
})


test_that("Z controls interval", {
  int1 <- c(-3, -1)
  int2 <- c(-1, 1)
  int3 <- c(1, 3)
  n1 <- 1
  n2 <- 2
  n3 <- 3
  Z <- c(runif(n1, int1[1], int1[2]), runif(n2, int2[1], int2[2]),
         runif(n3, int3[1], int3[2]))
  data <- generateData(Z = Z) # !NULL Z, NULL true_precision
  expect_equal(data$interval, c(rep(1, n1), rep(2, n2), rep(3, n3)))
  expect_equal(Z, data$Z)
})

test_that("Passing both Z and true_precision gives an error", {
  expect_error(generateData(Z = data$Z, true_precision = data$true_precision))
})

test_that("Passing true_precision leaves Z null", {
  expect_null(generateData(true_precision = data$true_precision)$Z)
})

test_that("Z controls precision matrices, precision matrices from Z control data", {

  # interval 1
  int1 <- c(-3, -1)
  n1 <- 100
  Z <- c(runif(n1, int1[1], int1[2]))
  data_2 <- generateData(Z = Z)
  expect_equal(1, length(unique(data_2$true_precision)))
  expect_equal(data$true_precision[[1]], data_2$true_precision[[1]])
  out <- covdepGE(data_2$X, data_2$Z, tau = 1e5, prog_bar = F)
  expect_equal(1, out$model_details$num_unique)
  graph <- (data_2$true_precision[[1]] != 0) * 1 - diag(ncol(data_2$X))
  expect_equal(graph, out$graphs$unique_graphs$graph1$graph)

  # interval 3
  int1 <- c(1, 3)
  Z <- c(runif(n1, int1[1], int1[2]))
  data_2 <- generateData(Z = Z)
  expect_equal(1, length(unique(data_2$true_precision)))
  expect_equal(data$true_precision[[nrow(data$X)]], data_2$true_precision[[1]])
  out <- covdepGE(data_2$X, data_2$Z, tau = 1e5, prog_bar = F)
  expect_equal(1, out$model_details$num_unique)
  graph <- (data_2$true_precision[[1]] != 0) * 1 - diag(ncol(data_2$X))
  expect_equal(graph, out$graphs$unique_graphs$graph1$graph)
})

test_that("Precision matrices control data", {

  # interval 1
  n1 <- 100
  prec <- rep(list(data$true_precision[[1]]), n1)
  data_2 <- generateData(true_precision = prec)
  expect_equal(1, length(unique(data_2$true_precision)))
  expect_equal(data$true_precision[[1]], data_2$true_precision[[1]])
  out <- suppressWarnings(covdepGE(data_2$X, prog_bar = F))
  expect_equal(1, out$model_details$num_unique)
  graph <- (data_2$true_precision[[1]] != 0) * 1 - diag(ncol(data_2$X))
  expect_equal(graph, out$graphs$unique_graphs$graph1$graph)

  # interval 3
  prec <- rep(list(data$true_precision[[nrow(data$X)]]), n1)
  data_2 <- generateData(true_precision = prec)
  expect_equal(1, length(unique(data_2$true_precision)))
  expect_equal(data$true_precision[[nrow(data$X)]], data_2$true_precision[[1]])
  out <- suppressWarnings(covdepGE(data_2$X, prog_bar = F))
  expect_equal(1, out$model_details$num_unique)
  graph <- (data_2$true_precision[[1]] != 0) * 1 - diag(ncol(data_2$X))
  expect_equal(graph, out$graphs$unique_graphs$graph1$graph)
})

# inclusionCurve

out <- covdepGE(data$X, data$Z, ssq = 0.5, sbsq = 0.5, pip = 0.1, prog_bar = F)
inc_curve <- inclusionCurve(out, 1, 2)

test_that("Error is thrown when out is not of class covdepGE", {
  expect_error(inclusionCurve(NULL, 1, 2))
})

test_that("Save the default inclusion curve plot", {
  suppressWarnings(vdiffr::expect_doppelganger("inc_curve", inc_curve))
  expect_true(T)
})

test_that("Verify that vdiffr::expect_doppelganger is working", {
  vdiffr::expect_doppelganger("inc_curve", inclusionCurve(out, 1, 2))
})

test_that("Different vertices give different plots", {
  expect_failure(vdiffr::expect_doppelganger(
    "inc_curve", inclusionCurve(out, 1, 3)))
})

test_that("Different vertices give different plots", {
  expect_failure(vdiffr::expect_doppelganger(
    "inc_curve", inclusionCurve(out, 1, 3)))
})

test_that("Different line types", {
  expect_failure(vdiffr::expect_doppelganger(
    "inc_curve", inclusionCurve(out, 1, 2, line_type = "dotted")))
})

test_that("Different line size", {
  expect_failure(vdiffr::expect_doppelganger(
    "inc_curve", inclusionCurve(out, 1, 2, line_size = 1)))
})

test_that("Different line color", {
  expect_failure(vdiffr::expect_doppelganger(
    "inc_curve", inclusionCurve(out, 1, 2, line_color = "dodgerblue")))
})

test_that("Different line color", {
  expect_failure(vdiffr::expect_doppelganger(
    "inc_curve", inclusionCurve(out, 1, 2, line_color = "dodgerblue")))
})

test_that("Different point shape", {
  expect_failure(vdiffr::expect_doppelganger(
    "inc_curve", inclusionCurve(out, 1, 2, point_shape = 2)))
})

test_that("Different point size", {
  expect_failure(vdiffr::expect_doppelganger(
    "inc_curve", inclusionCurve(out, 1, 2, point_size = 2)))
})

test_that("Different point color", {
  expect_failure(vdiffr::expect_doppelganger(
    "inc_curve", inclusionCurve(out, 1, 2, point_color = "dodgerblue")))
})

test_that("Different point fill", {
  expect_failure(vdiffr::expect_doppelganger(
    "inc_curve", inclusionCurve(out, 1, 2, point_fill = "grey55")))
})

# matViz

adj_mat <- out$graphs$graphs[[1]]
a_viz <- matViz(adj_mat)
prec_mat <- out$graphs$inclusion_probs_sym[[1]]
p_viz <- matViz(prec_mat)
rcname_mat <- prec_mat
rownames(rcname_mat) <- colnames(rcname_mat) <- paste0("V", 1:ncol(rcname_mat))

test_that("Error is thrown when x is not of class matrix", {
  expect_error(matViz(1))
})

test_that("Save the default matViz plots", {
  suppressWarnings(vdiffr::expect_doppelganger("a", a_viz))
  suppressWarnings(vdiffr::expect_doppelganger("p", p_viz))
  expect_true(T)
})

test_that("Verify that vdiffr::expect_doppelganger is working", {
  vdiffr::expect_doppelganger("a", matViz(adj_mat))
  vdiffr::expect_doppelganger("p", matViz(prec_mat))
})

test_that("Row and column names are added back", {
  expect_failure(
    vdiffr::expect_doppelganger("p", matViz(rcname_mat)))
})

test_that("Different color 1", {
  expect_failure(
    vdiffr::expect_doppelganger("a", matViz(adj_mat, color1 = "grey90")))
  expect_failure(
    vdiffr::expect_doppelganger("p", matViz(prec_mat, color1 = "grey90")))
})

test_that("Different color 2", {
  expect_failure(
    vdiffr::expect_doppelganger("a", matViz(adj_mat, color2 = "dodgerblue")))
  expect_failure(
    vdiffr::expect_doppelganger("p", matViz(prec_mat, color2 = "dodgerblue")))
})

test_that("Different grid color", {
  expect_failure(
    vdiffr::expect_doppelganger("a", matViz(adj_mat, grid_color = "dodgerblue")))
  expect_failure(
    vdiffr::expect_doppelganger("p", matViz(prec_mat, grid_color = "dodgerblue")))
})

test_that("Include cell values", {
  vdiffr::expect_doppelganger(
    "a", matViz(adj_mat, incl_val = T, prec = 5, font_size = 7,
                font_color1 = "forestgreen", font_color2 = "grey90",
                font_thres = 0.2))
  expect_failure(
    vdiffr::expect_doppelganger("p", matViz(prec_mat, incl_val = T)))
})

test_that("Precision", {
  expect_failure(
    vdiffr::expect_doppelganger("p", matViz(prec_mat, incl_val = T,
                                            prec = 4)))
})

test_that("Font size", {
  expect_failure(
    vdiffr::expect_doppelganger("p", matViz(prec_mat, incl_val = T,
                                            font_size = 4)))
})

test_that("Font color 1", {
  expect_failure(
    vdiffr::expect_doppelganger("p", matViz(prec_mat, incl_val = T,
                                            font_color1 = "dodgerblue")))
})

test_that("Font color 2", {
  expect_failure(
    vdiffr::expect_doppelganger("p", matViz(prec_mat, incl_val = T,
                                            font_color2 = "dodgerblue")))
})

test_that("Font threshold", {
  expect_failure(
    vdiffr::expect_doppelganger("p", matViz(prec_mat, incl_val = T,
                                            font_thres = 0.05)))
})

# plot.covdepGE

plots <- plot(out)
plots2 <- plot(out, c("dodgerblue", "forestgreen"))

test_that("Save the default plot lists", {
  suppressWarnings(vdiffr::expect_doppelganger("plots_1", plots[[1]]))
  suppressWarnings(vdiffr::expect_doppelganger("plots_2", plots[[2]]))
  suppressWarnings(vdiffr::expect_doppelganger("plots2_1", plots2[[1]]))
  suppressWarnings(vdiffr::expect_doppelganger("plots2_2", plots2[[2]]))
  expect_true(T)
})

test_that("Verify that vdiffr::expect_doppelganger is working and that colors are recycled", {
  vdiffr::expect_doppelganger("plots_1", plot(out)[[1]])
  vdiffr::expect_doppelganger(
    "plots2_1", plot(out, c("dodgerblue", "forestgreen"))[[1]])
  vdiffr::expect_doppelganger(
    "plots2_2", plot(out, c("dodgerblue", "forestgreen"))[[2]])
})

test_that("Different graph_colors", {
  expect_failure(vdiffr::expect_doppelganger(
    "plots_2", plot(out, graph_colors = "dodgerblue")[[2]]))
})

test_that("No title summary", {
  expect_failure(vdiffr::expect_doppelganger(
    "plots_1", plot(out, title_sum = F)[[1]]))
})
