## _____________________________________________________________________________
## _____________________________covdepGE_checks_________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## function to check compatiibility of arguments of covdepGE
## -----------------------------ARGUMENTS---------------------------------------
## data_mat: n x (p + 1) matrix; data
##
## Z: n x p' matrix; extraneous covariates
##
## tau: scalar in (0, Inf) OR n x 1 vector, entries in (0, Inf)
##
## kde: logical; if T, use 2-step KDE methodology
##
## alpha: scalar in [0, 1]; global initialization value
##
## mu: scalar; global initialization value
##
## sigmasq: scalar in (0, Inf); variance hyperparameter for spike-and-slab.
##
## sigmabetasq_vec: n_sigma x 1 vector, entries in (0, Inf); candidate values
##
## var_min: scalar in (0, Inf); lower bound of sigmabetasq_vec
##
## var_max: scalar in (0, Inf); upper bound of sigmabetasq_vec
##
## n_sigma: scalar in {1, 2,...}
##
## pi_vec: n_pi x 1 vector, entries in [0, 1]; candidate values of pi
##
## norm: scalar in [1, Inf]; norm to use when calculating weights
##
## scale: logical; if T, center and scale extraneous covariates
##
## tolerance: scalar in (0, Inf); end the variational update loop
##
## max_iter: scalar in {1, 2,...}; end the variational update loop
##
## edge_threshold: scalar in (0, 1)
##
## sym_method: character in {"mean", "max", "min"}
##
## monitor_final_elbo logical scalar; if `T`, monitor the final model ELBO
##
## monitor_cand_elbo logical scalar; if `T`, monitor the non-convergent ELBO history
##
## monitor_period scalar in \eqn{{1, 2,..., max_iter}}; the periodicity of ELBO monitoring
##
## warnings: logical; if T, convergence and grid warnings will be displayed
##
covdepGE_checks <- function(data_mat, Z, tau, kde, alpha, mu, sigmasq,
                            sigmabetasq_vec, var_min, var_max, n_sigma, pi_vec,
                            norm, scale, tolerance, max_iter, edge_threshold,
                            sym_method, monitor_final_elbo, monitor_cand_elbo,
                            monitor_period, print_time, warnings){

  # ensure vector input for parameters that are expected to be vectors
  args_vector <- list(tau = tau, pi_vec = pi_vec)
  if (any(!sapply(args_vector, is.vector))){

    # get the name of the non-vector
    non_vector <- names(which(!sapply(args_vector, is.vector)))[1]
    stop(paste0(non_vector, " is of class ", class(args_vector[[non_vector]])[1],
                "; expected vector"))
  }

  # ensure scalar input for parameters that are expected to be scalars
  args_scalar <- list(kde = kde, alpha = alpha, mu = mu, sigmasq = sigmasq,
                      var_min = var_min, var_max = var_max, n_sigma = n_sigma,
                      norm = norm, scale = scale, tolerance = tolerance,
                      max_iter = max_iter, edge_threshold = edge_threshold,
                      sym_method = sym_method,
                      monitor_final_elbo = monitor_final_elbo,
                      monitor_cand_elbo = monitor_cand_elbo,
                      monitor_period = monitor_period, print_time = print_time,
                      warnings = warnings)
  if (any(!(sapply(args_scalar, length) == 1))){

    # get the name of the non-scalar
    non_scalar <- names(which(!(sapply(args_scalar, length) == 1)))[1]
    stop(paste0(non_scalar, " is of length ",
                length(args_scalar[[non_scalar]]), "; expected scalar value"))
  }

  # ensure numeric input for parameters that are expected to be numeric
  args_numeric <- list(data_mat = data_mat, Z = Z, tau = tau, alpha = alpha,
                       mu = mu, sigmasq = sigmasq, var_min = var_min,
                       n_sigma = n_sigma, var_max = var_max, pi_vec = pi_vec,
                       norm = norm, tolerance = tolerance, max_iter = max_iter,
                       edge_threshold = edge_threshold,
                       monitor_period = monitor_period)
  if (any(!sapply(args_numeric, is.numeric))){

    # get the name of the non-numeric
    non_numeric <- names(which(!sapply(args_numeric, is.numeric)))[1]
    stop(paste0(non_numeric, " is of type ",
                typeof(args_numeric[[non_numeric]]), "; expected numeric"))
  }

  # ensure non-NA input for parameters that are to be non-NA
  args_nonNA <- list(data_mat = data_mat, Z = Z, tau = tau, kde = kde,
                     alpha = alpha, mu = mu, sigmasq = sigmasq,
                     var_min = var_min, var_max = var_max, n_sigma= n_sigma,
                     pi_vec = pi_vec, norm = norm, scale = scale,
                     tolerance = tolerance, max_iter = max_iter,
                     edge_threshold = edge_threshold, sym_method = sym_method,
                     monitor_final_elbo = monitor_final_elbo,
                     monitor_cand_elbo = monitor_cand_elbo,
                     monitor_period = monitor_period, print_time = print_time,
                     warnings = warnings)
  if (any(sapply(args_nonNA, function (x) any(is.na(x))))){

    # get the name of the NA
    na <- names(which(sapply(args_nonNA, function (x) any(is.na(x)))))[1]
    stop(paste0(na, " should have all non-NA entries"))
  }

  # ensure finite input for parameters that are expected to be finite
  args_finite <- list(data_mat = data_mat, Z = Z, tau = tau, mu = mu,
                      sigmasq = sigmasq, var_min = var_min, var_max = var_max,
                      n_sigma = n_sigma, tolerance = tolerance,
                      max_iter = max_iter, monitor_period = monitor_period)
  if (any(!sapply(args_finite, function (x) all(is.finite(x))))){

    # get the name of the non-finite
    non_finite <- names(which(!sapply(args_finite,
                                      function (x) all(is.finite(x)))))[1]
    stop(paste0(non_finite, " should have all finite entries"))
  }

  # ensure integer input for parameters that are expected to be integers
  args_integers <- list(n_sigma = n_sigma, max_iter = max_iter,
                        monitor_period = monitor_period)
  if (any(!(sapply(args_integers, function(x) x %% 1) == 0))){

    # get the name of the non-integer
    non_integer <- names(which(!(sapply(args_integers, function(x) x %% 1)
                                 == 0)))[1]
    stop(paste0(non_integer, " should be integer-valued"))
  }

  # ensure logical input for parameters that are expected to be logicals
  args_logical <- list(kde = kde, scale = scale,
                       monitor_final_elbo = monitor_final_elbo,
                       monitor_cand_elbo = monitor_cand_elbo,
                       print_time = print_time,
                       warnings = warnings)
  if (any(!sapply(args_logical, is.logical))){

    # get the name of the non-logical
    non_logical <- names(which(!sapply(args_logical, is.logical)))[1]
    stop(paste0(non_logical, " is of type ",
                typeof(args_logical[[non_logical]]), "; expected logical"))
  }

  # ensure character input for parameters that are expected to be characters
  args_character <- list(sym_method = sym_method)
  if (any(!sapply(args_character, is.character))){

    # get the name of the non-character
    non_character <- names(which(!sapply(args_character, is.character)))[1]
    stop(paste0(non_character, " is of type ",
                typeof(args_character[[non_character]]),
                "; expected character"))
  }

  # ensure positive input for parameters that are expected to be positive
  args_positive <- list(tau = tau, sigmasq = sigmasq, var_min = var_min,
                        var_max = var_max, n_sigma = n_sigma,
                        tolerance = tolerance, max_iter = max_iter,
                        edge_threshold = edge_threshold,
                        monitor_period = monitor_period)
  if (any(!sapply(args_positive, function (x) all(x > 0)))){

    # get the name of the non-positive
    non_positive <- names(which(!sapply(args_positive,
                                        function (x) all(x > 0))))[1]
    stop(paste0(non_positive, " should have all positive entries"))
  }

  # ensure [0, 1] input for parameters that are expected to be in the interval
  # [0, 1]
  args_01 <- list(alpha = alpha, pi_vec = pi_vec)
  if (any(!sapply(args_01, function (x) all(0 <= x & x <= 1)))){

    # get the name of the non-[0, 1]
    non_01 <- names(which(!sapply(args_01,
                                  function (x) all(0 <= x & x <= 1))))[1]
    stop(paste0(non_01, " should have all entries in the interval [0, 1]"))
  }

  # ensure that data_mat and Z are matrices
  data_mat <- tryCatch(as.matrix(data_mat),
                       error = function(msg){
                         stop(paste("data_mat should be of class matrix or a class that is coercible to a matrix;",
                                    class(data_mat), "is not coercible to a matrix"))})
  Z <- tryCatch(as.matrix(Z),
                error = function(msg){
                  stop(paste("Z should be of class matrix or a class that is coercible to a matrix;",
                             class(Z)[1], "is not coercible to a matrix"))})

  # ensure data_mat and Z have compatible dimensions
  n <- nrow(data_mat)
  if (n != nrow(Z)){
    stop(paste0("Number of rows in data_mat (", n,
                ") is not equal to the number of rows in Z (", nrow(Z), ")"))
  }

  # ensure that tau is either of length 1 or n
  if (!(length(tau) %in% c(1, n))){
    stop(paste0("tau should be of length 1 or ", n, ", not ", length(tau)))
  }

  # check that var_max is greater than var_min
  if (!(var_max > var_min)){
    stop("var_max should be greater than var_min")
  }

  # if the user has specified sigmabeta_sq, ensure that it is of the proper form
  if (!is.null(sigmabetasq_vec)){

    # check that all entries are numeric
    if (!is.numeric(sigmabetasq_vec)){
      stop(paste0("sigmabetasq_vec is of type ", typeof(sigmabetasq_vec),
                  "; expected numeric"))
    }

    # check that it is a vector
    if (!is.vector(sigmabetasq_vec)){
      stop(paste0("sigmabetasq_vec is of class ", class(sigmabetasq_vec)[1],
                  "; expected vector"))
    }

    # check that is does not have any NA entries
    if (any(is.na(sigmabetasq_vec))){
      stop("sigmabetasq_vec should have all non-NA entries")
    }

    # check that all entries are finite
    if (any(is.infinite(sigmabetasq_vec))){
      stop("sigmabetasq_vec should have all finite entries")
    }

    # check that entries are postive
    if (!all(sigmabetasq_vec > 0)){
      stop("sigmabetasq_vec should have all positive entries")
    }
  }

  # ensure that norm is greater than or equal to 1
  if (norm < 1){
    stop("norm should be greater than or equal to 1")
  }

  # if scale is true, ensure that the standard deviation of is covariates is positive
  if ((scale | kde) & isTRUE(all.equal(0, as.numeric(var(Z))))){
    stop("Constant Z; set `scale = F` and `kde = F`")
  }

  # ensure that edge_threshold is in (0, 1)
  if (!(0 < edge_threshold & edge_threshold < 1)){
    stop("edge_threshold should be in the interval (0, 1)")
  }

  # ensure that sym_method is mean, min, or max
  if (!(sym_method %in% c("mean", "min", "max"))){
    stop("sym_method should be one of \"mean\", \"min\", or \"max\"")
  }

  # ensure that monitor_period is not greater than max_iter
  if (monitor_period > max_iter){
    stop("monitor_period is greater than max_iter")
  }
}

## _____________________________________________________________________________
## _____________________________adjMat_checks___________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## function to check compatibility of arguments to gg_adjMat
## -----------------------------ARGUMENTS---------------------------------------
## out: list OR matrix; return of covdepGE function OR adjacency matrix
##
## l: scalar in {1, 2, ..., n}; individual index
##
## prob_shade: logical scalar; if T, entries will be shaded probability-wise
##
## color0: scalar; color for 0 entries
##
## color1: scalar; color for 1 entries
##
## grid_color: scalar; color of grid lines
##
## incl_probs: logical; display posterior inclusion probability
##
## prob_prec: scalar in {1, 2, ...}; number of decimal places to round
##
## font_size: scalar in (0, Inf); size of font if incl_probs = T
##
## font_color0: scalar; color of font for 0 entries if incl_probs = T
##
## font_color1: scalar; color of font for 1 entries if incl_probs = T
##
adjMat_checks <- function(out, l, prob_shade, color0, color1, grid_color,
                          incl_probs, prob_prec, font_size, font_color0,
                          font_color1){

  # ensure scalar input for parameters that are expected to be scalars
  args_scalar <- list(l = l, prob_shade = prob_shade, color0 = color0,
                      color1 = color1, grid_color = grid_color,
                      incl_probs = incl_probs, prob_prec = prob_prec,
                      font_size = font_size, font_color0 = font_color0,
                      font_color1 = font_color1)
  if (any(!(sapply(args_scalar, length) == 1))){

    # get the name of the non-scalar
    non_scalar <- names(which(!(sapply(args_scalar, length) == 1)))[1]
    stop(paste0(non_scalar, " is of length ", length(args_scalar[[non_scalar]]),
                "; expected scalar value"))
  }

  # ensure numeric input for parameters that are expected to be numeric
  args_numeric <- list(l = l, prob_prec = prob_prec, font_size = font_size)
  if (any(!sapply(args_numeric, is.numeric))){

    # get the name of the non-numeric
    non_numeric <- names(which(!sapply(args_numeric, is.numeric)))[1]
    stop(paste0(non_numeric, " is of type ", typeof(args_numeric[[non_numeric]]),
                "; expected numeric"))
  }

  # ensure non-NA input for parameters that are to be non-NA
  args_nonNA <- list(out = out, l = l, prob_shade = prob_shade, color0 = color0,
                     color1 = color1, grid_color = grid_color,
                     incl_probs = incl_probs, prob_prec = prob_prec,
                     font_size = font_size, font_color0 = font_color0,
                     font_color1 = font_color1)
  if (any(sapply(args_nonNA, function (x) any(is.na(x))))){

    # get the name of the NA
    na <- names(which(sapply(args_nonNA, function (x) any(is.na(x)))))[1]
    stop(paste0(na, " should have all non-NA entries"))
  }

  # ensure finite input for parameters that are expected to be finite
  args_finite <- list(l = l, prob_prec = prob_prec, font_size = font_size)
  if (any(!sapply(args_finite, function (x) all(is.finite(x))))){

    # get the name of the non-finite
    non_finite <- names(which(!sapply(args_finite,
                                      function (x) all(is.finite(x)))))[1]
    stop(paste0(non_finite, " should have all finite entries"))
  }

  # ensure logical input for parameters that are expected to be logicals
  args_logical <- list(prob_shade = prob_shade, incl_probs = incl_probs)
  if (any(!sapply(args_logical, is.logical))){

    # get the name of the non-logical
    non_logical <- names(which(!sapply(args_logical, is.logical)))[1]
    stop(paste0(non_logical, " is of type ",
                typeof(args_logical[[non_logical]]), "; expected logical"))
  }

  # ensure valid color input for parameters that are expected to be colors
  # citation:
  # https://stackoverflow.com/a/13290832/10965084
  args_colors <- list(color0 = color0, color1 = color1, grid_color = grid_color,
                      font_color0 = font_color0, font_color1 = font_color1)
  if (any(!sapply(args_colors, function(x){
    tryCatch(is.matrix(grDevices::col2rgb(x)), error = function(msg) F)}))){

    # get the name of the non-color
    non_color <- args_colors[!sapply(args_colors, function(x){
      tryCatch(is.matrix(grDevices::col2rgb(x)),
               error = function(msg) F)})][[1]]
    stop(paste0(non_color, " is not a valid color"))
  }

  # ensure positive input for parameters that are expected to be positive
  args_positive <- list(l = l, prob_prec = prob_prec, font_size = font_size)
  if (any(!sapply(args_positive, function (x) all(x > 0)))){

    # get the name of the non-positive
    non_positive <- names(which(!sapply(args_positive,
                                        function (x) all(x > 0))))[1]
    stop(paste0(non_positive, " should have all positive entries"))
  }

  # ensure out is a list OR a matrix
  if (!(is.list(out) | is.matrix(out))){
    stop(paste0("out is of class ", class(out)[1], "; expected list or matrix"))
  }

  # ensure out has inclusion_probs and graphs if it is a list
  if (!(all(c("inclusion_probs", "graphs") %in% names(out))) & is.list(out)){
    stop(paste0("out should be return of function covdepGE"))
  }

  # ensure out has numeric entries if it is a matrix
  if (!is.numeric(out) & is.matrix(out)){
    stop(paste0("out is of type ", typeof(out), "; expected numeric"))
  }

  # ensure l in 1,...,n
  n <- ifelse(is.list(out), length(out$graphs), 1)
  if (!(l %in% 1:n)){
    stop(paste0("l should be in 1, 2, ..., ", n))
  }

  # ensure prob_prec is an integer
  if ((prob_prec %% 1) != 0){
    stop("prob_prec should be integer-valued")
  }
}

## _____________________________________________________________________________
## _____________________________inclusionCurve_checks___________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## function to check compatibility of arguments to gg_inclusionCurve
## -----------------------------ARGUMENTS---------------------------------------
## out: list; return of covdepGE function
##
## col_idx1: scalar in {1, 2, ..., p + 1}; column index of the first variable
##
## col_idx2: scalar in {1, 2, ..., p + 1}; column index of the second variable
##
## line_type: scalar; ggplot2 line type
##
## line_size: scalar in (0, Inf); thickness of the interpolating line
##
## line_color: scalar; color of interpolating line
##
## point_shape: scalar; shape of the points
##
## point_size: scalar in (0, Inf); size of probability points
##
## point_color: scalar; color of probability points
##
## point_fill: scalar; fill of probability points
##
## sort: logical scalar; sort the subject indices for smooth inclusion curve
##
inclusionCurve_checks <- function(out, col_idx1, col_idx2, line_type, line_size,
                                  line_color, point_shape, point_size,
                                  point_color, point_fill, sort){

  # ensure scalar input for parameters that are expected to be scalars
  args_scalar <- list(col_idx1 = col_idx1, col_idx2 = col_idx2,
                      line_type = line_type, line_size = line_size,
                      line_color = line_color, point_shape = point_shape,
                      point_size = point_size, point_color = point_color,
                      point_fill = point_fill, sort = sort)
  if (any(!(sapply(args_scalar, length) == 1))){

    # get the name of the non-scalar
    non_scalar <- names(which(!(sapply(args_scalar, length) == 1)))[1]
    stop(paste0(non_scalar, " is of length ", length(args_scalar[[non_scalar]]),
                "; expected scalar value"))
  }

  # ensure numeric input for parameters that are expected to be numeric
  args_numeric <- list(col_idx1 = col_idx1, col_idx2 = col_idx2,
                       line_size = line_size, point_size = point_size)
  if (any(!sapply(args_numeric, is.numeric))){

    # get the name of the non-numeric
    non_numeric <- names(which(!sapply(args_numeric, is.numeric)))[1]
    stop(paste0(non_numeric, " is of type ", typeof(args_numeric[[non_numeric]]),
                "; expected numeric"))
  }

  # ensure non-NA input for parameters that are to be non-NA
  args_nonNA <- list(out = out, col_idx1 = col_idx1, col_idx2 = col_idx2,
                     line_type = line_type, line_size = line_size,
                     line_color = line_color, point_shape = point_shape,
                     point_size = point_size, point_color = point_color,
                     point_fill = point_fill, sort = sort)
  if (any(sapply(args_nonNA, function (x) any(is.na(x))))){

    # get the name of the NA
    na <- names(which(sapply(args_nonNA, function (x) any(is.na(x)))))[1]
    stop(paste0(na, " should have all non-NA entries"))
  }

  # ensure finite input for parameters that are expected to be finite
  args_finite <- list(col_idx1 = col_idx1, col_idx2 = col_idx2,
                      line_size = line_size, point_size = point_size)
  if (any(!sapply(args_finite, function (x) all(is.finite(x))))){

    # get the name of the non-finite
    non_finite <- names(which(!sapply(args_finite,
                                      function (x) all(is.finite(x)))))[1]
    stop(paste0(non_finite, " should have all finite entries"))
  }

  # ensure valid color input for parameters that are expected to be colors
  # citation:
  # https://stackoverflow.com/a/13290832/10965084
  args_colors <- list(line_color = line_color, point_color = point_color,
                      point_fill = point_fill)
  if (any(!sapply(args_colors, function(x){
    tryCatch(is.matrix(grDevices::col2rgb(x)), error = function(msg) F)}))){

    # get the name of the non-color
    non_color <- args_colors[!sapply(args_colors, function(x){
      tryCatch(is.matrix(grDevices::col2rgb(x)),
               error = function(msg) F)})][[1]]
    stop(paste0(non_color, " is not a valid color"))
  }

  # ensure positive input for parameters that are expected to be positive
  args_positive <- list(col_idx1 = col_idx1, col_idx2 = col_idx2,
                        line_size = line_size, point_size = point_size)
  if (any(!sapply(args_positive, function (x) all(x > 0)))){

    # get the name of the non-positive
    non_positive <- names(which(!sapply(args_positive,
                                        function (x) all(x > 0))))[1]
    stop(paste0(non_positive, " should have all positive entries"))
  }

  # ensure out is a list
  if (!is.list(out)){
    stop(paste0("out is of class ", class(out)[1], "; expected list"))
  }

  # ensure out has inclusion_probs and graphs
  if (!(all(c("inclusion_probs", "graphs") %in% names(out)))){
    stop(paste0("out should be return of function covdepGE"))
  }

  # ensure sort is a logical
  if (!is.logical(sort)){
    stop(paste0("sort is of type ", typeof(sort), "; expected logical"))
  }

  # ensure col_idx1 and col_idx2 are in 1,...,p + 1
  p <- ncol(out$graphs[[1]])
  if (!(col_idx1 %in% 1:p)){
    stop(paste0("col_idx1 should be in 1, 2, ..., ", p))
  }
  if (!(col_idx2 %in% 1:p)){
    stop(paste0("col_idx2 should be in 1, 2, ..., ", p))
  }
}

## _____________________________________________________________________________
## _____________________________gg_adjMats_checks_______________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## function to check compatibility of arguments to gg_adjMats
## -----------------------------ARGUMENTS---------------------------------------
## out: list; return of covdepGE function.
##
## graph_colors: g x 1 vector; g is the number of unique graphs from out. The
## v-th element vector is the color for the v-th unique graph
##
## seed: scalar in (-Inf, Inf); when colors is NULL, the RNG seed for selecting
## the color for each graph. 1 by default.
##
adjMats_checks <- function(out, seed, graph_colors){

  # ensure out is a list
  if (!is.list(out)){
    stop(paste0("out is of class ", class(out)[1], "; expected list"))
  }

  # ensure out has inclusion_probs and graphs
  if (!(all(c("inclusion_probs", "graphs") %in% names(out)))){
    stop(paste0("out should be return of function covdepGE"))
  }

  # if graph_colors is not null, run compatibility checks
  if(!is.null(graph_colors)){

    # ensure graph_colors is a vector
    if (!is.vector(graph_colors)){
      stop(paste0(graph_colors, " is of class ", class(args_vector[[graph_colors]])[1],
                  "; expected vector"))
    }

    # ensure that there are a sufficient number of colors in graph_colors
    if(length(graph_colors) < length(out$unique_graphs)){
      stop(paste0("graph_colors is of length ", length(graph_colors),
                  "; should have at least ", length(out$unique_graphs), " elements"))
    }
  }

  # ensure seed is numeric
  if (!is.numeric(seed)){
    stop(paste0("seed is of type ", typeof(seed), "; expected numeric"))
  }

  # ensure seed is finite
  if (!all(is.finite(seed))){
    stop(paste0("seed should have all finite entries"))
  }
}

