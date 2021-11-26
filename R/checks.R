## _____________________________________________________________________________
## _____________________________covdepGE_checks_________________________________
## _____________________________________________________________________________
## -----------------------------DESCRIPTION-------------------------------------
## Function to check compatiibility of arguments of covdepGE
## -----------------------------ARGUMENTS---------------------------------------
## data_mat: n x (p + 1) matrix; data
## Z: n x p' matrix; extraneous covariates
## tau: scalar in (0, Inf) OR n x 1 vector, entries in (0, Inf)
## kde: logical; if T, use 2-step KDE methodology
## alpha: scalar in [0, 1]; global initialization value
## mu: scalar; global initialization value
## sigmasq: scalar in (0, Inf); variance hyperparameter for spike-and-slab.
## sigmabetasq_vec: n_sigma x 1 vector, entries in (0, Inf); candidate values
## var_min: scalar in (0, Inf); lower bound of sigmabetasq_vec
## var_max: scalar in (0, Inf); upper bound of sigmabetasq_vec
## n_sigma: scalar in {1, 2,...}
## pi_vec: n_pi x 1 vector, entries in [0, 1]; candidate values of pi
## norm: scalar in [1, Inf]; norm to use when calculating weights
## scale: logical; if T, center and scale extraneous covariates
## tolerance: scalar in (0, Inf); end the variational update loop
## max_iter: scalar in {1, 2,...}; end the variational update loop
## edge_threshold: scalar in (0, 1)
## sym_method: string in {"mean", "max", "min"}
## print_time: logical; if T, function run time is printed
## warnings: logical; if T, convergence and grid warnings will be displayed
covdepGE_checks <- function(data_mat, Z, tau, kde, alpha, mu, sigmasq,
                            sigmabetasq_vec, var_min, var_max, n_sigma, pi_vec,
                            norm, scale, tolerance, max_iter, edge_threshold,
                            sym_method, print_time, warnings){

  # ensure that data_mat and Z are matrices
  data_mat <- tryCatch(as.matrix(data_mat),
                       error = function(msg){
                         stop(paste("data_mat should be of class matrix or a class that is coercible to a matrix;",
                                    class(data_mat), "is not coercible to a matrix"))})
  Z <- tryCatch(as.matrix(Z),
                error = function(msg){
                  stop(paste("Z should be of class matrix or a class that is coercible to a matrix;",
                             class(Z)[1], "is not coercible to a matrix"))})

  # check dimensions of data_mat and Z for compatibility
  n <- nrow(data_mat)
  if (n != nrow(Z)){
    stop(paste0("Number of rows in data_mat (", n,
                ") is not equal to the number of rows in Z (", nrow(Z), ")"))
  }

  # ensure numeric input for parameters that are expected to be numeric
  args_numeric <- list(data_mat = data_mat, Z = Z, tau = tau, alpha = alpha,
                       mu = mu, sigmasq = sigmasq, var_min = var_min,
                       n_sigma = n_sigma, var_max = var_max, pi_vec = pi_vec,
                       norm = norm, tolerance = tolerance, max_iter = max_iter,
                       edge_threshold = edge_threshold)
  if (any(!sapply(args_numeric, is.numeric))){

    # get the name of the non-numeric
    non_numeric <- names(which(!sapply(args_numeric, is.numeric)))[1]
    stop(paste0(non_numeric, " is of type ",
                typeof(args_numeric[[non_numeric]]),
                "; expected numeric"))
  }

  # ensure non-NA input for parameters that are to be non-NA
  args_nonNA <- list(data_mat = data_mat, Z = Z, tau = tau, kde = kde,
                     alpha = alpha, mu = mu, sigmasq = sigmasq,
                     var_min = var_min, var_max = var_max, n_sigma= n_sigma,
                     pi_vec = pi_vec, norm = norm, scale = scale,
                     tolerance = tolerance, max_iter = max_iter,
                     edge_threshold = edge_threshold, sym_method = sym_method,
                     print_time = print_time, warnings = warnings)
  if (any(sapply(args_nonNA, function (x) any(is.na(x))))){

    # get the name of the NA
    na <- names(which(sapply(args_nonNA, function (x) any(is.na(x)))))[1]
    stop(paste0(na, " should have all non-NA entries"))
  }

  # ensure finite input for parameters that are expected to be finite
  args_finite <- list(data_mat = data_mat, Z = Z, tau = tau, mu = mu,
                      sigmasq = sigmasq, var_min = var_min, var_max = var_max,
                      n_sigma = n_sigma, tolerance = tolerance,
                      max_iter = max_iter)
  if (any(!sapply(args_finite, function (x) all(is.finite(x))))){

    # get the name of the non-finite
    non_finite <- names(which(!sapply(args_finite,
                                      function (x) all(is.finite(x)))))[1]
    stop(paste0(non_finite, " should have all finite entries"))
  }

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
                      sym_method = sym_method, print_time = print_time,
                      warnings = warnings)
  if (any(!(sapply(args_scalar, length) == 1))){

    # get the name of the non-scalar
    non_scalar <- names(which(!(sapply(args_scalar, length) == 1)))[1]
    stop(paste0(non_scalar, " is of length ",
                length(args_scalar[[non_scalar]]),
                "; expected scalar value"))
  }

  # ensure integer input for parameters that are expected to be integers
  args_integers <- list(n_sigma = n_sigma, max_iter = max_iter)
  if (any(!(sapply(args_integers, function(x) x %% 1) == 0))){

    # get the name of the non-integer
    non_integer <- names(which(!(sapply(args_integers, function(x) x %% 1)
                                 == 0)))[1]
    stop(paste0(non_integer, " should be integer-valued"))
  }

  # ensure logical input for parameters that are expected to be logicals
  args_logical <- list(kde = kde, scale = scale, print_time = print_time,
                       warnings = warnings)
  if (any(!sapply(args_logical, is.logical))){

    # get the name of the non-logical
    non_logical <- names(which(!sapply(args_logical, is.logical)))[1]
    stop(paste0(non_logical, " is of type ",
                typeof(args_logical[[non_logical]]),
                "; expected logical"))
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
                        edge_threshold = edge_threshold)
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

  # check that norm is greater than or equal to 1
  if (norm < 1){
    stop("norm should be greater than or equal to 1")
  }

  # check that edge_threshold is in (0, 1)
  if (!(0 < edge_threshold & edge_threshold < 1)){
    stop("edge_threshold should be in the interval (0, 1)")
  }

  # check that sym_method is mean, min, or max
  if (!(sym_method %in% c("mean", "min", "max"))){
    stop("sym_method should be one of \"mean\", \"min\", or \"max\"")
  }
}
