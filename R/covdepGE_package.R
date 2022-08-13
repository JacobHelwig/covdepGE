#' @aliases covdepGE-package
#'
#' @details The core function, `covdepGE`, uses the weighted pseudo-likelihood
#' approach to estimate the conditional dependence structure of the data as a
#' function of an extraneous covariate. Inference is conducted efficiently via
#' a parallelized block mean-field variational approximation. Three choices
#' for hyperparameter specification are offered, the default being a hybrid
#' between model averaging and grid search.
#'
#' Additionally, the function `generateData` returns covariate dependent data
#' based on the data from the simulation study in (1). The functions
#' `inclusionCurve`, `matViz`, and `plot.covdepGE` enable visualization of the
#' estimates returned by the `covdepGE` function.
#'
#' @references
#' (1) Sutanoy Dasgupta, Peng Zhao, Prasenjit Ghosh, Debdeep Pati, and Bani
#' Mallick. An approximate Bayesian approach to covariate-dependent graphical
#' modeling. pages 1–59, 2022.
#'
"_PACKAGE"

## usethis namespace: start
#' @useDynLib covdepGE, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
