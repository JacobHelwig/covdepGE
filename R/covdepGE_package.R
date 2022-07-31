#' @aliases covdepGE-package
#'
#' @details The core function, `covdepGE`, uses the weighted psuedo-likelihood
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
#' (1) Dasgupta S., Zhao P., Ghosh P., Pati D., Mallick B., *An approximate
#' Bayesian approach to covariate dependent graphical modeling*, 2021
#'
"_PACKAGE"

## usethis namespace: start
#' @useDynLib covdepGE, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
