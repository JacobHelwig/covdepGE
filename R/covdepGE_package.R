#' @aliases covdepGE-package
#'
#' @details The conditional dependence structure (CDS) of a data matrix with
#' \eqn{p} variables can be modeled as an undirected graph with \eqn{p}
#' vertices, where two variables are connected if, and only if, the two
#' variables are dependent given the remaining variables in the data. Gaussian
#' graphical modeling (GGM) seeks to capture the CDS of the data under the
#' assumption that the data are normally distributed. This distributional
#' assumption is convenient for inference, as the CDS is given by the sparsity
#' structure of the precision matrix, where the precision matrix is defined as
#' the inverse covariance matrix of the data.
#'
#' There is extensive GGM literature and many R packages for GGM, however, all
#' make the restrictive assumption that the precision matrix is homogeneous
#' throughout the data, or that there exists a partition of homogeneous
#' subgroups. `covdepGE` avoids this strong assumption by utilizing information
#' sharing to model the CDS as varying continuously with an extraneous
#' covariate. Intuitively, this implies that observations having similar
#' extraneous covariate values will have similar precision matrices.
#'
#' To facilitate information sharing while managing complexity, `covdepGE` uses
#' an efficient variational approximation conducted under the novel weighted
#' pseudo-likelihood framework proposed by (1). `covdepGE` further accelerates
#' inference by employing parallelism and executing expensive iterative
#' computations in  C++. Additionally, `covdepGE` offers a principled,
#' data-driven approach for hyperparameter specification that only requires the
#' user to input data and extraneous covariates to perform inference. Finally,
#' `covdepGE` offers several wrappers around `ggplot2` for seamless
#' visualization of resulting estimates, such as `matViz`, `inclusionCurve`, and
#' the S3 method `plot.covdepGE`.
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
