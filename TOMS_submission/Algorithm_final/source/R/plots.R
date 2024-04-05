## -----------------------------------------------------------------------------
## Distributed under GPL (≥ 3) license
#'
#' @title Visualize a matrix
#' @export
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' @description Create a visualization of a matrix
## -----------------------------ARGUMENTS---------------------------------------
#' @param x matrix; matrix to be visualized
#'
#' @param color1 color; color for low entries. `"white"` by default
#'
#' @param color2 color; color for high entries. `"#500000"` by default
#'
#' @param grid_color color; color of grid lines. `"black"` by default
#'
#' @param incl_val logical; if `TRUE`, the value for each entry will be
#' displayed. `FALSE` by default
#'
#' @param prec positive integer; number of decimal places to round entries to if
#' `incl_val` is `TRUE`. `2` by default
#'
#' @param font_size positive numeric; size of font if `incl_val` is `TRUE`. `3`
#' by default
#'
#' @param font_color1 color; color of font for low entries if `incl_val` is
#' `TRUE`. `"black"` by default
#'
#' @param font_color2 color; color of font for high entries if `incl_val` is
#' `TRUE`. `"white"` by default
#'
#' @param font_thres numeric; values less than `font_thres` will be displayed
#' in `font_color1` if `incl_val` is `TRUE`. `mean(x)` by default
## -----------------------------RETURNS-----------------------------------------
#' @return Returns `ggplot2` visualization of matrix
## -----------------------------EXAMPLES----------------------------------------
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # get the data
#' set.seed(12)
#' data <- generateData()
#' X <- data$X
#' Z <- data$Z
#' interval <- data$interval
#' prec <- data$true_precision
#'
#' # get overall and within interval sample sizes
#' n <- nrow(X)
#' n1 <- sum(interval == 1)
#' n2 <- sum(interval == 2)
#' n3 <- sum(interval == 3)
#'
#' # visualize the distribution of the extraneous covariate
#' ggplot(data.frame(Z = Z, interval = as.factor(interval))) +
#'   geom_histogram(aes(Z, fill = interval), color = "black", bins = n %/% 5)
#'
#' # visualize the true precision matrices in each of the intervals
#'
#' # interval 1
#' matViz(prec[[1]], incl_val = TRUE) +
#'   ggtitle(paste0("True precision matrix, interval 1, observations 1,...,", n1))
#'
#' # interval 2 (varies continuously with Z)
#' cat("\nInterval 2, observations ", n1 + 1, ",...,", n1 + n2, sep = "")
#' int2_mats <- prec[interval == 2]
#' int2_inds <- c(5, n2 %/% 2, n2 - 5)
#' lapply(int2_inds, function(j) matViz(int2_mats[[j]], incl_val = TRUE) +
#'          ggtitle(paste("True precision matrix, interval 2, observation", j + n1)))
#'
#' # interval 3
#' matViz(prec[[length(prec)]], incl_val = TRUE) +
#'   ggtitle(paste0("True precision matrix, interval 3, observations ",
#'                  n1 + n2 + 1, ",...,", n1 + n2 + n3))
#'
#' # fit the model and visualize the estimated graphs
#' (out <- covdepGE(X, Z))
#' plot(out)
#'
#' # visualize the posterior inclusion probabilities for variables (1, 3) and (1, 2)
#' inclusionCurve(out, 1, 2)
#' inclusionCurve(out, 1, 3)
#' }
## -----------------------------------------------------------------------------
matViz <- function(x, color1 = "white", color2 = "#500000",
                   grid_color = "black", incl_val = FALSE, prec = 2,
                   font_size = 3, font_color1 = "black", font_color2 = "white",
                   font_thres = mean(x)){

  # verify that `x` is a matrix
  if (!is.matrix(x)){
    stop(paste0("x is of class ", class(x), " but should be of class matrix"))
  }

  # takes care of no visible bindings
  Var1 <- Var2 <- value <- graph <- NULL
  rm(list = c("Var1", "Var2", "value", "graph"))

  # save col and row names and remove them
  colnames <- rownames <- NULL
  colnames <- colnames(x)
  rownames <- rownames(x)
  x <- unname(x)

  # melt to long form
  long_graph <- reshape2::melt(x)

  # shade on a gradient according to value
  if (length(unique(c(x))) > 2){

    # if values are to be displayed, add an indicator to edges for the cells that
    # will be darker
    if (incl_val){
      long_graph$graph <- (long_graph$value >= font_thres) * 1
      long_graph$graph <- as.factor(long_graph$graph)
    }

    vis <- (ggplot2::ggplot(long_graph, ggplot2::aes(Var1, Var2, fill = value)) +
              ggplot2::geom_tile(color = grid_color) +
              ggplot2::scale_fill_gradient(low = color1, high = color2) +
              ggplot2::theme_classic() + ggplot2::xlab("") + ggplot2::ylab("") +
              ggplot2::ggtitle("") +
              ggplot2::theme(legend.title = ggplot2::element_blank(),
                             plot.title = ggplot2::element_text(hjust = 0.5)) +
              ggplot2::scale_x_continuous(breaks = 1:nrow(x)) +
              ggplot2::scale_y_continuous(breaks = 1:nrow(x)))

    # if probabilities are to be displayed, add them
    if (incl_val){
      vis <- (vis + ggplot2::geom_text(ggplot2::aes(label = round(
        value, prec), color = graph), show.legend = FALSE, size = font_size) +
          ggplot2::scale_color_manual(values = c(font_color1, font_color2)))
    }

  }else{

    # otherwise, x has less than 3 unique values; use binary shading

    # factor the edges - specifies a discrete scale to ggplot2
    long_graph$value <- as.factor(long_graph$value)

    vis <- (ggplot2::ggplot(long_graph, ggplot2::aes(Var1, Var2, fill = value)) +
              ggplot2::geom_tile(color = grid_color) +
              ggplot2::scale_fill_manual(values = c(color1, color2)) +
              ggplot2::theme_classic() + ggplot2::xlab("") + ggplot2::ylab("") +
              ggplot2::ggtitle("") +
              ggplot2::theme(legend.title = ggplot2::element_blank(),
                             plot.title = ggplot2::element_text(hjust = 0.5)) +
              ggplot2::scale_x_continuous(breaks = 1:nrow(x)) +
              ggplot2::scale_y_continuous(breaks = 1:nrow(x)))
  }

  # check to see if there are col names or row names to be added back
  if (!is.null(colnames)){
    vis <- suppressMessages(vis + ggplot2::scale_x_continuous(
      breaks = 1:nrow(x), labels = colnames))
  }
  if (!is.null(rownames)){
    vis <- suppressMessages(vis + ggplot2::scale_y_continuous(
      breaks = 1:nrow(x), labels = rownames))
  }

  return(vis)
}

## -----------------------------------------------------------------------------
## Distributed under GPL (≥ 3) license
#'
#' @title Plot PIP as a Function of Index
#' @export
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' @description Plot the posterior inclusion probability of an edge between two
#' variables as a function of observation index
## -----------------------------ARGUMENTS---------------------------------------
#' @param out object of class `covdepGE`; return of `covdepGE` function
#'
#' @param col_idx1 integer in \eqn{[1, p]}; column index of the first variable
#'
#' @param col_idx2 integer in \eqn{[1, p]}; column index of the second variable
#'
#' @param line_type linetype; `ggplot2` line type to interpolate the
#' probabilities. `"solid"` by default
#'
#' @param line_size positive numeric; thickness of the interpolating line.
#' `0.5` by default
#'
#' @param line_color color; color of interpolating line. `"black"` by default
#'
#' @param point_shape shape; shape of the points denoting observation-specific
#' inclusion probabilities; `21` by default
#'
#' @param point_size positive numeric; size of probability points. `1.5` by
#' default
#'
#' @param point_color color; color of probability points. `"#500000"` by default
#'
#' @param point_fill color; fill of probability points. Only applies to select
#' shapes. `"white"` by default
## -----------------------------RETURNS-----------------------------------------
#' @return Returns `ggplot2` visualization of inclusion probability curve
## -----------------------------EXAMPLES----------------------------------------
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # get the data
#' set.seed(12)
#' data <- generateData()
#' X <- data$X
#' Z <- data$Z
#' interval <- data$interval
#' prec <- data$true_precision
#'
#' # get overall and within interval sample sizes
#' n <- nrow(X)
#' n1 <- sum(interval == 1)
#' n2 <- sum(interval == 2)
#' n3 <- sum(interval == 3)
#'
#' # visualize the distribution of the extraneous covariate
#' ggplot(data.frame(Z = Z, interval = as.factor(interval))) +
#'   geom_histogram(aes(Z, fill = interval), color = "black", bins = n %/% 5)
#'
#' # visualize the true precision matrices in each of the intervals
#'
#' # interval 1
#' matViz(prec[[1]], incl_val = TRUE) +
#'   ggtitle(paste0("True precision matrix, interval 1, observations 1,...,", n1))
#'
#' # interval 2 (varies continuously with Z)
#' cat("\nInterval 2, observations ", n1 + 1, ",...,", n1 + n2, sep = "")
#' int2_mats <- prec[interval == 2]
#' int2_inds <- c(5, n2 %/% 2, n2 - 5)
#' lapply(int2_inds, function(j) matViz(int2_mats[[j]], incl_val = TRUE) +
#'          ggtitle(paste("True precision matrix, interval 2, observation", j + n1)))
#'
#' # interval 3
#' matViz(prec[[length(prec)]], incl_val = TRUE) +
#'   ggtitle(paste0("True precision matrix, interval 3, observations ",
#'                  n1 + n2 + 1, ",...,", n1 + n2 + n3))
#'
#' # fit the model and visualize the estimated graphs
#' (out <- covdepGE(X, Z))
#' plot(out)
#'
#' # visualize the posterior inclusion probabilities for variables (1, 3) and (1, 2)
#' inclusionCurve(out, 1, 2)
#' inclusionCurve(out, 1, 3)
#' }
## -----------------------------------------------------------------------------
inclusionCurve <- function(out, col_idx1, col_idx2, line_type = "solid",
                           line_size = 0.5, line_color = "black",
                           point_shape = 21, point_size = 1.5,
                           point_color = "#500000", point_fill = "white"){

  # verify that out is of class covdepGE
  if (class(out)[1] != "covdepGE"){
    stop(paste0("out is of class ", class(out)[1], "; expected covdepGE"))
  }

  # takes care of "no visible bindings" note
  idx <- prob <- NULL
  rm(list = c("idx", "prob"))

  # get the probabilities for each observation of an edge between the variables
  # corresponding to col_idx1 and col_idx2
  prob1_2 <- sapply(out$graphs$inclusion_probs_sym,
                    function(x) x[col_idx1, col_idx2])
  prob1_2 <- data.frame(idx = 1:length(prob1_2), prob = c(prob1_2))

  vis <- (ggplot2::ggplot(prob1_2, ggplot2::aes(idx, prob)) +
            ggplot2::geom_line(
              linetype = line_type, size = line_size, color = line_color) +
            ggplot2::geom_point(shape = point_shape, size = point_size,
                                color = point_color, fill = point_fill) +
            ggplot2::xlab("Observation Index") +
            ggplot2::ylab("Posterior Inclusion Probability") +
            ggplot2::ggtitle(latex2exp::TeX(paste(
              "Inclusion probability of an edge between $x_", col_idx1,
              "$ and $x_", col_idx2, "$"))) +
            ggplot2::theme_bw() +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::ylim(c(0, 1)))

  return(vis)
}

## -----------------------------------------------------------------------------
## Distributed under GPL (≥ 3) license
#'
#' @title Plot the Graphs Estimated by `covdepGE`
#' @export
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' @description Create a list of the unique graphs estimated by `covdepGE`
## -----------------------------ARGUMENTS---------------------------------------
#' @param x object of class `covdepGE`; return of `covdepGE` function
#'
#' @param graph_colors `NULL` OR vector; the \eqn{j}-th element is the color for
#' the \eqn{j}-th graph. If `NULL`, all graphs will be colored with `"#500000"`.
#' `NULL` by default
#'
#' @param title_sum logical; if `TRUE` the indices of the observations
#' corresponding to the graph will be included in the title. `TRUE` by default
#'
#' @param ... additional arguments will be ignored
## -----------------------------RETURNS-----------------------------------------
#' @return Returns list of `ggplot2` visualizations of unique graphs estimated
#' by `covdepGE`
## -----------------------------EXAMPLES----------------------------------------
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # get the data
#' set.seed(12)
#' data <- generateData()
#' X <- data$X
#' Z <- data$Z
#' interval <- data$interval
#' prec <- data$true_precision
#'
#' # get overall and within interval sample sizes
#' n <- nrow(X)
#' n1 <- sum(interval == 1)
#' n2 <- sum(interval == 2)
#' n3 <- sum(interval == 3)
#'
#' # visualize the distribution of the extraneous covariate
#' ggplot(data.frame(Z = Z, interval = as.factor(interval))) +
#'   geom_histogram(aes(Z, fill = interval), color = "black", bins = n %/% 5)
#'
#' # visualize the true precision matrices in each of the intervals
#'
#' # interval 1
#' matViz(prec[[1]], incl_val = TRUE) +
#'   ggtitle(paste0("True precision matrix, interval 1, observations 1,...,", n1))
#'
#' # interval 2 (varies continuously with Z)
#' cat("\nInterval 2, observations ", n1 + 1, ",...,", n1 + n2, sep = "")
#' int2_mats <- prec[interval == 2]
#' int2_inds <- c(5, n2 %/% 2, n2 - 5)
#' lapply(int2_inds, function(j) matViz(int2_mats[[j]], incl_val = TRUE) +
#'          ggtitle(paste("True precision matrix, interval 2, observation", j + n1)))
#'
#' # interval 3
#' matViz(prec[[length(prec)]], incl_val = TRUE) +
#'   ggtitle(paste0("True precision matrix, interval 3, observations ",
#'                  n1 + n2 + 1, ",...,", n1 + n2 + n3))
#'
#' # fit the model and visualize the estimated graphs
#' (out <- covdepGE(X, Z))
#' plot(out)
#'
#' # visualize the posterior inclusion probabilities for variables (1, 3) and (1, 2)
#' inclusionCurve(out, 1, 2)
#' inclusionCurve(out, 1, 3)
#' }
## -----------------------------------------------------------------------------
plot.covdepGE <- function(x, graph_colors = NULL, title_sum = TRUE, ...){

  # if no colors have been provided, set to default
  if(is.null(graph_colors)){
    graph_colors <- rep("#500000", length(x$graphs$unique_graphs))
  }else if(length(graph_colors) < x$model_details$num_unique){

    # otherwise, ensure that enough colors have been provided
    graph_colors <- c(matrix(graph_colors, x$model_details$num_unique, 1))
  }

  # create the titles for the plots
  titles <- paste("Graph", 1:length(x$graphs$unique_graphs))

  # check if the title should include the observation's summaries
  if(title_sum){
    obs_sum <- sapply(x$graphs$unique_graphs, `[[`, "ind_sum")
    titles <- paste0(titles, ", observations ", obs_sum)
  }

  # get the unique graphs
  unique_graphs <- lapply(x$graphs$unique_graphs, `[[`,"graph")

  # create a visualization for each graph and store it in a list
  graph_viz <- lapply(1:length(unique_graphs), function(gr_idx) matViz(
    unique_graphs[[gr_idx]], color2 = graph_colors[gr_idx]) + ggplot2::ggtitle(
      titles[gr_idx]))

  return(graph_viz)
}
