## -----------------------------------------------------------------------------
#' @title matViz
#' @export
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' @description Create a visualization of a matrix
## -----------------------------ARGUMENTS---------------------------------------
#' @param x `matrix`; `matrix` to be visualized
#'
#' @param color1 color; color for low entries. `"white"` by default
#'
#' @param color2 color; color for high entries. `"#500000"` by default
#'
#' @param shade `logical`; if `T`, then entries will be shaded according to value
#' on a gradient ranging from `color_low` for the least value to `color_high`
#' for the greatest value
#'
#' @param grid_color color; color of grid lines. `"black"` by default
#'
#' @param incl_val logical; if `T`, the value for each entry will be displayed.
#' `F` by default
#'
#' @param prec positive integer; number of decimal places to round entries to if
#' `incl_val` is `T`. `2` by default
#'
#' @param font_size positive `numeric`; size of font if `incl_val` is `T`. `3`
#' by default
#'
#' @param font_color1 color; color of font for low entries if `incl_val` is `T`.
#' `"black"` by default
#'
#' @param font_color2 color; color of font for high entries if `incl_val` is
#'  `T`. `"white"` by default
#'
#' @param font_thres `numeric`; values less than `font_thres` will be displayed
#' in `font_color1` if `incl_val` is `T`. `mean(x)` by default
## -----------------------------RETURNS-----------------------------------------
#' @return Returns `ggplot2` visualization of `matrix`
## -----------------------------EXAMPLES----------------------------------------
#' @examples
## -----------------------------------------------------------------------------
matViz <- function(x, color1 = "white", color2 = "#500000", shade = F,
                   grid_color = "black", incl_val = F, prec = 2, font_size = 3,
                   font_color1 = "black", font_color2 = "white",
                   font_thres = mean(x)){

  # verify that `x` is a matrix
  if (!is.matrix(x)){
    stop(paste0("x is of class ", class(x), " but should be of class matrix"))
  }

  # takes care of no visible bindings
  Var1 <- Var2 <- value <- NULL
  rm(list = c("Var1", "Var2", "value"))

  # save col and row names and remove them
  colnames <- rownames <- NULL
  colnames <- colnames(x)
  rownames <- rownames(x)
  x <- unname(x)

  # melt to long form
  long_graph <- reshape2::melt(x)

  # shade on a gradient according to value
  if (shade){

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
        value, prec), color = graph), show.legend = F, size = font_size) +
          ggplot2::scale_color_manual(values = c(font_color1, font_color2)))
    }

  }else{

    # otherwise, use binary shading; verify that x has 2 or less unique
    # values
    if (length(unique(as.vector(x))) > 2){
      stop("Set `shade = T` if x has greater than 2 unique values")
    }

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
#' @title inclusionCurve
#' @export
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' @description Create a visualization of the probabilities of inclusion of an
#' edge between two variables across all `n` individuals
## -----------------------------ARGUMENTS---------------------------------------
#' @param out object of class `covdepGE`; return of `covdepGE` function
#'
#' @param col_idx1 integer in `[1, p]`; column index of the first variable
#'
#' @param col_idx2 integer in `[1, p]`; column index of the second variable
#'
#' @param line_type linetype; `ggplot2` line type to interpolate the
#' probabilities. `"solid"` by default
#'
#' @param line_size positive `numeric`; thickness of the interpolating line.
#' `0.5` by default
#'
#' @param line_color color; color of interpolating line. `"black"` by default
#'
#' @param point_shape shape; shape of the points denoting individual-specific
#' inclusion probabilities; `21` by default
#'
#' @param point_size positive `numeric`; size of probability points. `1.5` by
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

  # get the probabilities for each individual of an edge between the variables
  # corresponding to col_idx1 and col_idx2
  prob1_2 <- sapply(out$graphs$inclusion_probs_sym,
                    function(x) x[col_idx1, col_idx2])
  prob1_2 <- data.frame(idx = 1:length(prob1_2), prob = c(prob1_2))

  vis <- (ggplot2::ggplot(prob1_2, ggplot2::aes(idx, prob)) +
            ggplot2::geom_line(
              linetype = line_type, size = line_size, color = line_color) +
            ggplot2::geom_point(shape = point_shape, size = point_size,
                                color = point_color, fill = point_fill) +
            ggplot2::xlab("Subject index") +
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
#' @title plot.covdepGE
#' @export
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' @description Given the return value of `covdepGE` function, create a list of
#'  visualizations of the adjacency matrix for each of the unique graphs
## -----------------------------ARGUMENTS---------------------------------------
#' @param x object of class `covdepGE`; return of `covdepGE` function
#'
#' @param graph_colors `NULL` OR `vector` of length `g`; `g` is the number of
#' unique graphs from `x`. The `v`-th element is the color for the `v`-th unique
#' graph. If `NULL`:
#'
#' `graph_colors <- rep("#500000", length(out$unique_graphs))`
#'
#' `NULL` by default
#'
#' @param title_sum `logical`; whether the indices of the individuals
#' corresponding to the graph should be included in the title. `T` by default
#'
#' @param ... additional arguments will be ignored
## -----------------------------RETURNS-----------------------------------------
#' @return Returns list of `ggplot2` visualizations of unique adjacency matrices
#' estimated by `covdepGE`
## -----------------------------EXAMPLES----------------------------------------
#' @examples
## -----------------------------------------------------------------------------
plot.covdepGE <- function(x, graph_colors = NULL, title_sum = T, ...){

  # verify that out is of class covdepGE
  if (class(x)[1] != "covdepGE"){
    stop(paste0("x is of class ", class(x)[1], "; expected covdepGE"))
  }

  # if no colors have been provided, set to default
  if(is.null(graph_colors)){
    graph_colors <- rep("#500000", length(x$graphs$unique_graphs))
  }

  # create the titles for the plots
  titles <- paste("Graph", 1:length(x$graphs$unique_graphs))

  # check if the title should include the individual's summaries
  if(title_sum){
    indv_sum <- sapply(x$graphs$unique_graphs, `[[`, "individuals_summary")
    titles <- paste0(titles, ", Individuals ", indv_sum)
  }

  # get the unique graphs
  unique_graphs <- lapply(x$graphs$unique_graphs, `[[`,"graph")

  # create a visualization for each graph and store it in a list
  graph_viz <- lapply(1:length(unique_graphs), function(gr_idx) matViz(
    unique_graphs[[gr_idx]], color2 = graph_colors[gr_idx]) + ggplot2::ggtitle(
      titles[gr_idx]))

  return(graph_viz)
}
