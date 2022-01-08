## -----------------------------------------------------------------------------
#' @title gg_adjMat
#' @export
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' @description Create a visualization of an adjacency matrix
## -----------------------------ARGUMENTS---------------------------------------
#' @param object of class `covdepGE` OR `matrix`; return of `covdepGE` function
#' OR an adjacency matrix
#'
#' @param l integer in \eqn{{1, 2, ..., n}}; individual index for which the
#' adjacency matrix will be visualized. Ignored if `out` is a
#' `matrix`. `1` by default
#'
#' @param prob_shade logical; if `T`, then entries will be shaded according to
#' posterior inclusion probabilities on a gradient ranging from `color0`
#' (probability 0) to `color1` (probability 1); if `F`, binary coloring is used.
#' Ignored if `out` is a `matrix`. `T` by default
#'
#' @param color0 color; color for 0 entries. `"white"` by default
#'
#' @param color1 color; color for 1 entries. `"#500000"` by default
#'
#' @param grid_color color; color of grid lines. `"black"` by default
#'
#' @param incl_probs logical; if `T`, the posterior inclusion probability
#' will be displayed for each entry. Ignored if `out` is a `matrix` with 2 or
#' less unique entries. `T` by default
#'
#' @param prob_prec positive integer; number of decimal places to round
#' probabilities to if `incl_probs = T`. Ignored if `out` is a `matrix` with 2
#' or less unique entries. `2` by default
#'
#' @param font_size positive numeric; size of font if `incl_probs = T`. Ignored
#' if `out` is a `matrix` with 2 or less unique entries. `3` by default
#'
#' @param font_color0 color; color of font for 0 entries if `incl_probs = T`.
#' Ignored if `out` is a `matrix` with 2 or less unique entries. `"black"` by
#' default
#'
#' @param font_color1 color; color of font for 1 entries if `incl_probs = T`.
#' Ignored if `out` is a `matrix` with 2 or less unique entries. `"white"` by
#' default
## -----------------------------RETURNS-----------------------------------------
#' @return Returns `ggplot2` visualization of adjacency matrix
## -----------------------------EXAMPLES----------------------------------------
#' @examples
#' set.seed(1)
#' n <- 100
#' p <- 4
#'
#' # generate the extraneous covariate
#' Z_neg <- sort(runif(n / 2) * -1)
#' Z_pos <- sort(runif(n / 2))
#' Z <- c(Z_neg, Z_pos)
#' summary(Z)
#'
#' # create true covariance structure for 2 groups: positive Z and negative Z
#' true_graph_pos <- true_graph_neg <- matrix(0, p + 1, p + 1)
#' true_graph_pos[1, 2] <- true_graph_pos[2, 1] <- 1
#' true_graph_neg[1, 3] <- true_graph_neg[3, 1] <- 1
#'
#' # visualize the true covariance structures
#' (gg_adjMat(true_graph_neg) +
#'     ggplot2::ggtitle("True graph for individuals with negative Z"))
#' (gg_adjMat(true_graph_pos, color1 = "steelblue") +
#'     ggplot2::ggtitle("True graph for individuals with positive Z"))
#'
#' # generate the covariance matrices as a function of Z
#' sigma_mats_neg <- lapply(Z_neg, function(z) z * true_graph_neg + diag(p + 1))
#' sigma_mats_pos <- lapply(Z_pos, function(z) z * true_graph_pos + diag(p + 1))
#' sigma_mats <- c(sigma_mats_neg, sigma_mats_pos)
#'
#' # generate the data using the covariance matrices
#' data_mat <- t(sapply(sigma_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))
#'
#' # visualize the sample correlation
#' gg_adjMat(abs(cor(data_mat[1:(n / 2), ])) - diag(p + 1))
#' gg_adjMat(abs(cor(data_mat[(n / 2 + 1):n, ])) - diag(p + 1),
#'           color1 = "dodgerblue")
#'
#' # estimate the covariance structure
#' out <- covdepGE(data_mat, Z)
#'
#' # analyze results
#' gg_adjMat(out, 1)
#' gg_adjMat(out, 50, color1 = "tomato")
#' gg_adjMat(out, 54, color1 = "steelblue")
#' gg_adjMat(out, 100, color1 = "dodgerblue")
#'
#' gg_inclusionCurve(out, 1, 2)
#' gg_inclusionCurve(out, 1, 3, point_color = "dodgerblue")
## -----------------------------------------------------------------------------
gg_adjMat <- function(out, l = 1, prob_shade = T, color0 = "white",
                      color1 = "#500000", grid_color = "black", incl_probs = T,
                      prob_prec = 2, font_size = 3, font_color0 = "black",
                      font_color1 = "white"){

  # takes care of no visible bindings
  Var1 <- Var2 <- value <- NULL
  rm(list = c("Var1", "Var2", "value"))

  # run compatibility checks
  adjMat_checks(out, l, prob_shade, color0, color1, grid_color, incl_probs,
                prob_prec, font_size, font_color0, font_color1)

  if (prob_shade & is.list(out)){

    # out is return of covdepGE; color the matrix using a gradient corresponding
    # to inclusion probabilities

    # get the posterior inclusion probability matrix for individual l
    probs <- out$inclusion_probs[[l]]

    # melt to long form
    long_probs <- reshape2::melt(probs)

    # if probabilities are to be displayed, get the adjacency matrix for
    # individual l and add the long version to long_probs for font coloring
    if (incl_probs){
      graph <- out$graphs[[l]]
      long_graph <- reshape2::melt(graph)
      long_graph$value <- as.factor(long_graph$value)
      long_probs$graph <- long_graph$value
    }

    vis <- (ggplot2::ggplot(long_probs, ggplot2::aes(Var1, Var2, fill = value)) +
              ggplot2::geom_tile(color = grid_color) +
              ggplot2::scale_fill_gradient(low = color0, high = color1,
                                           limits = c(0, 1)) +
              ggplot2::theme_classic() + ggplot2::xlab("") + ggplot2::ylab("") +
              ggplot2::ggtitle(paste("Posterior inclusion probabilties for individual",
                                     l)) +
              ggplot2::theme(legend.title = ggplot2::element_blank(),
                             plot.title = ggplot2::element_text(hjust = 0.5)) +
              ggplot2::scale_x_continuous(breaks = 1:nrow(graph)) +
              ggplot2::scale_y_continuous(breaks = 1:nrow(graph)))

    # if probabilities are to be displayed, add them
    if (incl_probs){
      vis <- (vis + ggplot2::geom_text(ggplot2::aes(label = round(value, prob_prec),
                                                    color = graph), show.legend = F,
                                       size = font_size) +
                ggplot2::scale_color_manual(values = c(font_color0, font_color1)))

    }


  }else if (is.list(out)){

    # out is return of covdepGE; color the matrix using a binary pallete

    # get the adjacency matrix for individual l
    graph <- out$graphs[[l]]

    # melt to long form
    long_graph <- reshape2::melt(graph)

    # factor the edges - specifies a discrete scale to ggplot2
    long_graph$value <- as.factor(long_graph$value)

    # if probabilities are to be displayed, get the posterior inclusion
    # probabilities for individual l and add them to the long_graph
    if (incl_probs){
      probs <- out$inclusion_probs[[l]]
      long_probs <- reshape2::melt(probs)
      long_graph$probs <- long_probs$value
    }

    vis <- (ggplot2::ggplot(long_graph, ggplot2::aes(Var1, Var2, fill = value)) +
              ggplot2::geom_tile(color = grid_color) +
              ggplot2::scale_fill_manual(values = c(color0, color1)) +
              ggplot2::theme_classic() + ggplot2::xlab("") + ggplot2::ylab("") +
              ggplot2::ggtitle(paste("Adjacency matrix for individual", l)) +
              ggplot2::theme(legend.title = ggplot2::element_blank(),
                             plot.title = ggplot2::element_text(hjust = 0.5)) +
              ggplot2::scale_x_continuous(breaks = 1:nrow(graph)) +
              ggplot2::scale_y_continuous(breaks = 1:nrow(graph)))


    # if probabilities are to be displayed, add them to the graph
    if (incl_probs){
      vis <- (vis + ggplot2::geom_text(ggplot2::aes(label = round(probs, prob_prec),
                                                    color = value), show.legend = F,
                                       size = font_size) +
                ggplot2::scale_color_manual(values = c(font_color0, font_color1)))
    }


  }else{

    # out is an adjacency matrix

    # save col and row names and remove them
    colnames <- rownames <- NULL
    colnames <- colnames(out)
    rownames <- rownames(out)
    out <- unname(out)

    # melt to long form
    long_graph <- reshape2::melt(out)

    # if the adjacency matrix has greater than 2 values, color using a gradient
    if (length(unique(as.vector(out))) > 2){

      # if probabilities are to be displayed, add an indicator to edges for the
      # squares that will be darker
      if (incl_probs){
        long_graph$graph <- ifelse(long_graph$value >
                                     mean(setdiff(long_graph$value, 0)), 1, 0)
        long_graph$graph <- as.factor(long_graph$graph)
      }

      vis <- (ggplot2::ggplot(long_graph, ggplot2::aes(Var1, Var2, fill = value)) +
                ggplot2::geom_tile(color = grid_color) +
                ggplot2::scale_fill_gradient(low = color0, high = color1) +
                ggplot2::theme_classic() + ggplot2::xlab("") + ggplot2::ylab("") +
                ggplot2::ggtitle("") +
                ggplot2::theme(legend.title = ggplot2::element_blank(),
                               plot.title = ggplot2::element_text(hjust = 0.5)) +
                ggplot2::scale_x_continuous(breaks = 1:nrow(out)) +
                ggplot2::scale_y_continuous(breaks = 1:nrow(out)))

      # if probabilities are to be displayed, add them
      if (incl_probs){
        vis <- (vis + ggplot2::geom_text(ggplot2::aes(label = round(value, prob_prec),
                                                      color = graph),
                                         show.legend = F, size = font_size) +
                  ggplot2::scale_color_manual(values = c(font_color0, font_color1)))
      }
    }else{

      # out has 2 or less unique values; color using a binary pallete

      # factor the edges - specifies a discrete scale to ggplot2
      long_graph$value <- as.factor(long_graph$value)

      vis <- (ggplot2::ggplot(long_graph, ggplot2::aes(Var1, Var2, fill = value)) +
                ggplot2::geom_tile(color = grid_color) +
                ggplot2::scale_fill_manual(values = c(color0, color1)) +
                ggplot2::theme_classic() + ggplot2::xlab("") + ggplot2::ylab("") +
                ggplot2::ggtitle("") +
                ggplot2::theme(legend.title = ggplot2::element_blank(),
                               plot.title = ggplot2::element_text(hjust = 0.5)) +
                ggplot2::scale_x_continuous(breaks = 1:nrow(out)) +
                ggplot2::scale_y_continuous(breaks = 1:nrow(out)))
    }

    # check to see if there are col names or row names to be added
    if (!is.null(colnames)){
      vis <- suppressMessages(vis + ggplot2::scale_x_continuous(breaks = 1:nrow(out),
                                                                labels = colnames))
    }
    if (!is.null(rownames)){
      vis <- suppressMessages(vis + ggplot2::scale_y_continuous(breaks = 1:nrow(out),
                                                                labels = rownames))
    }
  }

  return(vis)
}

## -----------------------------------------------------------------------------
#' @title gg_inclusionCurve
#' @export
## -----------------------------------------------------------------------------
## -----------------------------DESCRIPTION-------------------------------------
#' @description Create a visualization of the probabilities of inclusion of an
#' edge between two variables across all \eqn{n} individuals
## -----------------------------ARGUMENTS---------------------------------------
#' @param out object of class `covdepGE`; return of `covdepGE` function
#'
#' @param col_idx1 integer in \eqn{{1, 2, ..., p + 1}}; column index of the
#' first variable
#'
#' @param col_idx2 integer in \eqn{{1, 2, ..., p + 1}}; column index of the
#' second variable
#'
#' @param line_type linetype; `ggplot2` line type to interpolate the
#' probabilities. `"solid"` by default
#'
#' @param line_size positive numeric; thickness of the interpolating line. `0.5`
#' by default
#'
#' @param line_color color; color of interpolating line. `"black"` by default
#'
#' @param point_shape shape; shape of the points denoting individual-specific
#' inclusion probabilities; `21` by default
#'
#' @param point_size positive numeric; size of probability points. `1.5` by
#' default
#'
#' @param point_color color; color of probability points. `"#500000"` by default
#'
#' @param point_fill color; fill of probability points. Only applies to select
#' shapes. `"white"` by default
#'
#' @param sort logical; if `T`, applies a piori sorting algorithm to re-order
#' subject indices according to weights such for \eqn{j} in \eqn{1,...,n - 1},
#' subject \eqn{j + 1} is the most similar in terms of their extraneous
#' covariate to subject \eqn{j}. Helps visualize the continuity with which the
#' edge probabilities are modeled with respect to the the extraneous covariates
## -----------------------------RETURNS-----------------------------------------
#' @return Returns `ggplot2` visualization of inclusion probability curve
## -----------------------------EXAMPLES----------------------------------------
#' @examples
#' set.seed(1)
#' n <- 100
#' p <- 4
#'
#' # generate the extraneous covariate
#' Z_neg <- sort(runif(n / 2) * -1)
#' Z_pos <- sort(runif(n / 2))
#' Z <- c(Z_neg, Z_pos)
#' summary(Z)
#'
#' # create true covariance structure for 2 groups: positive Z and negative Z
#' true_graph_pos <- true_graph_neg <- matrix(0, p + 1, p + 1)
#' true_graph_pos[1, 2] <- true_graph_pos[2, 1] <- 1
#' true_graph_neg[1, 3] <- true_graph_neg[3, 1] <- 1
#'
#' # visualize the true covariance structures
#' (gg_adjMat(true_graph_neg) +
#'     ggplot2::ggtitle("True graph for individuals with negative Z"))
#' (gg_adjMat(true_graph_pos, color1 = "steelblue") +
#'     ggplot2::ggtitle("True graph for individuals with positive Z"))
#'
#' # generate the covariance matrices as a function of Z
#' sigma_mats_neg <- lapply(Z_neg, function(z) z * true_graph_neg + diag(p + 1))
#' sigma_mats_pos <- lapply(Z_pos, function(z) z * true_graph_pos + diag(p + 1))
#' sigma_mats <- c(sigma_mats_neg, sigma_mats_pos)
#'
#' # generate the data using the covariance matrices
#' data_mat <- t(sapply(sigma_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))
#'
#' # visualize the sample correlation
#' gg_adjMat(abs(cor(data_mat[1:(n / 2), ])) - diag(p + 1))
#' gg_adjMat(abs(cor(data_mat[(n / 2 + 1):n, ])) - diag(p + 1),
#'           color1 = "dodgerblue")
#'
#' # estimate the covariance structure
#' out <- covdepGE(data_mat, Z)
#'
#' # analyze results
#' gg_adjMat(out, 1)
#' gg_adjMat(out, 50, color1 = "tomato")
#' gg_adjMat(out, 54, color1 = "steelblue")
#' gg_adjMat(out, 100, color1 = "dodgerblue")
#'
#' gg_inclusionCurve(out, 1, 2)
#' gg_inclusionCurve(out, 1, 3, point_color = "dodgerblue")
## -----------------------------------------------------------------------------
gg_inclusionCurve <- function(out, col_idx1, col_idx2, line_type = "solid",
                              line_size = 0.5, line_color = "black",
                              point_shape = 21, point_size = 1.5,
                              point_color = "#500000", point_fill = "white",
                              sort = F){

  # takes care of "no visible bindings" note
  idx <- prob <- NULL
  rm(list = c("idx", "prob"))

  # run compatibility checks
  inclusionCurve_checks(out, col_idx1, col_idx2, line_type, line_size,
                        line_color, point_shape, point_size, point_color,
                        point_fill, sort)

  # get the probabilities for each individual of an edge between the variables
  # corresponding to col_idx1 and col_idx2
  prob1_2 <- unlist(lapply(out$inclusion_probs, function(x) x[col_idx1, col_idx2]))

  if (sort){

    # indv_sorted will be a vector of subject indices sorted such that the
    # individuals corresponding to neighboring indices are most similar to each
    # other (have the highest weighting)
    indv_sorted <- c(1, rep(NA, nrow(out$weights) - 1))

    for (l in 2:nrow(out$weights)){
      # get the index of the last individual that was added to ind_sorted
      last_indv_idx <- indv_sorted[l - 1]

      # find the individual whose weight is largest with respect to the last
      # individual, excluding those who are already in indv_sorted
      next_indv_idx <- setdiff(order(out$weights[ , last_indv_idx], decreasing = T),
                               indv_sorted)[1]

      # add this individual to the sorted list
      indv_sorted[l] <- next_indv_idx
    }

    # sort prob1_2 according to indv_sorted
    prob1_2 <- prob1_2[indv_sorted]

    # find the index of the largest discontinuity
    disc_lg_ind <- which.max(abs(diff(prob1_2)))

    # if the largest discontinuity is larger than the discontinuity between
    # the current first point and the right hand point of discontinuity, move
    # the points following the discontinuity to the front
    disc_left1 <- abs(diff(prob1_2[c(1, disc_lg_ind + 1)]))
    if (disc_left1 < max(abs(diff(prob1_2)))){

      # move all probabilities following the discontinuity to the front
      prob1_2 <- c(prob1_2[(disc_lg_ind + 1):length(prob1_2)],
                   prob1_2[1:disc_lg_ind])
    }
  }

  prob1_2 <- data.frame(idx = 1:length(prob1_2), prob = prob1_2)

  vis <- (ggplot2::ggplot(prob1_2, ggplot2::aes(idx, prob)) +
            ggplot2::geom_line(linetype = line_type, size = line_size,
                               color = line_color) +
            ggplot2::geom_point(shape = point_shape, size = point_size,
                                color = point_color, fill = point_fill) +
            ggplot2::xlab("Subject index") +
            ggplot2::ylab("Posterior Inclusion Probability") +
            ggplot2::ggtitle(latex2exp::TeX(
              paste("Inclusion probability of an edge between $x_",
                    col_idx1, "$ and $x_", col_idx2, "$"))) +
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
#' @description Given the return value of covdepGE function, create a list of
#'  visualizations of the adjacency matrix for each of the unique graphs
## -----------------------------ARGUMENTS---------------------------------------
#' @param out object of class covdepGE; return of `covdepGE` function
#'
#' @param graph_colors: vector of length g; g is the number of unique graphs
#' from out. The v-th element vector is the color for the v-th unique graph
#'
#' @param seed numeric; when colors is NULL, the RNG seed for selecting the
#' color for each graph. 1 by default.
## -----------------------------RETURNS-----------------------------------------
#' @return Returns list of `ggplot2` visualizations of unique graph adjacency
#' matrices
## -----------------------------EXAMPLES----------------------------------------
#' @examples
#' set.seed(1)
#' n <- 100
#' p <- 4
#'
#' # generate the extraneous covariate
#' Z_neg <- sort(runif(n / 2) * -1)
#' Z_pos <- sort(runif(n / 2))
#' Z <- c(Z_neg, Z_pos)
#' summary(Z)
#'
#' # create true covariance structure for 2 groups: positive Z and negative Z
#' true_graph_pos <- true_graph_neg <- matrix(0, p + 1, p + 1)
#' true_graph_pos[1, 2] <- true_graph_pos[2, 1] <- 1
#' true_graph_neg[1, 3] <- true_graph_neg[3, 1] <- 1
#'
#' # visualize the true covariance structures
#' (gg_adjMat(true_graph_neg) +
#'     ggplot2::ggtitle("True graph for individuals with negative Z"))
#' (gg_adjMat(true_graph_pos, color1 = "steelblue") +
#'     ggplot2::ggtitle("True graph for individuals with positive Z"))
#'
#' # generate the covariance matrices as a function of Z
#' sigma_mats_neg <- lapply(Z_neg, function(z) z * true_graph_neg + diag(p + 1))
#' sigma_mats_pos <- lapply(Z_pos, function(z) z * true_graph_pos + diag(p + 1))
#' sigma_mats <- c(sigma_mats_neg, sigma_mats_pos)
#'
#' # generate the data using the covariance matrices
#' data_mat <- t(sapply(sigma_mats, MASS::mvrnorm, n = 1, mu = rep(0, p + 1)))
#'
#' # visualize the sample correlation
#' gg_adjMat(abs(cor(data_mat[1:(n / 2), ])) - diag(p + 1))
#' gg_adjMat(abs(cor(data_mat[(n / 2 + 1):n, ])) - diag(p + 1),
#'           color1 = "dodgerblue")
#'
#' # estimate the covariance structure
#' out <- covdepGE(data_mat, Z)
#'
#' # analyze results
#' plot(out)
#'
#' gg_inclusionCurve(out, 1, 2)
#' gg_inclusionCurve(out, 1, 3, point_color = "dodgerblue")
## -----------------------------------------------------------------------------
plot.covdepGE <- function(out, graph_colors = NULL, seed = 1){

  # compatibility checks
  plot_checks(out, seed, graph_colors)

  # if no colors have been provided, select some randomly
  if(is.null(graph_colors)){
    set.seed(seed)
    graph_colors <- colors()[sample(1:length(colors()), length(out$unique_graphs))]
  }

  # get the summary of individuals corresponding to each graph
  indv_sum <- sapply(out$unique_graphs, `[[`, "individuals_summary")

  # get the unique graphs
  unique_graphs <- lapply(out$unique_graphs, `[[`,"graph")

  # create a visualization for each graph and store it in a list
  graph_viz <- lapply(1:length(unique_graphs), function(gr_idx) gg_adjMat(
    unique_graphs[[gr_idx]], color1 = graph_colors[gr_idx]) + ggplot2::ggtitle(
      paste("Individuals", indv_sum[gr_idx])))

  return(graph_viz)
}
