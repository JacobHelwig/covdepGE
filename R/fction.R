#' Soft-thresholding function
#'
#' @param x A scalar to be thresholded
#' @param lambda A parameter for thresholding a non-negative scalar
#'
#' @return Returns x soft-thresholded by lambda
#' @export
#'
#' @examples soft(3, 5)
soft <- function(x, lambda){
  return(sign(x) * max(abs(x) - lambda), 0)
}
