#' Compute index of center point in a (sub) data frame
#'
#' @param I  indices of points for which the center is to be computed
#' @param D  pairwise distance object of the full data frame
#' 
#' @return   the index of the center of df[I,]

center <- function(I, D) {
  subD <- as.matrix(expand.grid(I, I))
  subD <- matrix(as.matrix(D)[subD], nrow = length(I))
  c <- I[which.min(apply(subD, 1, max))]
  return(c)
}
