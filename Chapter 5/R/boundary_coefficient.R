#' Compute the boundary coefficients for a weighted graph
#'
#' @param G             the graph for which to compute the boundary measure of each node
#' @param pairD         pairwise distance object of G, computed if missing
#' @param constrain_mem whether to compute the entire distance matrix on G as an intermediate step to speed up the (currently inefficient) computation,
#'                      at the cost of more memory usage, standard is FALSE
#'
#' @return a vector of length |V(G)| containing the boundary coefficients

boundary_coefficient <- function(G, pairD=NULL, constrain_mem=FALSE){
  if(is.null(pairD)) pairD <- bounded_hop_pairG(G, constrain_mem=constrain_mem)
  hop1 <- as.matrix(summary(as_adj(G)))[,c(1, 2)]
  A <- spam(x=list(indices=hop1, values=pairD[hop1]))
  oneOvA <- 1 / A
  BC <- (apply.spam(A, 1, sum) * apply.spam(oneOvA, 1, sum) - diag(oneOvA %*% pairD^2%*% oneOvA) / 2) / (degree(G)^2)
  return(BC)
}
