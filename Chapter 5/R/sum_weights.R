#' Compute softmax weights for a graph using a predefined function
#'
#' @param G  Graph containing the edges for which the softmax values have to be computed
#' @param f  vector of function values of length |V(G)|
#' 
#' @return a vector containing the sum weight for each edge in G, ordered as E(G) 

sum_weights <- function(G, f){
  W <- get.edges(G, E(G))
  W <- f[W[,1]] + f[W[,2]]
  return(W)
}