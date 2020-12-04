#' Compute the weighted distances between nodes if there unweighted distance is less than a given integer k
#' This code is more memory efficient than computing the full paiwise distance matrix and discarding unnecessary nodes,
#' but may run slower due to the optimized implementation of igraphs distances function in C
#'
#' @param G             the graph for which to compute the
#' @param k             maximum allowd hops between pairs of nodes for which the weighted distance is to be computed
#' @param constrain_mem whether to compute the entire distance matrix on G as an intermediate step to speed up the (currently inefficient) computation,
#'                      at the cost of more memory usage, standard is FALSE
#'
#' @return a sparse matrix storing the weighted distance D_uv between nodes u and v if the unweighted distance
#'         between u and v is less than k

bounded_hop_pairG <- function(G, k=2, constrain_mem=FALSE){
  to <- ego(G, order=k)
  i <- as.integer(unlist(sapply(1:length(to), function(n) rep(n, length(to[[n]])))))
  to <- unlist(to)
  if(!constrain_mem) return(spam(x=list(i=i, j=to, values=distances(G)[cbind(i, to)])))
  return(spam(x=list(i=i, j=to, values=unlist(lapply(1:length(V(G)), function(v) distances(G, v=v, to=to[[v]]))))))
}
