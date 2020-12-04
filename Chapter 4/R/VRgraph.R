#' Construct Vietoris-Rips graph from a point cloud data set
#'
#' @param D   pairwise distance object to be used for constructing the complex
#' @param eps upper bound on edge lengths
#' 
#' @return    the Vietoris-Rips graph defined by df and distance parameter eps

VRgraph <- function(D, eps){
  adj <- which(as.matrix(D) < eps, arr.ind = TRUE)
  adj <- adj[which(adj[,1] > adj[,2]),]
  vr <- setNames(rep(list(integer(0)), tail(adj, 1)[2]), 1:tail(adj, 1)[2])
  vr[unique(adj[,2])] <- split(adj[,1], adj[,2])
  vr <- as.undirected(graph_from_adj_list(vr))
  weights <- numeric(length(E(vr)))
  for(e in seq(E(vr))){
    edge <- get.edges(vr, e)
    weights[e] <- D[attr(D, "Size") * (edge[1] - 1) - edge[1] * (edge[1] - 1) / 2 +
                           edge[2] - edge[1]]
  }
  E(vr)$weight <- weights
  return(vr)
}