#' Construct Vietoris-Rips graph from a point cloud data set by connecting all
#' nodes with a certain given distance eps
#'
#' @param X   either a data.frame/matrix containing observations,
#'            or a pairwise distance object to be used for constructing the complex
#' @param eps (strict) upper bound on edge lengths
#'
#' @return the Vietoris-Rips graph defined by X and distance parameter eps

VRgraph <- function(X, eps){
  if(is.data.frame(X) | is.matrix(X)){
    adj <- dplyr::bind_rows(lapply(1:(nrow(X) - 1), function(n){
      d <- pdist::pdist(X[n,], X[(n + 1):nrow(X),])@dist
      to <- which(d < eps)
      return(data.frame(from=rep(n, length(to)), to=((n + 1):nrow(X))[to], weight=d[to]))
    }))
    G <- graph.data.frame(adj, directed=FALSE)
    V(G)$name <- rownames(X)[as.integer(V(G)$name)]
  }
  else if(class(X) == "dist"){
    adj <- which(as.matrix(X) < eps, arr.ind=TRUE)
    adj <- adj[adj[,1] < adj[,2],]
    G <- graph.data.frame(adj, directed=FALSE)
    V(G)$name <- attr(X, "Labels")[as.integer(V(G)$name)]
    E(G)$weight <- as.matrix(X)[adj]
  }
  else stop("X must be of either of class data.frame/matrix, or of class dist")
  return(G)
}
