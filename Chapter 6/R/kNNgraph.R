#' Construct a undirected kNN graph from a distance object by connecting each node
#' to its k closest neighbors
#'
#' @param X either a data.frame/matrix containing observations,
#'          or a pairwise distance object to be used for constructing the complex
#' @param k number of neighbors in the kNN graph, standard is 10
#'
#' @return  undirected k-nearest neighbor igraph object

kNNgraph <- function(X, k=10){
  if(is.data.frame(X) | is.matrix(X)){
    if(k > nrow(X) - 1) stop("k must be at least equal to nrow(X) - 1")
    kNN <- FNN::get.knn(X, k)
    adj <- cbind(from=rep(1:nrow(X), each=k), to=as.integer(t(kNN[[1]])), weight=as.numeric(t(kNN[[2]])))
    I <- adj[,1] > adj[,2]
    adj[I, c(1, 2)] <- adj[I, c(2, 1)]
    adj <- adj[!duplicated(adj[,c(1, 2)]),]
    G <- graph.data.frame(adj, directed=FALSE)
    V(G)$name <- rownames(X)[as.integer(V(G)$name)]
  }
  else if(class(X) == "dist"){
    if(k > attr(X, "Size") - 1) stop("k must be at least equal to nrow(X) - 1")
    adj <- cbind(from=rep(1:attr(X, "Size"), each=k),
                  to=as.integer(apply(as.matrix(X), 1, function(r) order(r)[2:(k + 1)])))
    I <- adj[,1] > adj[,2]
    adj[I, c(1, 2)] <- adj[I, c(2, 1)]
    adj <- adj[!duplicated(adj[,c(1, 2)]),]
    G <- graph.data.frame(adj, directed=FALSE)
    V(G)$name <- attr(X, "Labels")[as.integer(V(G)$name)]
    E(G)$weight <- as.matrix(X)[adj]
  }
  else stop("X must be of either of class data.frame/matrix, or of class dist")
  return(G)
}
