#' Find a backbone in a graph G through furthest point sampling and connecting them by shortest paths
#'
#' @param G graph from which backbone needs to be extracted
#' @param k the number of points must be sampled from each component of G
#'
#' @return the pruned graph

furthest_point_backbone <- function(G, k){
  if(sum(k) < 0) stop("number of points sampled must be at least 0")
  old <- Sys.time()
  if(sum(k) == 0){
    new <- Sys.time() - Old
    print(paste("Backbone obtained in", round(new, 3), attr(new, "units")))
    return(induced_subgraph(G, NULL))
  }
  C <- components(G)
  if(length(k) == 1) k <- rep(k, C$no)
  else if(length(k) != C$no) stop("length of k must be 1 or the number of components")
  VtoAdd <- character(0)
  EtoAdd <- E(G)[0]
  for(i in 1:C$no){
    I <- which(C$membership == i)
    if(length(I) == 1){
      if(k[i] > 0) VtoAdd <- cbind(VtoAdd, V(G)$name[I])
    }
    else if (k[i] == 1){
      VtoAdd <- cbind(VtoAdd, V(G)$name[I][which.min(apply(distances(G, I, I), 1, max))])
    }
    else if(k[i] > 1){
      D <- distances(G, I, I)
      u <- apply(D, 1, which.max)
      v <- which.max(D[cbind(1:length(I), u)])
      u <- u[v]
      l <- 2
      P <- shortest_paths(G, V(G)$name[I][u], V(G)$name[I][v], output="both")
      EtoAdd <- c(EtoAdd, P$epath[[1]])
      closestTreePoint <- character(length(I))
      while(l < k[i]){
        closestPathPoint <- names(P$vpath[[1]])[apply(D[,names(P$vpath[[1]])], 1, which.min)]
        if(l == 2) closestTreePoint <- closestPathPoint
        else{
          newClosest <- which(D[cbind(V(G)$name[I], closestPathPoint)] < D[cbind(V(G)$name[I], closestTreePoint)])
          closestTreePoint[newClosest] <- closestPathPoint[newClosest]
        }
        u <- which.max(D[cbind(V(G)$name[I], closestTreePoint)])
        v <- closestTreePoint[u]
        if(V(G)$name[I][u]==v) break
        P <- shortest_paths(G, V(G)$name[I][u], v, output="both")
        EtoAdd <- c(EtoAdd, P$epath[[1]])
        l <- l + 1
      }
    }
  }
  FPSB <- subgraph.edges(G, EtoAdd)
  FPSB <- add_vertices(FPSB, length(VtoAdd), name=VtoAdd)
  new <- Sys.time() - old
  print(paste("Backbone obtained in", round(new, 3), attr(new, "units")))
  return(FPSB)
}
