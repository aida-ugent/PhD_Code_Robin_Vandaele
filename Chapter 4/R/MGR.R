#' Reconstruct graph topology based on pointwise edge-branch assignment
#'
#' @param D pairwise distance object of data frame
#' @param r distance parameter
#' 
#' @return Reconstructed graph

MGR <- function(D, r){
  vr <- VRgraph(D, 4 * r / 3) 
  mgl <- integer(attr(D, "Size"))
  for(i in seq(attr(D, "Size"))){
    S <- as.matrix(D)[as.matrix(expand.grid(i, seq(attr(D, "Size"))))]
    deg <- components(induced_subgraph(vr, which((S < 5 * r / 3) & (S >= r))))$no
    if(deg == 2) mgl[i] <- 1
  }
  B <- as.matrix(D)[as.matrix(expand.grid(which(mgl == 0), seq(attr(D, "Size"))))]
  B <- which(apply(matrix(as.matrix(D)[as.matrix(
    expand.grid(which(mgl == 0), seq(attr(D, "Size"))))], 
    nrow = length(B)), 2, min) < 2 * r)
  mgl[B] <- 0
  vr <- VRgraph(D, 2 * r)
  IV <-  which(mgl == 0)
  IE <-  which(mgl == 1)
  V <- components(induced_subgraph(vr, IV))
  E <- components(induced_subgraph(vr, IE))
  adj <- rep(list(integer(0)), V$no)
  if(E$no > 0){
    for(i in seq(V$no)){
      I <- IV[which(V$membership == i)]
      for(j in i:V$no){
        J <- IV[which(V$membership == j)]
        if(i < j){
          for(k in seq(E$no)){
            K <- IE[which(E$membership == k)]
            dik <- min(as.matrix(D)[as.matrix(expand.grid(I, K))])
            djk <- min(as.matrix(D)[as.matrix(expand.grid(J, K))])
            if(max(dik, djk) < 2 * r){
              adj[[i]] <- c(adj[[i]], j)
              break
            }
          }
        }
      }
    }
    for(i in seq(E$no)){
      n <- 0
      I <- IE[which(E$membership == i)]
      for(j in seq(V$no)){
        J <- IV[which(V$membership == j)]
        if(min(as.matrix(D)[as.matrix(expand.grid(I, J))]) < 2 * 3){
          n <- n + 1
          if(n == 2) break
          v <- j
        }
      }
      if(n == 1) adj[[v]] <- c(adj[[v]], v)
      else if(n == 0) adj[[length(adj) + 1]] <- length(adj) + 1
    }
  }
  mgr <- as.undirected(graph_from_adj_list(adj))
  mgr$assignment <- factor(mgl)
  return(mgr)
}