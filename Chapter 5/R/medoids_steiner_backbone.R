#' Find a backbone in a graph G through building a Steiner tree in G through k-medoids
#'
#' @param G     graph from which backbone needs to be extracted
#' @param k     the number of terminal nodes that must be selected for each component of G
#'
#' @return the pruned graph

medoids_steiner_backbone <- function(G, k){
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
  terminals <- list()
  for(i in 1:C$no){
    I <- which(C$membership == i)
    if(k[i] > 0){
      if(length(I) == k[i]) terminals[[i]] <- V(G)$name[I]
      else  terminals[[i]] <- V(G)$name[I][cluster::pam(as.dist(distances(G, I, I)), k[i])$id.med]
    }
  }
  new <- Sys.time() - old
  print(paste("Terminal nodes obtained in", round(new, 3), attr(new, "units")))
  old <- Sys.time()
  MSB <- do.call(union, lapply(terminals, function(t) SteinerNet::steinertree("KB", terminals=t, graph=G)[[2]]))
  MSB$terminals <- unlist(terminals)
  new <- Sys.time() - old
  print(paste("Steiner tree approximated in", round(new, 3), attr(new, "units")))
  return(MSB)
}
