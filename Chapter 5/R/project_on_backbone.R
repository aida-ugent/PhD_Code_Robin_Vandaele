#' project a given graph G on a subgraph B
#'
#' @param G  igraph object
#' @param B  subgraph of G
#' @param dG pairwise distance matrix of G, computed if not provided
#'
#' @return  new graph GprojB, the projection of G on B

project_on_backbone <- function(G, B, dG=NULL){
  VtoAdd <- setdiff(names(V(G)), names(V(B)))
  if(is.null(dG)) dG <- distances(G, VtoAdd, names(V(B)))
  connectTo <- names(V(B))[apply(dG[VtoAdd, names(V(B))], 1, which.min)]
  EtoAdd <- unique(unlist(lapply(1:length(VtoAdd), function(n)
    shortest_paths(G, VtoAdd[n], connectTo[n], output="epath")$epath[[1]])))
  GprojB <- add.edges(add_vertices(B, length(VtoAdd), name=VtoAdd), t(ends(G, EtoAdd)))
  if(!is.null(E(B)$weight)) E(GprojB)$weight[(length(E(B)) + 1):(length(E(GprojB)))] <- E(G)$weight[EtoAdd]
  GprojB$addedV <- VtoAdd
  GprojB$connectTo <- connectTo
  return(GprojB)
}
