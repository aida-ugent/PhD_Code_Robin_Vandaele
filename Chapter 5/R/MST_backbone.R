#' Find a backbone in a graph G through furthes point sampling and connecting them by paths
#'
#' @param G     graph from which backbone needs to be extracted
#' @param k     the number of leaves must be included for each component of G
#' @param vcost the type of vertex cost associated to the MST to maximize (must be "degree" or "betweenness")
#'
#' @return the pruned graph

MST_backbone <- function(G, k, vcost="betweenness"){
  if(!vcost %in% c("degree", "betweenness")) stop("vcost must be 'degree', or 'betweenness', standard is 'betweenness'")
  old <- Sys.time()
  MST <- mst(G)
  new <- Sys.time() - old
  print(paste("MST obtained in", round(new, 3), attr(new, "units")))
  if(vcost=="betweeness") MSTB <- prune(MST, leaves=k, vcost=betweenness(MST, weights=NA), method="search")
  else MSTB <- prune(MST, leaves=k, vcost=degree(MST), method="search")
  new <- Sys.time() - old
  print(paste("Backbone obtained in", round(new, 3), attr(new, "units")))
  return(MSTB$B)
}
