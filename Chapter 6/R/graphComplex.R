#' Convert a igraph object to a simplicial complex
#'
#' @param G the input graph
#'
#' @return a simplicial complex S

graphComplex <- function(G){
  S <- c(lapply(V(G), function(v) as.integer(v)), split(t(ends(G, E(G), names=FALSE)), rep(1:length(E(G)), each=2)))
  return(S)
}
