#' Compute the degree of centrality for each node in a graph
#'
#' @param G the input graph
#'
#' @return a numeric vector containing the degrees of centrality

degreeOfCentrality <- function(G){
  vertexEccentricities <- apply(distances(G), 1, max)
  m <- min(vertexEccentricities )
  M <- max(vertexEccentricities )
  degreeOfCentrality <- (M - vertexEccentricities) / (M - m)
  return(degreeOfCentrality)
}
