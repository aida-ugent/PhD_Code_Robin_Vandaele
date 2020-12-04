#' Prune all leaves in a graph either a given amount of times or until there are no more than a specified number of leaves
#'
#' @param G      graph to be pruned
#' @param times  the number of times all leaves must be pruned, standard set to 1
#' @param leaves (maximum) number of leaves to be kept, 'times' is ignored if provided
#'
#' @return the pruned graph

hard_prune <- function(G, times=1, leaves=NA){
  if(times < 0) stop("'times' must be at least 0")
  if(times == 0) return(G)
  leaves <- max(as.integer(leaves), 0, na.rm = TRUE)
  if(leaves >= 2) times <- Inf
  to_prune <- which(degree(G) == 1)
  it <- 0
  while(length(to_prune) > leaves){
    it <- it + 1
    G <- delete_vertices(G, to_prune)
    to_prune <- which(degree(G) == 1)
    if(it >= times) break
  }
  return(G)
}
