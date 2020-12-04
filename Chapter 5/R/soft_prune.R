#' Prune leaves in a graph either a given amount of times or until there are no more than a specified number of leaves
#' At each iteration, the leaves representing the fewest amount of nodes are pruned
#' Each node is set to represent only itself at the start of the iterations
#' When leaves are pruned, they are consequently represented by their sole neighbor
#'
#' @param G      graph to be pruned
#' @param times  the number of times all leaves must be pruned, standard set to 1
#' @param leaves (maximum) number of leaves to be kept, 'times' is ignored if provided
#' 
#' @return the pruned graph

soft_prune <- function(G, times=1, leaves=NA){
  if(times < 0) stop("'times' must be at least 0")
  if(times == 0) return(G)
  if(times == 1 & is.na(leaves)) return(hard_prune(G))
  leaves <- max(as.integer(leaves), 0, na.rm = TRUE)
  if(leaves >= 2) times <- Inf 
  to_prune <- which(degree(G) == 1)
  no_leaves <- length(to_prune)
  NoR <- rep(1, length(V(G)))
  names(NoR) <- V(G)$name
  it <- 0
  while(no_leaves > leaves){
    it <- it + 1
    for(u in to_prune){
      v <- neighbors(G, u)$name
      NoR[v] <- NoR[v] + NoR[V(G)[u]$name] 
    }
    G <- delete_vertices(G, to_prune)
    if(it >= times) break
    to_prune <- which(degree(G) == 1)
    no_leaves <- length(to_prune)
    if(length(to_prune) > 0)  to_prune <- to_prune[NoR[V(G)$name[to_prune]] == min(NoR[V(G)$name[to_prune]])]
  }
  return(G)
}