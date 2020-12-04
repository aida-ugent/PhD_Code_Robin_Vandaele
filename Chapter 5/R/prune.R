#' Prune a graph a given number of times by choosing the number of times/leaves and used method
#'
#' @param G      graph to be pruned
#' @param times  times the graph must be pruned, only used if 'method' is 'hard' or 'soft', standard set to 1
#'               may also be a vector with size the number of components of 'the (constructed) graph 'G'
#' @param leaves (maximum) number of leaves to be included in the subtree, standard set to NA, 'times' is ignored if provided
#'               may also be a vector with size the number of components of 'the (constructed) graph 'G'
#' @param method used method to perform the pruning, standard set to 'hard':
#'                    - 'hard': iteratively prune all leaves in the graph
#'                    - 'soft': iteratively prune the leaves representing the fewest nodes
#'                    - 'search': search for a set of leaves resulting in a high cost tree optimal if G itself is a tree)
#' @param weight optional reweighting vector of the edges of G (used only with search prune)
#' @param vcost  optional vertex cost for which the sum is to be optimized (used only with search prune, ignores weight if provided)
#'
#' @return A list containing the following items
#'                - B: the pruned graph derived from 'G'
#'                - membership: the membership of the 'X' points/nodes according to the components of the (constructed) graph
#'                - (if method is 'search') cost: a data.frame corresponding to the obtained cost according to the number of leaves and component
#'                - (if method is 'search') full_cost: a vector of the full cost of each component
#'                - (if method is 'search') includedE: a list of lists (one for each connected component) containing the added edges of 'G' at each iteration

prune <- function(G, times=1, leaves=NA, method="hard", weight=NULL, vcost=NULL, preprune=TRUE){
  if(!method %in% c("hard", "soft", "search")) stop("method must be either 'hard', 'soft', or 'search'")
  if(!is.null(vcost) & is.null(names(vcost))) names(vcost) <- names(G)
  C <- components(G)
  if(length(leaves) == 1) leaves <- rep(leaves, C$no)
  else if(length(leaves) != C$no) stop("Length of 'leaves' must be 1 or the number of components of G")
  if(length(times) == 1) times <- rep(times, C$no)
  else if(length(times) != C$no) stop("Length of 'times' must be 1 or the number of components of G")
  nodes_to_add <- character()
  edges_included <- integer()
  if(method=="search"){
    cost <- data.frame(integer(0), numeric(0), integer(0))
    full_cost <- numeric(C$no)
    if(!is.null(vcost)) includedV <- rep(list(NULL), C$no)
    else includedE <- rep(list(NULL), C$no)
  }
  for(c in seq(C$no)){
    H <- induced.subgraph(G, which(C$membership == c))
    if(method == "hard") H <- list(B=hard_prune(H, times=times[c], leaves=leaves[c]))
    else if(method == "soft") H <- list(B=soft_prune(H, times=times[c], leaves=leaves[c]))
    else H <- search_prune(H, leaves=leaves[c], weight=weight[get.edge.ids(H, t(ends(H, E(H))))], vcost=vcost[names(V(H))], preprune=preprune)
    if(!is.null(H$B)){
      nodes_to_add <- c(nodes_to_add, V(H$B)$name[which(degree(H$B) == 0)])
      edges_included <- c(edges_included, get.edge.ids(G, t(ends(H$B, E(H$B)))))
      if(!is.null(H$cost)) cost <- rbind(cost, cbind(H$cost, c))
      if(!is.null(H$full_cost)) full_cost[c] <- H$full_cost
      if(!is.null(H$includedV)) includedV[[c]] <- H$includedV
      if(!is.null(H$includedE)) includedE[[c]] <- H$includedE
    }
  }
  pruned <- list(B = subgraph.edges(G, edges_included))
  pruned$membership <- C$membership
  if(length(nodes_to_add) > 0){
    pruned$B <- add_vertices(pruned$B, length(nodes_to_add))
    V(pruned$B)$name[seq(length(V(pruned$B)) - length(nodes_to_add) + 1, length(V(pruned$B)))] <- nodes_to_add
  }
  if(method=="search"){
    colnames(cost)[3] <- "component"
    cost[,3] <- factor(cost[,3])
    cost <- cost[order(cost[,1]),]
    rownames(cost) <- 1:nrow(cost)
    pruned$cost <- cost
    if(!is.null(full_cost)) pruned$full_cost <- full_cost
    if(!is.null(vcost)) pruned$includedV <- includedV
    else pruned$includedE <- includedE
  }
  return(pruned)
}
