#' Get a new backbone with a different amount of leaves
#' This method only works when the backbone is returned using the 'search' method to prune
#'
#' @param backbone a backbone object to be converted to a new backbone object with a different amoubt of leaves
#' @param leaves   vector of integers, denoting the number of leaves to be included for each component treated separately in the pine
#' @param G        the original graph on top of the backbone, standard set to NULL
#' @param assign   whether to produce a branch assignment of the resulting new backbone, standard set to FALSE
#' @param prob     whether to assign nodes in G to paths according to a probability measure
#'                 if FALSE, then nodes are assigned to exactly one path, standard is TRUE
#' @param col      whether to provide a coloring scheme for the paths and nodes, standard set to TRUE
#' @param stdize   whether to standardize the cost of each component according to its full cost, standard set to FALSE
#'
#' @return a backbone object with the new number of leaves

get_new_leaves <- function(backbone, leaves, G=NULL, assign=FALSE, prob=TRUE, col=TRUE, stdize=FALSE){
  if(stdize) adjust <- rep(1, length(backbone$includedV))
  else adjust <- backbone$full_cost
  if(!is.null(backbone$includedV)){
    if(length(leaves)==1 & length(backbone$includedV) > 1){
      comp_rows <- lapply(1:length(backbone$includedV), function(c) which(backbone$cost[,"component"]==c))
      poss_leaves <- 1:leaves
      leave_rows <- lapply(poss_leaves, function(l) which(backbone$cost[,"leaves"]==(if(l == 1) 0 else l)))
      names(leave_rows) <- poss_leaves
      leavecombs <- get_leave_combs(length(backbone$includedV), leaves)
      cost <- apply(leavecombs, 1, function(r){
        sum(sapply(1:length(r), function(i){
          cr <- intersect(comp_rows[[i]], leave_rows[[as.character(r[i])]])
          if(length(cr) == 0) return(0)
          backbone$cost[cr, "cost"] * adjust[i]
        }))
      })
      leaves <- leavecombs[which.max(cost),]
      rm(comp_rows, leave_rows, leavecombs, cost)
    }
    else if(length(leaves) != length(backbone$includedV))
      stop("length of leaves must 1 or the number of components in the original graph")
    nodes <- character()
    for(c in seq(length(leaves))){
      if(is.na(leaves[c])) nodes <- c(nodes, unlist(backbone$includedV[[c]][min(length(backbone$includedV[[c]]), 2)]))
      else if(leaves[c] == 1) nodes <- c(nodes, unlist(backbone$includedV[[c]][1]))
      else if(leaves[c] > 1) nodes <- c(nodes, unlist(backbone$includedV[[c]][2:leaves[c]]))
    }
    backbone$B <- induced.subgraph(backbone$pine, nodes)
  }
  else{
    if(length(leaves)==1 & length(backbone$includedE) > 1){
      comp_rows <- lapply(1:length(backbone$includedE), function(c) which(backbone$cost[,"component"]==c))
      if(leaves %in% c(2, 3)) poss_leaves <- c(0, leaves)
      else poss_leaves <- c(0, 2:(leaves - 2), leaves)
      leave_rows <- lapply(poss_leaves, function(l) which(backbone$cost[,"leaves"]==l))
      names(leave_rows) <- poss_leaves
      leavecombs <- get_leave_combs(length(backbone$includedE), leaves)
      leavecombs <- leavecombs[!apply(leavecombs, 1, function(r) any(r==1)),]
      cost <- apply(leavecombs, 1, function(r){
        sum(sapply(1:length(r), function(i){
          cr <- intersect(comp_rows[[i]], leave_rows[[as.character(r[i])]])
          if(length(cr) == 0) return(0)
          backbone$cost[cr, "cost"] * adjust[i]
        }))
      })
      leaves <- leavecombs[which.max(cost),]
      rm(comp_rows, leave_rows, leavecombs, cost)
    }
    else if(length(leaves) != length(backbone$includedE))
      stop("length of leaves must be the number of components in the original graph")
    edges <- character()
    for(c in seq(length(leaves))){
      if(is.na(leaves[c]) & length(backbone$includedE[[c]]) > 1) edges <- c(edges, unlist(backbone$includedE[[c]][[1]]))
      else if(leaves[c] > 1) edges <- c(edges, unlist(backbone$includedE[[c]][1:leaves[c]]))
    }
    backbone$B <- subgraph.edges(backbone$pine, get.edge.ids(backbone$pine, edges))
  }
  if(assign | !is.null(backbone$branch)) backbone <- decompose(backbone, G=G, prob=prob, col=col)
  return(backbone)
}
