#' Extract a subtree of a connected pine by iteratively including the shortest path to the point providing the best cost increase
#'
#' @param G      connected pine to be pruned
#' @param leaves number of leaves to be kept, standard set to NA, in which case the search continues until either all nodes are included,
#'               or the addition of further nodes would be random
#' @param weight optional reweighting vector of the edges of G
#' @param vcost  optional vertex cost for which the sum is to be optimized (ignores weight if provided)
#'
#' @return a list storing the following items:
#'                - B: the pruned graph
#'                  encoded as an 'igraph' object
#'                - cost: a data.frame corresponding to the obtained cost according to the number of leaves
#'                - (if method is 'search') full_cost: a vector of the full cost of each component
#'                - (if vcost is not provided) includedE: a list containing the added edges in the tree at each iteration
#'                - (if vcost is provided) includedV: a list containing the added nodes in the tree at each iteration

search_prune <- function(G, leaves=NA, weight=NULL, vcost=NULL, preprune=TRUE){
  if(length(V(G)) == 0) stop("G must not be the empty graph")
  if(length(E(G)) != length(V(G)) - 1 | any(is.na(bfs(G, 1, unreachable=FALSE)$order))) stop("G is not a connected tree")
  if(is.na(leaves)) leaves <- Inf
  if(is.null(vcost)){
    if(is.null(weight)){
      if(is.null(E(G)$weight)) full_cost <- length(E(G))
      else full_cost <- sum(E(G)$weight)
    }
    else{
      if(!is.null(weight) & any(weight<0)) stop("weights must be nonnegative")
      full_cost <- sum(weight)
    }
  }
  else{
    if(any(vcost < 0)){
      print(vcost)
      stop("vcost must be nonnegative")
    }
    full_cost <- sum(vcost)
  }
  pruned_once <- FALSE
  ends <- V(G)$name[degree(G) == 1]
  to <- NULL
  if(preprune & length(V(G)) > 2 & ((is.null(weight) & is.null(E(G)$weight) & is.null(vcost)) | length(unique(vcost[ends])) == 1)){
    pruned_once <- TRUE
    value <- if(is.null(vcost)) 1 else vcost[ends[1]]
    to <- setdiff(V(G)$name, ends)
    H <- delete_vertices(G, ends)
    ends <- V(H)$name[degree(H) == 1]
    rm(H)
  }
  if((length(V(G)) == 2 & !pruned_once) | length(ends) == 0 | leaves < 2){ # special cases
    l <- 0
    includedE <- NULL
    includedV <- NULL
    if(is.null(vcost) | full_cost == 0) cost <- 0
    else{
      i <- which.max(vcost[if(is.null(to)) V(G) else V(G)[to]])
      includedV <- list((if(is.null(to)) V(G)$name else V(G)[to]$name)[i])
      cost <- (vcost[if(is.null(to)) V(G) else V(G)[to]])[i] / full_cost
    }
    if(leaves > 1 & length(ends) > 0){
      l <- c(l, 2)
      cost <- c(cost, 1)
      if(is.null(vcost)) includedE <- list(integer(0), t(ends(G, E(G))))
      else includedV[[2]] <- V(G)$name
    }
    return(list(B=if(is.null(vcost)) subgraph.edges(G, get.edge.ids(G, t(ends(G, includedE[[length(includedE)]])))) else
      induced.subgraph(G, includedV[[length(includedV)]]), cost=data.frame("leaves"=l, "cost"=cost), full_cost=full_cost,
                includedV=includedV, includedE=includedE))
  }
  if(is.null(vcost)){
    store_nodes <- FALSE
    vcost <- integer(length(V(G)))
    D <- distances(G, v=ends, to=if(is.null(to)) V(G) else to, weight=weight)
  } else{
    store_nodes <- TRUE
    D <- (distances(G, v=ends, to=if(is.null(to)) V(G) else to, weight=sum_weights(G, vcost)) +
            sapply(vcost[if(is.null(to)) V(G) else to], function(c) c + vcost[ends])) / 2
  }
  if(is.null(names(vcost))) names(vcost) <- names(V(G))
  rowMaxInd <- apply(D[ends, ends], 1, which.max)
  u <- which.max(D[ends, ends][cbind(seq(length(ends)), rowMaxInd)])
  v <- ends[rowMaxInd[u]]
  u <- ends[u]
  P <- shortest_paths(G, from=u, to=v, output="both")
  if(store_nodes){
    i <- which.max(vcost[if(is.null(to)) V(G) else V(G)[to]])
    cost <- c((vcost[if(is.null(to)) V(G) else V(G)[to]])[i] / full_cost, D[u, v] / full_cost)
    includedE <- NULL
    includedV <- list((if(is.null(to)) V(G)$name else V(G)[to]$name)[i], P$vpath[[1]]$name)
  } else{
    cost <- c(0, D[u, v] / full_cost)
    includedE <- list(integer(0), t(ends(G, P$epath[[1]])))
    includedV <- NULL
  }
  l <- 2
  add_bif <- FALSE
  leaves_to_add <- setdiff(rownames(D)[!is.na(D[,u])], c(u, v))
  if(length(leaves_to_add) > 0 & l < leaves){
    add_bif <- TRUE
    if(length(leaves_to_add) == 1) closest_point <- names(P$vpath[[1]][which.min(D[leaves_to_add, names(P$vpath[[1]])])])
    else closest_point <- names(P$vpath[[1]][apply(D[leaves_to_add, names(P$vpath[[1]])], 1, which.min)])
    names(closest_point) <- leaves_to_add
    while(l < leaves & length(leaves_to_add) > 0){
      x <- leaves_to_add[which.max(D[cbind(leaves_to_add, closest_point[leaves_to_add])] - vcost[closest_point[leaves_to_add]])]
      y <- closest_point[x]
      P <- shortest_paths(G, from=x, to=y,  output="both")
      if(l == 2) cost <- c(cost[1], max(D[u, y], D[v, y]) / full_cost, cost[2])
      l <- l + 1
      cost <- c(cost, cost[length(cost)] + (D[x, y] - vcost[y]) / full_cost)
      if(store_nodes) includedV[[l]] <- P$vpath[[1]]$name[-length(P$vpath[[1]])]
      else includedE[[l]] <- t(ends(G, P$epath[[1]]))
      leaves_to_add <- setdiff(leaves_to_add, x)
      for(x in leaves_to_add){
        y <- names(P$vpath[[1]][which.min(D[x, names(P$vpath[[1]])] - vcost[P$vpath[[1]]])])
        if(D[x, y] - vcost[y] < D[x, closest_point[x]] - vcost[closest_point[x]]) closest_point[x] <- y
      }
    }
  }
  if(pruned_once){
    cost <- c(cost, cost[length(cost)] + value / full_cost)
    l <- l + 1
  }
  if(store_nodes) G <- induced.subgraph(G, unlist(includedV)[-1])
  else G <- subgraph.edges(G, get.edge.ids(G, unlist(includedE)))
  return(list(B=G, cost=data.frame("leaves"=c(0, if(add_bif) 1 else NULL, 2:l), "cost"=cost),
              full_cost=full_cost, includedV=includedV, includedE=includedE))
}
