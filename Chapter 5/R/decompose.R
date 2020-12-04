#' Decompose a subforest Fst in G in longest paths not passing through multifurcation points (branches) and assign all nodes to branches
#' The decomposition starts from an arbitrary leave within each tree of the forest Fst
#' If an edge in Fst connects two multifurcation points, its assignment will be randomly one of its neighboring edges in Fst
#'
#' @param backbone a backbone object, as returned from the 'backbone_pipe' function
#' @param G        the original graph, only needed if not stored in 'backbone'
#' @param prob     whether to assign nodes in G to paths according to a probability measure
#'                 if FALSE, then nodes are assigned to exactly one path (hard assignment), standard is TRUE
#' @param col      whether to provide a coloring scheme for the paths and nodes, standard set to TRUE
#'
#' @return the backbone object storing the following items:
#'                - branch: a vector clustering the edges of backbone$B into paths
#'                - prob: a (sparse) n x k matrix, where n = |V(G)| and k is the number of resulting paths in backbone$B,
#'                  corresponding to a probabilistic assignment of the nodes of G according to these paths (if prob)
#'                  or a one-hot-encoding (if !prob)
#'                - (if col) palette: a vector of hex code colors, one for each path
#'                - (if col) col: a vector of |V(G)| hex code colors, interpolating between the palette colors according to prob

decompose <- function(backbone, G=NULL, prob=TRUE, col=TRUE){
  if(is.null(backbone$G) & is.null(G)) stop("backbone does not contain the original graph, please provide G")
  distToB <- distances(if(is.null(G)) backbone$G else G, v=V(backbone$B)$name)
  decompV <- integer(length(V(backbone$B)))
  degrees <- degree(backbone$B)
  arm_nodes <- which(degrees < 3)
  decompV[arm_nodes] <- as.integer(components(induced_subgraph(backbone$B, arm_nodes))$membership)
  if(length(arm_nodes) != length(V(backbone$B))){
    decompV[-arm_nodes] <- sapply(which(decompV == 0), function(v)
      decompV[arm_nodes[which.min(distToB[v, V(backbone$B)$name[arm_nodes]])]])
  }
  decompE <- apply(get.edges(backbone$B, E(backbone$B)), 1, function(e){
    c <- decompV[e[degrees[e] %in% c(1, 2)]][1]
    if(is.na(c)) c <- decompV[e[1]]
    return(c)
  })
  BA <- sapply(V(backbone$pine), function(v) decompV[which.min(distToB[,v])])
  if(prob){
    BA_prob <- matrix(numeric(length(unique(BA)) * ncol(distToB)), nrow = ncol(distToB))
    for(i in seq(ncol(distToB))){
      if(!is.null(backbone$G)) f <- prop.table(table(BA[neighbors(backbone$G, V(backbone$G)[i])]))
      else f <- prop.table(table(BA[neighbors(G, V(G)[i])]))
      BA_prob[i, as.integer(names(f))] <- f
    }
  }
  else BA_prob <- as.matrix(diag(max(BA))[BA,])
  backbone$branch <- decompE
  backbone$prob <- as.spam(BA_prob)
  if(col){
    set.seed(42)
    backbone$palette <- randomcoloR::distinctColorPalette(ncol(BA_prob))
    backbone$col <- rgb(t(sapply(seq(nrow(BA_prob)), function(i)
      apply((BA_prob[i,] * t(sapply(backbone$palette, function(c) col2rgb(c))) / 255), 2, sum))))
    names(backbone$col) <- if(is.null(G)) names(V(BCB$G)) else names(V(G))
  }
  return(backbone)
}
