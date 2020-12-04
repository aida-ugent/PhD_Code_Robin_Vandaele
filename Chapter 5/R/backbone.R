#' @param X             data for which a backbone is to be computed, either in an 'igraph' format, or a 'matrix'
#'                      or 'data.frame' format, or a 'dist' object, in the latter cases a proximity graph is constructed
#' @param k             number of neighbors to construct a kNN graph, standard set to 10
#' @param eps           distance parameter for Vietoris-Rips graph, must be provided if X is not in 'igraph' format and type is 'Rips'
#' @param type          used if X is not in 'igraph' format, must be either 'NN' if kNN graph is to be constructed,
#'                      or 'Rips' if Vietoris-Rips graph is to be constructed, standard is 'NN'
#' @param f             optional function values to compute f-pine, boundary coefficients are computed and used if missing
#' @param times         optional number of times the pine must be pruned, only used if method is 'hard' or 'soft', standard set to 1
#'                      may also be a vector with size the number of components of the (constructed) graph
#' @param max_leaves    optional upper bound on the number of leaves to be included in the backbone
#'                      may also be a vector with size the number of components of the (constructed) graph
#'                      if method is 'search', a number of leaves of at most 'max_leaves' is estimated
#'                      standard set to NA (which is equivalent to +infinity)
#' @param leaves        optional number of leaves to be included in the backbone, 'max_leaves' and 'times' are ignored if provided
#'                      may also be a vector with size the number of components of the (constructed) graph
#'                      if not provided: the 'max_leaves' parameter is used in stead,
#'                      note that the specified number of leaves may not be attainable, especially when using the 'hard' pruning method
#' @param method        used method to perform the pruning, standard set to 'hard':
#'                         - 'hard': iteratively prune all leaves in the graph
#'                         - 'soft': iteratively prune the leaves representing the fewest nodes
#'                         - 'search': search for a set of k leaves resulting in the highest cost tree with such many leaves
#' @param stdize        whether to standardize the cost of each component according to its full cost, standard set to FALSE
#' @param assign        whether or not to produce a branch assignment, standard is FALSE
#' @param prob          (if assign) whether to assign nodes in G to paths according to a probability measure
#'                      if FALSE, then nodes are assigned to exactly one path, standard is TRUE
#' @param col           (if assign) whether to provide a coloring scheme for the paths and nodes, standard set to TRUE
#' @param bb_only       if TRUE then only the backbone subgraph is returned, standard is FALSE
#' @param constrain_mem whether to compute the entire distance matrix on G as an intermediate step to speed up the (currently inefficient) computation
#'                      for the boundary coefficients, the cost of more memory usage, standard is FALSE
#'
#' @return if bb_only is TRUE, a list storing the following items is returned:
#'                - (if 'method' is 'search' or either 'leaves' or 'max_leaves' is provided) B: the backbone in 'X'
#'                  encoded as an 'igraph' object
#'                - (if constructed) G: the constructed proximity graph from 'X'
#'                - (if computed) f: the boundary coefficients of the (constructed) graph
#'                - pine: an f-pine in the (constructed) graph
#'                - membership: the membership of the 'X' points/nodes according to the components of the (constructed) graph
#'                - (if method is 'search') cost: a X.frame corresponding to the obtained cost according to the number of leaves and component
#'                - (if method is 'search') full_cost: a vector of the full cost of each component
#'                - (if method is 'search') includedE: a list of lists (one for each connected component) containing the added edges of 'G' at each iteration
#'                - (if 'assign') branch: a vector clustering the edges of B into paths
#'                - (if 'assign') prob: a n x k matrix, where n = |'X'| and k is the number of resulting paths in B,
#'                  corresponding to a probabilistic assignment of the nodes/points of 'X' according to these paths (if 'prob')
#'                  or a one-hot-encoding (if !'prob')
#'                - (if 'col') palette: a vector of hex code colors, one for each path
#'                - (if 'col') col: a vector of |X| hex code colors, interpolating between the palette colors according to prob

backbone <- function(X, k=10, eps=NULL, type="NN", pairG=NULL, f=NULL, times=1, max_leaves=NA, leaves=NA, method="search",
                     stdize=FALSE, preprune=TRUE, vcost="betweenness", assign=FALSE, prob=TRUE, col=TRUE, bb_only=FALSE, constrain_mem=FALSE){
  Old <- Sys.time()
  store_G <- FALSE
  G <- NULL
  if(!method %in% c("hard", "soft", "search")) stop("method must be 'hard', 'soft', or 'search', standard is 'search'")
  if(!vcost %in% c("degree", "betweenness")) stop("vcost must be 'degree', or 'betweenness', standard is 'betweenness'")
  if(!is_igraph(X)){
    if(!(is.data.frame(X) | is.matrix(X) | class(X)=="dist")) stop("X must be either an igraph object, a data frame,
                                                                   a matrix, or distance object")
    if(!type %in% c("NN", "Rips")) stop("type must be 'NN' or 'Rips', standard is 'NN'")
    if(type == "Rips" & is.null(eps)) stop("type is 'Rips', please provide a value for 'eps'")
    store_G <- TRUE
    old <- Sys.time()
    if(type == "NN") G <- kNNgraph(X, k)
    else G <- VRgraph(X, eps)
    new <- Sys.time() - old
    print(paste("Proximity graph constructed in", round(new, 3), attr(new, "units")))
  }
  else if(is.null(V(X)$name)) V(X)$name <- as.character(seq(length(V(X))))
  store_f <- FALSE
  if(is.null(f)){
    store_f <- TRUE
    old <- Sys.time()
    f <- boundary_coefficient(if(is.null(G)) X else G, constrain_mem=constrain_mem)
    new <- Sys.time() - old
    print(paste("Boundary coefficients obtained in", round(new, 3), attr(new, "units")))
  }
  old <- Sys.time()
  pine <- mst(if(is.null(G)) X else G, weights = sum_weights(if(is.null(G)) X else G, f))
  new <- Sys.time() - old
  print(paste("f-pine obtained in", round(new, 3), attr(new, "units")))
  old <- Sys.time()
  if(any(is.na(leaves))){
    if(vcost == "degree") backbone <- prune(pine, leaves=max_leaves + 1, times=times, method=method, vcost=degree(pine), preprune=preprune)
    else backbone <- prune(pine, leaves=max_leaves + 1, times=times, method=method, vcost=betweenness(pine, weights=NA), preprune=preprune)
    backbone$pine <- pine
    rm(pine)
    if(method != "search"){
      new <- Sys.time() - old
      print(paste("Backbone obtained in", round(new, 3), attr(new, "units")))
    }
    else{
      backbone <- get_new_leaves(backbone, leaves=get_elbow(backbone$cost)[,1], G=if(is.null(G)) X else G, assign=FALSE, stdize=stdize)
      new <- Sys.time() - old
      print(paste("Leaves tuned and backbone obtained in", round(new, 3), attr(new, "units")))
    }
  }
  else{
    if(vcost == "degree") backbone <- prune(pine, leaves=leaves, method=method, vcost=degree(pine), preprune=preprune)
    else backbone <- prune(pine, leaves=leaves, method=method, vcost=betweenness(pine, weights=NA), preprune=preprune)
    backbone$pine <- pine
    if(length(leaves) == 1 & length(levels(factor(backbone$membership))) > 1)
      backbone <- get_new_leaves(backbone, leaves=leaves, G=if(is.null(G)) X else G, assign=FALSE, stdize=stdize)
    new <- Sys.time() - old
    print(paste("Backbone obtained in", round(new, 3), attr(new, "units")))
  }
  if(bb_only){
    if(method=="search") return(backbone$B)
    else return(backbone)
  }
  if(assign){
    old <- Sys.time()
    backbone <- decompose(backbone, G=if(is.null(G)) X else G, prob=prob, col=col)
    new <- Sys.time() - old
    print(paste("Branch assignment obtained in", round(new, 3), attr(new, "units")))
  }
  if(store_f) backbone$f <- f
  backbone$G <- G
  new <- Sys.time() - Old
  print("--------------------------------------------")
  print(paste("Backbone pipeline conducted in", round(new, 3), attr(new, "units")))
  return(backbone)
}
