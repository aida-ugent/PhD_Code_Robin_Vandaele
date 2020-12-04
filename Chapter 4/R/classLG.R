#' Classify points in a dataframe byusing their local and global topology
#' according to discrete graph metric defined by hop-count distances
#'
#' @param g an igraph object constructed from a point cloud data set 
#'          approaching a 1D stratified space
#' @param r a parameter denoting how far each point can see
#' 
#' @return a 2D dataframe with columns for degree & lower bound on cycles
 
classLG <- function(g, r){
  lgt <- data.frame(degree = integer(length(V(g))), cycles = integer(length(V(g))))
  to_class = as.integer(V(g))
  while(length(to_class) > 0){
    bfs <- bfs(g, root = to_class[1], dist = TRUE, unreachable = FALSE)
    deg <- components(induced.subgraph(g, which(bfs$dist == r)))$no
    cyc <- deg - components(induced.subgraph(g, which(bfs$dist >= r)))$no
    for(i in c(to_class[1], neighbors(g, to_class[1]))){
      lgt[i,] <- c(deg, cyc)
      to_class <- to_class[which(to_class != i)] 
    }
  }
  return(lgt)
}