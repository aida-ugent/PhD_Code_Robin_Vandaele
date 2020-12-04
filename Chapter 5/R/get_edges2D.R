#' Given a two-dimensional dataframe df and a graph G = (V, E), where E is a list of edges
#' between nodes of df marked by their indexes, return a four-dimensional dataframe containing
#' the coordinates of the endpoints of each segment defined by df and E, suitable to plot with
#' geom_segment (ggplot2)
#'
#' @param df two-dimensional dataframe containing the points that must be connected through
#'              segments
#' @param G  Either a graph containing the edges (u, v) marking a segments between the u-th and
#'           v-th point of df, or a list of edges i
#'
#' @return   a four-dimensional dataframe containing the coordinates of the endpoints of the
#'           edges to plot

get_edges2D <- function(df, G){
  if(length(E(G)) == 0) return(NULL)
  if(ncol(df) != 2) stop('df must be two-dimensional')
  df.E <- ends(G, E(G))
  df.E <- data.frame(df[df.E[,1],], df[df.E[,2],])
  colnames(df.E) <- c("x1", "y1", "x2", "y2")
  rownames(df.E) <- 1:nrow(df.E)
  return(df.E)
}
