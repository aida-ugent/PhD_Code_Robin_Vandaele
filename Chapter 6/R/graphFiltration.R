#' Compute a filtration from a graph G
#'
#' @param S        a graph represented by an igraph or cmplx object
#' @param f        vertex valued function defining the filtration
#' @param sublevel if FALSE, superlevel filtration is computed in stead, standard is TRUE
#'
#' @return a filtration

graphFiltration <- function(S, f, sublevel=TRUE){
  if(class(S)=="igraph") S <- graphComplex(S)
  FltFun <- TDA::funFiltration(FUNvalues=f, cmplx=S, sublevel=sublevel)
  return(FltFun)
}
