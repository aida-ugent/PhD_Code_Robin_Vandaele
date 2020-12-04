#' Compute 0-dimensional persistence of a graph
#'
#' @param FltFun   optional filtration for which to compute persistence, if missing, 
#'                 the parameters for computing the filtration must be provided
#' @param S        optional graph represented by an igraph or cmplx object
#' @param f        optional vertex valued function defining the filtration
#' @param sublevel if FALSE and FltFun is missing, superlevel filtration is computed in stead, standard is TRUE
#'
#' @return a filtration

graphPersistence <- function(FltFun=NULL, S=NULL, f=NULL, sublevel=TRUE){
  if(is.null(FltFun) & (is.null(S) | is.null(f))) stop("Either filtration, or graph and vertex-valued function, must be provided")
  if(is.null(FltFun)) FltFun <- graphFiltration(S, f, sublevel)
  diag <- TDA::filtrationDiag(filtration=FltFun, maxdimension=0, library="Dionysus")$diagram
  return(diag)
}
