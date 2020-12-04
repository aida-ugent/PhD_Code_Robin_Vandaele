#' Check for any gaps in the obtained backbone using persistent homology
#'
#' @param backbone  a backbone object for which gaps need to be identified
#' @param G         the original graph on top of the backbone, standard set to NULL
#'
#' @return a backbone object storing information to investigate whether there are any gaps present

check_for_cycles <- function(backbone, G=NULL, maxscale=1){
  old <- Sys.time()
  persistence <- TDA::ripsDiag(X=distances(if(is.null(G)) backbone$G else G, v=V(backbone$B)$name, to=V(backbone$B)$name),
                               maxdimension=1, maxscale=Inf, dist="arbitrary", location=TRUE, library="Dionysus")
  new <- Sys.time() - old
  print(paste("Persistent homology computation performed in", round(new, 3), attr(new, "units")))
  backbone$diagram <- persistence$diagram
  backbone$repCycles <- lapply(persistence$cycleLocation[persistence$diagram[,"dimension"]==1], function(CL){
    return(cbind(V(backbone$B)$name[CL[,1]], V(backbone$B)$name[CL[,2]]))
  })
  return(backbone)
}
