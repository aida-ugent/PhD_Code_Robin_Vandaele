#' Compute possible leave configurations across multiple components for a given number of leaves
#'
#' @param c number of components
#' @param l number of leaves
#'
#' @return a matrix off c columns where each row is a possible combination of leaves across the c components

get_leave_combs <- function(c, l){
  if(l < 1) return(rep(0, c))
  else if(l < 2) return(diag(c))
  if(c==1) return(l)
  leavecombs <- t(sapply(c(0:l), function(k){
    if(k < l) return(cbind(k, get_leave_combs(c - 1, l - k)))
    return(c(k, rep(0, c - 1)))
  }))
  if(c >  2) leavecombs <- do.call("rbind", leavecombs)
  colnames(leavecombs) <- NULL
  return(leavecombs)
}
