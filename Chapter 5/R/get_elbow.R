#' Compute the elbow(s) for curve(s), the function is specifically designed to process the curves outputed from the search pruning
#'
#' @param C     data frame or matrix containing the x-coordinates (col1) and y-coordinates (col2) of the curve(s)
#'              optionally, a third column may be provided to group the coordinates according to different curves (col3)
#' @param alpha optional threshold between 0 and 1, denoting a sensitivity to adding new leaves
#'              if not provided, the elbow is computed by minimizing the absolute central differences
#'
#' @return integer vector containing the computed elbow index for each curve

get_elbow <- function(C, alpha=NA){
  elbows <- sapply(levels(C[,"component"]), function(c){
    I <- which(C[,"component"] == c)
    if(length(I) <= 3) return(c(NA, NA))
    m <- min(C[I, "leaves"])
    v <- rep(NA, max(C[I, "leaves"]) - m + 1)
    v[C[I, 1] - m + 1] <- C[I, "cost"]
    I <- which(is.na(v))
    while(length(I) > 0){
      i <- I[1] - 1
      j <- I[1] + 1
      while(is.na(v[j])) j <- j + 1
      invl <- 1 / (j - i)
      t <- seq(invl, 1 - invl, invl)
      v[(i + 1):(j - 1)] <- v[i] * (1 - t) + v[j] * t
      I <- which(is.na(v))
    }
    i <- which.min(diff(diff(v))[-1]) + 1
    return(c(i + m, v[i + 1]))
  })
  return(data.frame(leaves = as.integer(elbows[1,]), cost = elbows[2,]))
}
