#' Construct Vietoris-Rips graph from a point cloud data set by connecting all
#' nodes with a certain given distance eps
#'
#' @param trajectory The trajectory
#' @param dimred     The dimensionality reduction matrix or function which will run the dimensionality reduction
#' @param theta,phi	 The angles defining the viewing direction. theta gives the azimuthal direction and phi the colatitude.
#'
#' @return the Vietoris-Rips graph defined by X and distance parameter eps

plotTrajectory3D <- function(trajectory, dimred, theta=0, phi=0){
  if(is.character(dimred)) dimred <- dyndimred::dimred(trajectory$expression(), method=dimred, ndim=3)
  grouping <- factor(trajectory$grouping)
  set.seed(42)
  plot3D::scatter3D(x=dimred[,1], y=dimred[,2], z=dimred[,3], theta=theta, phi=phi, ticktype="detailed", bty="b2", pch=19, 
                    col.var=as.integer(grouping), col=randomcoloR::distinctColorPalette(length(levels(grouping))),
                    colkey=list(at=seq(min(newDimred[,3]), max(newDimred[,3]), length.out=length(levels(grouping))), 
                                labels=paste0("G", 1:length(levels(grouping))), side=1, length=0.5))
  avgCellPerGroup <- t(sapply(levels(grouping), function(l) colMeans(dimred[trajectory$grouping==l,])))
  for(idx in 1:nrow(trajectory$milestone_network)){
    u <- avgCellPerGroup[trajectory$milestone_network[idx, "from"][[1]],]
    v <- avgCellPerGroup[trajectory$milestone_network[idx, "to"][[1]],]
    plot3D::scatter3D(x=c(u[1], v[1]), y=c(u[2], v[2]), z=c(u[3], v[3]), col="black", pch=19, cex=2, add=TRUE)
    plot3D::scatter3D(x=c(u[1], v[1]), y=c(u[2], v[2]), z=c(u[3], v[3]), col="black", type="l", lwd=5, add=TRUE)
  }
}