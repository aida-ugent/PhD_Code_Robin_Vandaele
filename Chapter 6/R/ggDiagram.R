#' Construct a ggplot version of a TDA persistence diagram plot
#'
#' @param D          persistence diagram (TDA)
#' @param size       size of the points in the diagram
#' @param lim        optional numeric vector of length 2 denoting the (equal) limits of the x and y-axis
#' @param hideLegend whether to hide the legend denoting the corresponding dimensions of the topological features, standard is FALSE
#'
#' @output a ggplot version of the persistence diagram

ggDiagram <- function(diagram, size=4, lim=NULL, hideLegend=FALSE){
  if(is.null(lim))
    lim <- c(min(diagram[,c("Birth", "Death")]), max(diagram[abs(diagram[,"Birth"]) < Inf & abs(diagram[,"Death"]) < Inf, c("Birth", "Death")]))
  roundBreaks <- ceiling(abs(min(0, log(lim[2] - lim[1], base=10) - 1)))
  df.diagram <- data.frame("dimension"=integer(nrow(diagram)),
                           "Birth"=numeric(nrow(diagram)),
                           "Death"=numeric(nrow(diagram)))
  df.diagram[,"dimension"] <- factor(sapply(diagram[,"dimension"], function(dim) paste0("H", dim)))
  for(col in c("Birth", "Death")) df.diagram[,col] <- diagram[,col]
  if(hideLegend) legend.position="none" else legend.position="right"
  if(Inf %in% diagram){
    ggplot(df.diagram, aes(x=Birth, y=Death, col=dimension, shape=dimension)) +
      geom_point(size=size) +
      scale_x_continuous(limits=c(0.95 * lim[1], 1.05 * lim[2]), breaks=round(seq(lim[1], lim[2], length.out=4), roundBreaks)) +
      scale_y_continuous(limits=c(0.95 * lim[1], 1.05 * lim[2]), breaks=round(seq(lim[1], lim[2], length.out=4), roundBreaks)) +
      geom_abline(intercept=0, slope=1) +
      geom_hline(yintercept=Inf, linetype="dashed") +
      theme(legend.title=element_blank(), legend.position=legend.position, panel.background=element_rect(fill="white"), 
            axis.line=element_line(colour="black"), panel.grid.major=element_line(colour="grey"), panel.grid.minor=element_line(colour="lightgrey"))
  }
  else{
    ggplot(df.diagram, aes(x=Birth, y=Death, col=dimension, shape=dimension)) +
      geom_point(size=size) +
      scale_x_continuous(limits=c(0.95 * lim[1], 1.05 * lim[2]), breaks=round(seq(lim[1], lim[2], length.out=4), roundBreaks)) +
      scale_y_continuous(limits=c(0.95 * lim[1], 1.05 * lim[2]), breaks=round(seq(lim[1], lim[2], length.out=4), roundBreaks)) +
      geom_abline(intercept=0, slope=1) +
      theme(legend.position=legend.position, legend.title=element_blank(), panel.background=element_rect(fill="white"), 
            axis.line=element_line(colour="black"), panel.grid.major=element_line(colour="grey"), panel.grid.minor=element_line(colour="lightgrey"))
  }
}
