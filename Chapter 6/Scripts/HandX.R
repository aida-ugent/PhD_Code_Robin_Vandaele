# set working directory

setwd("/home/robin/Documents/LeaveEstimator")

# Load libraries

library("ggplot2") # plotting
library("ggpubr") # plotting in grid
library("igraph") # working with graphs

# Source functions

for(f in list.files("R")) source(file.path("R", f))

# Define topologies

Xshape <- list(list(start=c(0,0), end=c(1,1)),
               list(start=c(0,0), end=c(-1,1)),
               list(start=c(0,0), end=c(1,-1)),
               list(start=c(0,0), end=c(-1,-1)))
Hshape <- list(list(start=c(-1/10,0), end=c(1/10,0)),
               list(start=c(-1/10,0), end=c(-1,1)),
               list(start=c(-1/10,0), end=c(-1,-1)),
               list(start=c(1/10,0), end=c(1,1)),
               list(start=c(1/10,0), end=c(1,-1)))

shapes <- list("H"=Hshape, "X"=Xshape)

# Sample uniformly from topologies

set.seed(69)

npoints <- 300 # exact number of points may slightly differ due to rounding
cleanDataSets <- lapply(shapes, function(shape){
  cleanData <- data.frame(x=numeric(), y=numeric())
  pointsPerBranch <- sapply(shape, function(l) norm(l$start - l$end, type="2"))
  pointsPerBranch <- round(pointsPerBranch / sum(pointsPerBranch) * npoints)
  for(idx in 1:length(shape)){
    branch <- shape[[idx]]
    t <- runif(pointsPerBranch[[idx]])
    pointsOnBranchX <- t * branch$start[1] + (1 - t) * branch$end[1]
    pointsOnBranchY <- t * branch$start[2] + (1 - t) * branch$end[2]
    cleanData <- rbind(cleanData, data.frame(x=pointsOnBranchX, y=pointsOnBranchY))
  }
  return(cleanData)
})

# Plot clean data

cleanLim = c(min(unlist(cleanDataSets)), max(unlist(cleanDataSets)))
plotsOfCleanData <- lapply(cleanDataSets, function(cleanData){
  ggplot(cleanData, aes(x=x, y=y)) +
    geom_point() +
    xlim(cleanLim) +
    ylim(cleanLim) +
    coord_fixed()
})
ggpubr::ggarrange(plotlist=plotsOfCleanData)

# Add random noise and rotation, center to (0, 0)

nDistortions <- 1
sigma <- 0.005
noisyDataSets <- lapply(cleanDataSets, function(cleanData){
  lapply(seq_len(nDistortions), function(distortionIdx){
    noise <- MASS::mvrnorm(nrow(cleanData), mu=c(0,0), Sigma=diag(ncol(cleanData))*sigma)
    noisyData <- data.frame(cleanData + noise)
    colnames(noisyData) <- colnames(cleanData)
    return(noisyData)
  })
})
names(noisyDataSets) <- names(shapes)

# Plot noisy data

noisyLim <- c(min(unlist(noisyDataSets)), max(unlist(noisyDataSets)))
plotsOfNoisyData <- lapply(1:length(noisyDataSets), function(idx){
  noisyDatasForShape <- noisyDataSets[[idx]]
  groundTruth <- data.frame(do.call("rbind", lapply(shapes[[idx]], function(l) c(l$start, l$end))))
  colnames(groundTruth) <- c("x1", "y1", "x2", "y2")
  lapply(noisyDatasForShape, function(noisyData){
    ggplot(noisyData, aes(x=x, y=y)) +
      geom_point() +
      geom_segment(data=groundTruth, aes(x=x1, y=y1, xend=x2, yend=y2), col="red", size=2) +
      xlim(noisyLim) +
      ylim(noisyLim) +
      coord_fixed() +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
  })
})
ggpubr::ggarrange(plotlist=do.call("c", plotsOfNoisyData))

# Persistent homology of X topology using a Vietoris-Rips filtration

filtration <- TDA::ripsFiltration(noisyDataSets$X[[1]], maxdimension=1, maxscale=max(dist(noisyDataSets$X[[1]])))
diag <- TDA::filtrationDiag(filtration, maxdimension=1, library="Dionysus", printProgress=TRUE)
TDA::plot.diagram(diag[["diagram"]], diagLim=c(0,  2.5))
legend(2, 1, legend=c("H0", "H1"), col=c("black", "red"), pch=c(19, 2), pt.lwd=2, box.lty=0)

# Persistent homology of X topology through proximity graph

G <- kNNgraph(noisyDataSets$X[[1]])
EG <- get_edges2D(noisyDataSets$X[[1]], G)
eccG <- apply(distances(G), 1, max)
ggplot(noisyDataSets$X[[1]][names(V(G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="red", size=0.5) +
  geom_point(aes(col=eccG)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="Ecc") +
  coord_fixed()

diagG <- graphPersistence(S=G, f=eccG, sublevel=FALSE)
TDA::plot.diagram(diagG, diagLim=c(1.5,  3))
legend(2, 1, legend=c("H0"), col=c("black"), pch=c(19, 2), pt.lwd=2, box.lty=0)
