# set working directory

setwd("/home/robin/Documents/LeaveEstimator")

# Load libraries

library("ggplot2") # plotting
library("igraph") # working with graphs

# Source functions

for(f in list.files("R")) source(file.path("R", f))

# Define ground-truth topology

Ishape <- list(start=c(-1,0), end=c(1,0))

# Sample from topologies

set.seed(69)

npoints <- 500
t <- runif(npoints)
t[1:3] <- c(0, 1/2, 1) # make sure there is a point a each critical point
cleanData <- data.frame(
  x=t * Ishape$start[1] + (1 - t) * Ishape$end[1],
  y=t * Ishape$start[2] + (1 - t) * Ishape$end[2]
)

ggplot(cleanData, aes(x=x, y=y)) +
  geom_point()

# Add a small amount of random noise, construct a proximity graph and the eccentricity function
# Compare to ground-truth topology

groundTruthData <- cleanData + matrix(rep(c(0,1/4), npoints), ncol=2, byrow=TRUE)
groundTruthGraph <- emstreeR::ComputeMST(groundTruthData)[,c("from", "to", "distance")]
colnames(groundTruthGraph)[3] <- "weight"
groundTruthGraph <- graph_from_data_frame(groundTruthGraph, directed=FALSE)
EgroundTruthGraph <- get_edges2D(groundTruthData, groundTruthGraph)
eccGroundTruth <- apply(distances(groundTruthGraph), 1, max)

sigma <- 0.0001
noisyData <- cleanData + MASS::mvrnorm(nrow(cleanData), mu=c(0,0), Sigma=diag(ncol(cleanData))*sigma) + 
  matrix(rep(c(0,-1/4), npoints), ncol=2, byrow=TRUE)
proximityGraph <- VRgraph(noisyData, 2 * max(dist(noisyData)))
EproximityGraph <- get_edges2D(noisyData, proximityGraph)
eccProximityGraph <- apply(distances(proximityGraph), 1, max)

ggplot(data=NULL, aes(x=x, y=y)) +
  geom_segment(data=EgroundTruthGraph, aes(x=x1, y=y1, xend=x2, yend=y2)) +
  geom_segment(data=EproximityGraph, aes(x=x1, y=y1, xend=x2, yend=y2), size=0.125, alpha=0.5) +
  geom_point(data=groundTruthData[V(groundTruthGraph)$name,], aes(col=-eccGroundTruth), size=2) +
  geom_point(data=noisyData[V(proximityGraph)$name,], aes(col=-eccProximityGraph), size=0.125) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="-Ecc") +
  theme_bw() +
  theme(text=element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

diagEccGT <- graphPersistence(S=groundTruthGraph, f=-eccGroundTruth)
diagEccPG <- graphPersistence(S=proximityGraph, f=-eccProximityGraph)
TDA::plot.diagram(rbind(diagEccGT, c(1, diagEccPG[,c("Birth", "Death")])), diagLim=c(-2,  0))
legend(-1, -1, legend=c("H0, Ground Truth", "H0, Proximity Graph"), col=c("black", "red"), 
       pch=c(19, 2), pt.lwd=2, box.lty=0, title=expression(bold("-Ecc")))

# Comparison for arbitrary functions

f <- rep(-1, npoints)
maxf <- which.min(apply(groundTruthData, 1, function(r) norm(r - c(0, 1/4), type="2")))
f[maxf] <- 1

g <- rep(-1, npoints)
maxg <- which.min(apply(noisyData, 1, function(r) norm(r - c(0, -1/4), type="2")))
g[maxg] <- 1
proximityGraph2 <- VRgraph(noisyData, 0.1)
EproximityGraph2 <- get_edges2D(noisyData, proximityGraph2)

ggplot(data=NULL, aes(x=x, y=y)) +
  geom_segment(data=EgroundTruthGraph, aes(x=x1, y=y1, xend=x2, yend=y2)) +
  geom_segment(data=EproximityGraph2, aes(x=x1, y=y1, xend=x2, yend=y2), size=0.125) +
  geom_point(data=groundTruthData, aes(col=f), size=2) +
  geom_point(data=noisyData, aes(col=g), size=0.125) +
  geom_point(data=rbind(groundTruthData[maxf,], noisyData[maxg,]), size=3.5) +
  geom_point(data=rbind(groundTruthData[maxf,], noisyData[maxg,]), size=2.5, aes(col=1)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="f, g") +
  theme_bw() +
  theme(text=element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

diagFGT <- graphPersistence(S=groundTruthGraph, f=f[as.integer(V(groundTruthGraph)$name)])
diagGPG <- graphPersistence(S=proximityGraph2, f=g[as.integer(V(proximityGraph2)$name)])
TDA::plot.diagram(rbind(diagFGT, c(1, diagGPG[,c("Birth", "Death")])), diagLim=c(-1,  1.5))
legend(0.25, 0.25, legend=c("H0, Ground Truth", "H0, Proximity Graph"), col=c("black", "red"), 
       pch=c(19, 2), pt.lwd=2, box.lty=0, title=expression(bold("f, g")))
