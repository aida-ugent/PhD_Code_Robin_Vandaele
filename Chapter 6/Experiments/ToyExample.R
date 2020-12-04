# set working directory

setwd("/home/robin/Documents/LeaveEstimator")

# Load libraries

library("ggplot2") # plotting
library("ggpubr") # plotting in grid
library("igraph") # working with graphs

# Source functions

for(f in list.files("R")) source(file.path("R", f))

# Define topologies

Ishape <- list(list(start=c(-1,-1), end=c(0,0)),
               list(start=c(0,0), end=c(1, 1)))
Yshape <- list(list(start=c(0,0), end=c(1,1)),
               list(start=c(0,0), end=c(-1,1)),
               list(start=c(0,0), end=c(0,-sqrt(2))))
Xshape <- list(list(start=c(0,0), end=c(1,1)),
               list(start=c(0,0), end=c(-1,1)),
               list(start=c(0,0), end=c(1,-1)),
               list(start=c(0,0), end=c(-1,-1)))
Hshape <- list(list(start=c(-1/2,0), end=c(0,0)),
               list(start=c(1/2,0), end=c(0,0)),
               list(start=c(-1/2,0), end=c(-1,1)),
               list(start=c(-1/2,0), end=c(-1,-1)),
               list(start=c(1/2,0), end=c(1,1)),
               list(start=c(1/2,0), end=c(1,-1)))

shapes <- list("I"=Ishape, "Y"=Yshape, "X"=Xshape, "H"=Hshape)

# Compute the true diagrams

trueGraphs <- lapply(shapes, function(shape){
  G <- t(sapply(shape, function(shape) do.call("c", shape)))
  weight <- apply(G, 1, function(r) norm(r[1:2] - r[3:4], type="2"))
  G <- graph_from_data_frame(data.frame(from=(apply(matrix(G[,1:2], ncol=2), 1, function(r) paste(r[1], r[2]))),
                                        to=apply(matrix(G[,3:4], ncol=2), 1, function(r) paste(r[1], r[2]))),
                             directed=FALSE)
  E(G)$weight <- weight
  G$DC <- degreeOfCentrality(G)
  return(G)
})
names(trueDiagrams) <- names(shapes)

trueDiagrams <- lapply(trueGraphs, function(G) graphPersistence(S=G, f=G$DC))

# Sample uniformly from topologies

set.seed(42)

npoints <- 600 # exact number of points may slightly differ due to rounding
cleanDataSets <- lapply(names(shapes), function(name){
  shape <- shapes[[name]]
  cleanData <- data.frame(x=numeric(), y=numeric(), DC=numeric())
  pointsPerBranch <- sapply(shape, function(l) norm(l$start - l$end, type="2"))
  pointsPerBranch <- round(pointsPerBranch / sum(pointsPerBranch) * npoints)
  for(idx in 1:length(shape)){
    branch <- shape[[idx]]
    t <- runif(pointsPerBranch[[idx]])
    pointsOnBranchX <- t * branch$start[1] + (1 - t) * branch$end[1]
    pointsOnBranchY <- t * branch$start[2] + (1 - t) * branch$end[2]
    DC <- t * trueGraphs[[name]]$DC[paste(branch$start[1], branch$start[2])] +
      (1 - t) * trueGraphs[[name]]$DC[paste(branch$end[1], branch$end[2])] 
    cleanData <- rbind(cleanData, data.frame(x=pointsOnBranchX, y=pointsOnBranchY, DC=DC))
  }
  return(cleanData)
})
names(cleanDataSets) <- names(shapes)

# Plot clean data

plotsOfCleanData <- lapply(names(cleanDataSets), function(name){
  ggplot(cleanDataSets[[name]], aes(x=x, y=y, col=DC)) +
    geom_point(size=2.5, col="black") +
    geom_point(size=2) +
    scale_colour_gradientn(colours=topo.colors(7)) +
    ggtitle(paste(name, "ground truth")) +
    labs(col="Normalized Centrality") +
    theme_bw() +
    theme(plot.title=element_text(hjust = 0.5),
          text=element_text(size=15),
          legend.key.width=unit(3, "line"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    coord_fixed()
})

ggpubr::ggarrange(plotlist=plotsOfCleanData, nrow=2, ncol=2)

# Add random noise and rotation, center to (0, 0)

nDistortions <- 3
sigma <- 0.005
noisyDataSets <- lapply(cleanDataSets, function(cleanData){
  lapply(seq_len(nDistortions), function(distortionIdx){
    noise <- MASS::mvrnorm(nrow(cleanData), mu=c(0,0), Sigma=diag(2)*sigma)
    rotation <- mixAK::rRotationMatrix(1, 2)
    noisyData <- data.frame(scale(t(apply(cleanData[,1:2], 1, function(row) rotation %*% row)) + noise, scale=FALSE))
    colnames(noisyData) <- colnames(cleanData)[1:2]
    return(noisyData)
  })
})
names(noisyDataSets) <- names(shapes)

# Plot noisy data

noisyLim <- c(min(unlist(noisyDataSets)), max(unlist(noisyDataSets)))
plotsOfNoisyData <- lapply(1:length(noisyDataSets), function(idx1){
  lapply(1:length(noisyDataSets[[idx1]]), function(idx2){
    ggplot(noisyDataSets[[idx1]][[idx2]], aes(x=x, y=y)) +
      geom_point(size=0.25) +
      xlim(noisyLim) +
      ylim(noisyLim) +
      coord_fixed()
  })
})

ggpubr::ggarrange(plotlist=do.call("c", plotsOfNoisyData), common.legend=TRUE, nrow=length(shapes), ncol=nDistortions)

# Construct MST and degree of centrality for each noisy data set

mstOfDataSets <- lapply(noisyDataSets, function(noisyDatasForShape){
  lapply(noisyDatasForShape, function(noisyData){
    G <- emstreeR::ComputeMST(noisyData[,1:2], verbose=FALSE)[,c("from", "to", "distance")]
    G <- data.frame(G)
    colnames(G)[3] <- "weight"
    G <- graph_from_data_frame(G, directed=FALSE)
    G$degreeOfCentrality <- degreeOfCentrality(G)
    return(G)
  })
})

# Plot clean data and MST with degrees of centrality for each noisy data set

plotsOfMST <- lapply(1:length(noisyDataSets), function(idx1){
  lapply(1:length(noisyDataSets[[idx1]]), function(idx2){
    euclideanEdges <- get_edges2D(noisyDataSets[[idx1]][[idx2]], mstOfDataSets[[idx1]][[idx2]])
    ggplot(noisyDataSets[[idx1]][[idx2]][names(V(mstOfDataSets[[idx1]][[idx2]])),], aes(x=x, y=y)) +
      geom_segment(data=euclideanEdges, aes(x=x1, y=y1, xend=x2, yend=y2), color="red", size=0.5) +
      geom_point(size=1.5, col="black") +
      geom_point(size=1, aes(col=mstOfDataSets[[idx1]][[idx2]]$degreeOfCentrality)) +
      ggtitle(paste(names(noisyDataSets)[[idx1]], "sample", idx2)) +
      scale_colour_gradientn(colours=topo.colors(7)) +
      labs(col="Normalized Centrality") +
      xlim(noisyLim) +
      ylim(noisyLim) +
      theme_bw() +
      theme(plot.title=element_text(hjust = 0.5),
            text=element_text(size=15),
            legend.key.width=unit(3, "line"),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      coord_fixed()
  })
})

ggpubr::ggarrange(plotlist=do.call("c", lapply(1:length(shapes), function(idx) c(list(plotsOfCleanData[[idx]]), plotsOfMST[[idx]]))), 
                  common.legend=TRUE, nrow=length(shapes), ncol=nDistortions + 1)

# Compute persistence diagrams of the sublevel filtrations induced by the degrees of centrality on the MSTs

diagramsOfMST <- list()
progress <- 0
print("Computing persistence diagrams...")
for(idx1 in 1:length(mstOfDataSets)){
  diagramsOfMST[[idx1]] <- list()
  for(idx2 in 1:length(mstOfDataSets[[idx1]])){
    print(paste0("progress: ", 100 * progress / (length(mstOfDataSets) * length(mstOfDataSets[[1]])), "%"))
    G <- mstOfDataSets[[idx1]][[idx2]]
    diagramsOfMST[[idx1]][[idx2]] <- graphPersistence(S=G, f=G$degreeOfCentrality)
    progress <- progress + 1
  }
  if(idx1==length(mstOfDataSets)) print(paste0("progress: 100%"))
}

# Plot persistence diagrams of the sublevel filtration induced by the degrees of centrality on the MSTs

colors <- randomcoloR::distinctColorPalette(length(noisyDataSets))
par(mfrow=c(length(noisyDataSets), length(noisyDataSets[[1]]) + 1), mar=c(3,4,2,2))
for(idx1 in 1:length(noisyDataSets)){
  TDA::plot.diagram(trueDiagrams[[idx1]], main=paste(names(shapes)[idx1], "ground truth"), diagLim=c(0,1.1), col=colors[idx1])
  for(idx2 in 1:length(noisyDataSets[[idx1]])){
    TDA::plot.diagram(diagramsOfMST[[idx1]][[idx2]], main=paste(names(shapes)[idx1], "sample", idx2), diagLim=c(0,1.1), col=colors[idx1])
  }
}
dev.off()

# Compute and view pairwise bottleneck distances for all persistence diagrams (MST)

pairwiseBottleNecksMST <- matrix(numeric((length(noisyDataSets)*(length(noisyDataSets[[1]]) + 1))^2),
                                 nrow=length(noisyDataSets)*(length(noisyDataSets[[1]]) + 1))
for(idx1 in 1:(nrow(pairwiseBottleNecksMST) -  1)){
  for(idx2 in idx1:(nrow(pairwiseBottleNecksMST))){
    shapeIndex1 <- ceiling(idx1 / (nDistortions + 1))
    dataIndex1 <- idx1 - ((shapeIndex1 - 1) * (nDistortions + 1))
    if(dataIndex1 == 1) dgm1 <- trueDiagrams[[shapeIndex1]]
    else dgm1 <- diagramsOfMST[[shapeIndex1]][[dataIndex1 - 1]]
    shapeIndex2 <- ceiling(idx2 / (nDistortions + 1))
    dataIndex2 <- idx2 - ((shapeIndex2 - 1) * (nDistortions + 1))
    if(dataIndex2 == 1) dgm2 <- trueDiagrams[[shapeIndex2]]
    else dgm2 <- diagramsOfMST[[shapeIndex2]][[dataIndex2 - 1]]
    pairwiseBottleNecksMST[idx1, idx2] <- TDA::bottleneck(dgm1, dgm2, dimension=0)
    pairwiseBottleNecksMST[idx2, idx1] <- pairwiseBottleNecksMST[idx1, idx2]
  }
}

ggplot(reshape2::melt(pairwiseBottleNecksMST), aes(Var1, Var2, fill=value)) +
  geom_raster() +
  scale_x_continuous(breaks=1 + seq(0, nrow(pairwiseBottleNecksMST) - length(shapes), 
                                    length.out=length(shapes)), labels=names(shapes)) +
  scale_y_continuous(breaks=1 + seq(0, nrow(pairwiseBottleNecksMST) - length(shapes), 
                                   length.out=length(shapes)), 
                     trans="reverse", labels=names(shapes)) +
  xlab("") +
  ylab("") +
  labs(fill="BN Dist.") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5),
        text=element_text(size=15)) +
  coord_fixed()

# View MDS embedding of pairwise bottleneck distances (MST)

fitBottleNecksMST <- data.frame(cmdscale(pairwiseBottleNecksMST, k=2))
fitBottleNecksMST[,3] <- rep(names(noisyDataSets), each=length(noisyDataSets[[1]]) + 1)
colnames(fitBottleNecksMST) <- c("x", "y", "shape")
ggplot(fitBottleNecksMST, aes(x=x, y=y, col=shape)) +
  geom_point(size=4.5, col="white") +
  geom_point(size=4) +
  geom_point(data=fitBottleNecksMST[1 + seq(0, nrow(fitBottleNecksMST) - length(shapes), 
                                        length.out=length(shapes)),], size=5, col="black") +
  geom_point(data=fitBottleNecksMST[1 + seq(0, nrow(fitBottleNecksMST) - length(shapes), 
                                            length.out=length(shapes)),], size=4) +
  theme_bw() +
  theme(plot.title=element_text(hjust = 0.5),
        text=element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_color_manual(values=colors[order(names(noisyDataSets))])

# Construct Rips graph and degree of centrality for each noisy data set

ripsGraphOfDataSets <- lapply(noisyDataSets, function(noisyDatasForShape){
  lapply(noisyDatasForShape, function(noisyData){
    G <- VRgraph(noisyData, eps=0.25)
    G$degreeOfCentrality <- degreeOfCentrality(G)
    return(G)
  })
})

# Plot Rips graphs with degrees of centrality for each noisy data set

plotsOfRipsGraphs <- lapply(1:length(noisyDataSets), function(idx1){
  lapply(1:length(noisyDataSets[[idx1]]), function(idx2){
    euclideanEdges <- get_edges2D(noisyDataSets[[idx1]][[idx2]], ripsGraphOfDataSets[[idx1]][[idx2]])
    ggplot(noisyDataSets[[idx1]][[idx2]][names(V(ripsGraphOfDataSets[[idx1]][[idx2]])),], aes(x=x, y=y)) +
      geom_segment(data=euclideanEdges, aes(x=x1, y=y1, xend=x2, yend=y2), color="red", size=0.025) +
      geom_point(size=0.25, aes(col=ripsGraphOfDataSets[[idx1]][[idx2]]$degreeOfCentrality)) +
      scale_colour_gradientn(colours=topo.colors(7)) +
      labs(col="Degree of Centrality") +
      xlim(noisyLim) +
      ylim(noisyLim) +
      coord_fixed()
  })
})
ggpubr::ggarrange(plotlist=do.call("c", plotsOfRipsGraphs), common.legend=TRUE)

# Compute persistence diagrams of the sublevel filtrations induced by the degrees of centrality on the Rips Graphs

diagramsOfRipsGraphs <- list()
progress <- 0
print("Computing persistence diagrams...")
for(idx1 in 1:length(ripsGraphOfDataSets)){
  diagramsOfRipsGraphs[[idx1]] <- list()
  for(idx2 in 1:length(ripsGraphOfDataSets[[idx1]])){
    print(paste0("progress: ", 100 * progress / (length(ripsGraphOfDataSets) * length(ripsGraphOfDataSets[[1]])), "%"))
    G <- ripsGraphOfDataSets[[idx1]][[idx2]]
    diagramsOfRipsGraphs[[idx1]][[idx2]] <- graphPersistence(S=G, f=G$degreeOfCentrality)
    progress <- progress + 1
  }
  if(idx1==length(mstOfDataSets)) print(paste0("progress: 100%"))
}

# Plot persistence diagrams of the sublevel filtration induced by the degrees of centrality on the Rips Graphs

par(mfrow=c(length(noisyDataSets), length(noisyDataSets[[1]])), mar=c(3,4,2,2))
for(idx1 in 1:length(noisyDataSets)){
  for(idx2 in 1:length(noisyDataSets[[idx1]])){
    TDA::plot.diagram(diagramsOfRipsGraphs[[idx1]][[idx2]])
  }
}
dev.off()

# Compute and view pairwise bottleneck distances for all persistence diagrams (Rips Graph)

pairwiseBottleNecksRipsGraphs <- matrix(numeric((length(noisyDataSets)*length(noisyDataSets[[1]]))^2),
                                        nrow=length(noisyDataSets)*length(noisyDataSets[[1]]))
for(idx1 in 1:(nrow(pairwiseBottleNecksRipsGraphs) -  1)){
  for(idx2 in idx1:(nrow(pairwiseBottleNecksRipsGraphs))){
    shapeIndex1 <- floor((idx1 - 1) / length(noisyDataSets[[1]])) + 1 
    dataIndex1 <- idx1 - (shapeIndex1 - 1) * length(noisyDataSets[[1]])
    shapeIndex2 <- floor((idx2 - 1) / length(noisyDataSets[[1]])) + 1
    dataIndex2 <- idx2 - (shapeIndex2 - 1) * length(noisyDataSets[[1]])
    pairwiseBottleNecksRipsGraphs[idx1, idx2] <- TDA::bottleneck(diagramsOfRipsGraphs[[shapeIndex1]][[dataIndex1]], 
                                                                 diagramsOfRipsGraphs[[shapeIndex2]][[dataIndex2]], 
                                                                 dimension=0)
    pairwiseBottleNecksRipsGraphs[idx2, idx1] <- pairwiseBottleNecksRipsGraphs[idx1, idx2]
  }
}

ggplot(reshape2::melt(pairwiseBottleNecksRipsGraphs), aes(Var1, Var2, fill=value)) +
  geom_raster() +
  scale_y_continuous(trans="reverse") +
  xlab("") +
  ylab("") +
  labs(fill="BN Dist.") +
  coord_fixed()

# View MDS embedding of pairwise bottleneck distances (Rips Graph)

fitBottleNecksRipsGraphs <- data.frame(cmdscale(pairwiseBottleNecksRipsGraphs, k=2))
fitBottleNecksRipsGraphs[,3] <- rep(names(noisyDataSets), each=length(noisyDataSets[[1]]))
colnames(fitBottleNecksRipsGraphs) <- c("x", "y", "shape")
ggplot(fitBottleNecksRipsGraphs, aes(x=x, y=y, col=shape)) +
  geom_point(size=3) +
  coord_fixed()

# Construct kNN graph and degree of centrality for each noisy data set

kNNGraphOfDataSets <- lapply(noisyDataSets, function(noisyDatasForShape){
  lapply(noisyDatasForShape, function(noisyData){
    G <- kNNgraph(noisyData, k=10)
    G$degreeOfCentrality <- degreeOfCentrality(G)
    return(G)
  })
})

# Plot kNN graphs with degrees of centrality for each noisy data set

plotsOfKNNGraphs <- lapply(1:length(noisyDataSets), function(idx1){
  lapply(1:length(noisyDataSets[[idx1]]), function(idx2){
    euclideanEdges <- get_edges2D(noisyDataSets[[idx1]][[idx2]], kNNGraphOfDataSets[[idx1]][[idx2]])
    ggplot(noisyDataSets[[idx1]][[idx2]][names(V(kNNGraphOfDataSets[[idx1]][[idx2]])),], aes(x=x, y=y)) +
      geom_segment(data=euclideanEdges, aes(x=x1, y=y1, xend=x2, yend=y2), color="red", size=0.025) +
      geom_point(size=0.25, aes(col=kNNGraphOfDataSets[[idx1]][[idx2]]$degreeOfCentrality)) +
      scale_colour_gradientn(colours=topo.colors(7)) +
      labs(col="Degree of Centrality") +
      xlim(noisyLim) +
      ylim(noisyLim) +
      coord_fixed()
  })
})
ggpubr::ggarrange(plotlist=do.call("c", plotsOfKNNGraphs), common.legend=TRUE)

# Compute persistence diagrams of the sublevel filtrations induced by the degrees of centrality on the Rips Graphs

diagramsOfKNNGraphs <- list()
progress <- 0
print("Computing persistence diagrams...")
for(idx1 in 1:length(kNNGraphOfDataSets)){
  diagramsOfKNNGraphs[[idx1]] <- list()
  for(idx2 in 1:length(kNNGraphOfDataSets[[idx1]])){
    print(paste0("progress: ", 100 * progress / (length(kNNGraphOfDataSets) * length(kNNGraphOfDataSets[[1]])), "%"))
    G <- kNNGraphOfDataSets[[idx1]][[idx2]]
    diagramsOfKNNGraphs[[idx1]][[idx2]] <- graphPersistence(S=G, f=G$degreeOfCentrality)
    progress <- progress + 1
  }
  if(idx1==length(kNNGraphOfDataSets)) print(paste0("progress: 100%"))
}

# Plot persistence diagrams of the sublevel filtration induced by the degrees of centrality on the KNN Graphs

par(mfrow=c(length(noisyDataSets), length(noisyDataSets[[1]])), mar=c(3,4,2,2))
for(idx1 in 1:length(noisyDataSets)){
  for(idx2 in 1:length(noisyDataSets[[idx1]])){
    TDA::plot.diagram(diagramsOfKNNGraphs[[idx1]][[idx2]])
  }
}
dev.off()

# Compute and view pairwise bottleneck distances for all persistence diagrams (KNN Graph)

pairwiseBottleNecksKNNGraphs <- matrix(numeric((length(noisyDataSets)*length(noisyDataSets[[1]]))^2),
                                        nrow=length(noisyDataSets)*length(noisyDataSets[[1]]))
for(idx1 in 1:(nrow(pairwiseBottleNecksKNNGraphs) -  1)){
  for(idx2 in idx1:(nrow(pairwiseBottleNecksKNNGraphs))){
    shapeIndex1 <- floor((idx1 - 1) / length(noisyDataSets[[1]])) + 1 
    dataIndex1 <- idx1 - (shapeIndex1 - 1) * length(noisyDataSets[[1]])
    shapeIndex2 <- floor((idx2 - 1) / length(noisyDataSets[[1]])) + 1
    dataIndex2 <- idx2 - (shapeIndex2 - 1) * length(noisyDataSets[[1]])
    pairwiseBottleNecksKNNGraphs[idx1, idx2] <- TDA::bottleneck(diagramsOfKNNGraphs[[shapeIndex1]][[dataIndex1]], 
                                                                 diagramsOfKNNGraphs[[shapeIndex2]][[dataIndex2]], 
                                                                 dimension=0)
    pairwiseBottleNecksKNNGraphs[idx2, idx1] <- pairwiseBottleNecksKNNGraphs[idx1, idx2]
  }
}

ggplot(reshape2::melt(pairwiseBottleNecksKNNGraphs), aes(Var1, Var2, fill=value)) +
  geom_raster() +
  scale_y_continuous(trans="reverse") +
  xlab("") +
  ylab("") +
  labs(fill="BN Dist.") +
  coord_fixed()

# View MDS embedding of pairwise bottleneck distances (KNN Graph)

fitBottleNecksKNNGraphs <- data.frame(cmdscale(pairwiseBottleNecksKNNGraphs, k=2))
fitBottleNecksKNNGraphs[,3] <- rep(names(noisyDataSets), each=length(noisyDataSets[[1]]))
colnames(fitBottleNecksKNNGraphs) <- c("x", "y", "shape")
ggplot(fitBottleNecksKNNGraphs, aes(x=x, y=y, col=shape)) +
  geom_point(size=3) +
  coord_fixed()
