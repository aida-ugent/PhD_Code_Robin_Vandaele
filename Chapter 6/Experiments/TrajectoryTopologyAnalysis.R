# set working directory

setwd("/home/robin/Documents/LeaveEstimator")

# Load libraries

library("ggplot2") # plotting
library("igraph") # working with graphs

# Source functions

for(f in list.files("R")) source(file.path("R", f))

# Load data and labels

dimredType <- "dm_diffusionmap" # type of dimensionality reduction
dimredDim <- 2 # number of dimensions of the embedding
graphType <- "15NN" # type of proximity graph, kNN or MST

dataFolder <- file.path("Data", dimredType)
signaturesFolder <- file.path("Signatures", paste0(dimredType, dimredDim), graphType)
graphsFolder <- file.path("Graphs", paste0(dimredType, dimredDim), graphType)
landmarksFolder <- file.path("Landmarks", paste0(dimredType, dimredDim)) 
if(!file.exists(landmarksFolder)) dir.create(landmarksFolder)
landmarksFolder <- file.path(landmarksFolder, graphType) 
if(!file.exists(landmarksFolder)) dir.create(landmarksFolder)

dynverseIDs <- sapply(list.files(signaturesFolder), function(f){
  prefix <- paste(strsplit(f, split="-")[[1]][c(1, 2)], collapse="/")
  suffix <- stringr::str_remove(f, stringr::str_replace(paste0(prefix, "-"), "/", "-"))
  return(paste(prefix, suffix, sep="/"))
})

leaves <- integer()
types <- character()
signatures <- list()
for(dataID in list.files(signaturesFolder)){
  leaves <- c(leaves, as.integer(read.table(file.path(dataFolder, dataID, "leaves.csv"))))
  types <- c(types, strsplit(dataID, split="-")[[1]][1])
  signatures[[length(signatures) + 1]] <- as.matrix(read.table(file.path(signaturesFolder, dataID, "persistence.csv")))
}

# Compute and view pairwise bottleneck distances

pwBN <- matrix(numeric(length(signatures)^2), nrow=length(signatures))
for(idx1 in 1:(nrow(pwBN) -  1)){
  for(idx2 in idx1:(nrow(pwBN))){
    pwBN[idx1, idx2] <- TDA::bottleneck(signatures[[idx1]], signatures[[idx2]], dimension=0)
    pwBN[idx2, idx1] <- pwBN[idx1, idx2]
  }
}

# View MDS embedding of pairwise bottleneck distances

set.seed(42)
fitBN <- data.frame(cmdscale(pwBN, k=2))
fitBN[,1] <- -fitBN[,1] # optional rescaling for consistent plots 
fitBN[,3] <- as.factor(sapply(leaves, function(l) if(l==2) "linear" else "nonlinear"))
fitBN[,4] <- as.factor(types)
colnames(fitBN) <- c("MDS1", "MDS2", "trajectory", "type")
ggplot(fitBN, aes(x=MDS1, y=MDS2, col=trajectory, shape=type)) +
  geom_point(size=3) +
  ggtitle(paste0(dimredType, dimredDim)) +
  theme(plot.title=element_text(hjust=0.5)) +
  coord_fixed()

# Optional: define landmarks to investigate the nearby topologies

landmarks <- rbind(data.frame(x=-0.4, y=0, linear=FALSE),
                   data.frame(x=-0.4, y=0, linear=TRUE),
                   data.frame(x=-0.2, y=0, linear=TRUE),
                   data.frame(x=-0.1, y=0.025, linear=TRUE),
                   data.frame(x=-0.05, y=0.05, linear=TRUE),
                   data.frame(x=0, y=0.1, linear=FALSE),
                   data.frame(x=0.125, y=0, linear=TRUE),
                   data.frame(x=0.1, y=0.05, linear=FALSE),
                   data.frame(x=0, y=-0.125, linear=FALSE),
                   data.frame(x=0.2, y=-0.15, linear=FALSE),
                   data.frame(x=-0.1, y=-0.15, linear=FALSE),
                   data.frame(x=-0.1, y=-0.25, linear=FALSE))
#write.table(landmarks, file.path(landmarksFolder, "landmarks.csv"))

# Investigate landmark topologies

isLinear <- which(leaves == 2)
isNotLinear <- which(leaves > 2)
landmarks <- read.table(file.path(landmarksFolder, "landmarks.csv"))

IDsOfInterest <- apply(landmarks, 1, function(row){
  if(row[3]) return(isLinear[which.min(apply(fitBN[isLinear, 1:2], 1, function(r) norm(r - row[1:2], type="2")))])
  else return(isNotLinear[which.min(apply(fitBN[isNotLinear, 1:2], 1, function(r) norm(r - row[1:2], type="2")))])
})

ggplot(fitBN, aes(x=MDS1, y=MDS2, shape=trajectory, col=type)) +
  geom_point(size=3) +
  geom_point(data=fitBN[IDsOfInterest,], col="black", size=5) +
  geom_point(data=fitBN[IDsOfInterest,], size=3) + 
  theme(plot.title=element_text(hjust=0.5)) +
  coord_fixed()

listOfPlots <- list()
for(idx in IDsOfInterest){
  TrajectoryID <- list.files(signaturesFolder)[idx]
  df.TrajectoryOfInterest <- read.table(file.path(dataFolder, TrajectoryID, "dimred.csv"))
  G.TrajectoryOfInterest <- graph_from_data_frame(read.table(file.path(graphsFolder, TrajectoryID, "graph.csv")), directed=FALSE)
  EG <- get_edges2D(df.TrajectoryOfInterest[,1:2], G.TrajectoryOfInterest)
  listOfPlots[[length(listOfPlots) + 1]] <- dynplot::plot_dimred(dynbenchmark::load_dataset(dynverseIDs[idx]), dimred=as.matrix(df.TrajectoryOfInterest)) +
    geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), alpha=0.05) +
    ggtitle(paste0("MDS (", round(fitBN[idx, 1], 2), ", ", round(fitBN[idx, 2], 2), ")"))
  listOfPlots[[length(listOfPlots) + 1]] <- ggDiagram(signatures[[idx]], size=2, lim=c(0, 1), hideLegend=TRUE)
}

ggpubr::annotate_figure(ggpubr::ggarrange(plotlist=listOfPlots, nrow=6, ncol=4))

# Investigate how topological features relate to the performance of TI methods

metrics <- c("correlation", "him", "featureimp_wcor", "F1_branches")
metric <- metrics[1]
resultsTI <- readr::read_rds(dynbenchmark::result_file("benchmark_results_normalised.rds", experiment_id="06-benchmark"))$data_transformation
performanceTI <- as.matrix(resultsTI[resultsTI[,"metric"]==metric, c("dataset_id", "mean")])
performanceTI <- setNames(as.numeric(performanceTI[,"mean"]), performanceTI[,"dataset_id"])[dynverseIDs]

minPerformance <- min(performanceTI)
maxPerformance <- max(performanceTI)
minMDS2 <- min(fitBN[,2])
maxMDS2 <- max(fitBN[,2])
ggplot(fitBN, aes(x=MDS1, y=MDS2)) +
  geom_smooth(data=data.frame(x=fitBN[,1], y=scales::rescale(performanceTI, to=c(minMDS2, maxMDS2))), aes(x=x, y=y), method="loess", col="red") +
  geom_point(size=3, aes(shape=trajectory), col="white") +
  geom_point(size=3, aes(col=performanceTI, shape=trajectory)) +
  geom_point(data=fitBN[IDsOfInterest,], size=5, aes(shape=trajectory), col="black") +
  geom_point(data=fitBN[IDsOfInterest,], size=3, aes(col=performanceTI[IDsOfInterest], shape=trajectory)) +
  scale_y_continuous(sec.axis=sec_axis(~(maxPerformance - minPerformance) / (maxMDS2 - minMDS2) * (. - maxMDS2) + maxPerformance, name="TI Performance")) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  theme(plot.title=element_text(hjust=0.5), legend.position="top") +
  labs(col=paste0("TI Performance (", metric, ")"))
