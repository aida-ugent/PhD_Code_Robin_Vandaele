# set working directory

setwd("/home/robin/Documents/LeaveEstimator")

# Load libraries

library("ggplot2") # plotting
library("igraph") # working with graphs

# Source functions

for(f in list.files("R")) source(file.path("R", f))

# Load data and labels

maxNumberOfCells <- 7500 # upper bound on number of cells for signatures to be extracted
dimredType <- "dm_diffusionmap" # type of dimensionality reduction
dimredDim <- 2 # number of dimensions of the embedding
graphType <- "10NN" # type of proximity graph, kNN or MST

signaturesFolder <- file.path("Signatures", paste0(dimredType, dimredDim))
if(!file.exists(signaturesFolder)) dir.create(signaturesFolder)
signaturesFolder <- file.path(signaturesFolder, graphType)
if(!file.exists(signaturesFolder)) dir.create(signaturesFolder)

graphsFolder <- file.path("Graphs", paste0(dimredType, dimredDim), graphType)

old <- Sys.time() # get start time
for(dataID in list.files(graphsFolder)){
  thisFolder <- file.path(signaturesFolder, dataID)
  #if(!file.exists(thisFolder)){
    G <- graph_from_data_frame(read.table(file.path(graphsFolder, dataID, "graph.csv")), directed=FALSE)
    if(length(V(G)) <= maxNumberOfCells){
      #dir.create(thisFolder)
      dgm <- graphPersistence(S=G, f=degreeOfCentrality(G))
      #write.table(dgm, file.path(thisFolder, "persistence.csv"))
  #  }
  }
}
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
