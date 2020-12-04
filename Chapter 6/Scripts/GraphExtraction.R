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

folder <- file.path("Graphs", paste0(dimredType, dimredDim))
if(!file.exists(folder)) dir.create(folder)
folder <- file.path(folder, graphType)
if(!file.exists(folder)) dir.create(folder)

for(dataID in list.files(file.path("Data", dimredType))){
  thisFolder <- file.path(folder, dataID)
  if(!file.exists(thisFolder)){
    if(graphType=="MST"){
      G <- emstreeR::ComputeMST(read.table(file.path("Data", dimredType, dataID, "dimred.csv"))[,1:dimredDim], verbose=FALSE)[,c("from", "to", "distance")]
      colnames(G)[3] <- "weight"
      G <- graph_from_data_frame(G, directed=FALSE)
    }
    else{
      k <- readr::parse_number(graphType)
      G <- kNNgraph(read.table(file.path("Data", dimredType, dataID, "dimred.csv"))[,1:dimredDim], k=k)
    }
    if(components(G)$no == 1){
      dir.create(thisFolder)
      write.table(as_data_frame(G), file.path(thisFolder, "graph.csv"))
    }
  }
}
