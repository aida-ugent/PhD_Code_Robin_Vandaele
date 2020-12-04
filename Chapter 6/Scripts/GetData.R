# Set working directory

setwd("~/Documents/LeaveEstimator")

# Load libraries

library("igraph")

# Extract data

maxNumberOfCells <- 7500 # upper bound on number of cells for data to be extracted
dimredType <- "dm_diffusionmap" # type of dimensionality reduction
ndim <- 5 # maximal number of dimensions after embedding

folder <- file.path("Data", dimredType)
if(!file.exists(folder)) dir.create(folder)

for(id in dynbenchmark::list_datasets()$id){
  dataset <- dynbenchmark::load_datasets(id)
  G <- graph_from_data_frame(dataset$milestone_network, directed=FALSE)
  if(components(G)$no == 1 & length(V(G)) == length(E(G)) + 1){ # only extract data for trees
    thisDataFolder <- file.path(folder, gsub("/", "-", id))
    if(!file.exists(thisDataFolder)) {
      expression <- dataset$expression[[1]]()
      if(nrow(expression) <= maxNumberOfCells){
        
        dimred <- dyndimred::dimred(expression, dimredType, ndim=ndim)
        rownames(dimred) <- rownames(expression)
        
        dir.create(thisDataFolder)
        write.table(dimred, file.path(thisDataFolder, "dimred.csv"))
        write.table(sum(degree(G)==1), file.path(thisDataFolder, "leaves.csv"), row.names=FALSE, col.names=FALSE)   
      
      }
    }
  }
}
