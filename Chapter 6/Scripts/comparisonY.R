# set working directory

setwd("/home/robin/Documents/LeaveEstimator")

# Load libraries

library("ggplot2") # plotting
library("ggpubr") # plotting in grid
library("igraph") # working with graphs

# Source functions

for(f in list.files("R")) source(file.path("R", f))

# Load and plot noisy data

noisyData <- read.csv("Data/flare.csv")
colnames(noisyData) <- c("x", "y")
ggplot(noisyData, aes(x=x, y=y)) +
  geom_point() +
  coord_fixed() +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Persistent homology of X topology using a Vietoris-Rips filtration

D <- dist(noisyData)
maxD <- max(D)
old <- Sys.time()
filtration <- TDA::ripsFiltration(D, maxdimension=0, maxscale=maxD, dist="arbitrary")
diag <- TDA::filtrationDiag(filtration, maxdimension=0, library="Dionysus", printProgress=TRUE)
new <- Sys.time() - old # calculate difference
print(new)
p <-  par(mar=c(3.15, 3.15, 2, 2)); TDA::plot.diagram(diag[["diagram"]], diagLim=c(0,  maxD), 
                                                      main="VR filtration diagram"); p

# Persistent homology of X topology through proximity graph

G <- kNNgraph(noisyData, k=10)
EG <- get_edges2D(noisyData, G)
eccG <- apply(distances(G), 1, max)
ggplot(noisyData[names(V(G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="red", size=0.5, alpha=0.25) +
  geom_point(aes(col=-eccG)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="-Ecc") +
  theme_bw() +
  theme(text=element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_fixed()

old <- Sys.time()
diagG <- graphPersistence(S=G, f=-eccG)
new <- Sys.time() - old # calculate difference
print(new)
p <-  par(mar=c(3.15, 3.15, 2, 2)); TDA::plot.diagram(diagG, diagLim=c(min(-eccG),  max(-eccG)), 
                                                      main="Graph filtration diagram"); p

# Persistent homology of X topology using a functional Vietoris-Rips filtration

f <- eccG[as.character(1:nrow(noisyData))]
f <- (-f - min(-f)) / (max(-f) - min(-f))

old <- Sys.time()
filtration <- TDA::ripsFiltration(D, maxdimension=0, maxscale=maxD, dist="arbitrary")
new2 <- Sys.time() - old # calculate difference

vertices <- which(sapply(filtration$cmplx, function(s) length(s) == 1))
edges <- which(sapply(filtration$cmplx, function(s) length(s) == 2))
filtration$values[vertices] <-  sapply(filtration$cmplx[vertices], function(v) f[v])
filtration$values[edges] <- pmax(filtration$values[edges] / maxD,
                                 sapply(filtration$cmplx[edges], function(e) max(f[e[1]], f[e[2]])))
I <- order(filtration$values)
filtration$cmplx <- filtration$cmplx[I]
filtration$values <- filtration$values[I]

old <- Sys.time()
diag <- TDA::filtrationDiag(filtration, maxdimension=0, library="Dionysus", printProgress=TRUE)
new2 <- new2 + Sys.time() - old # calculate difference and add to previous timing
print(new2)
p <-  par(mar=c(3.15, 3.15, 2, 2)); TDA::plot.diagram(diag[["diagram"]], diagLim=c(min(filtration$values),  max(filtration$values)), 
                                                      main="Functional VR filtration diagram"); p
