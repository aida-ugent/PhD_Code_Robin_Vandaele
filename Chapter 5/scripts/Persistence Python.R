# Import R libraries

devtools::load_all() # source functions
library("ggplot2") # plotting

# Load data

df <- read.table("Data/Toy/Pikachu.csv", header=FALSE, sep=",")
colnames(df) <- c("x", "y")
ggplot(df, aes(x = x, y = y)) +
  geom_point(size = 3) +
  theme_bw() +
  coord_fixed()

# Compute persistent homology from data.frame or matrix containing point cloud observations

path <- "/home/robin/Documents/BCB/Python" # path to folder containing the Python script

Persistencedf <- persistence(df, path=path)

TDA::plot.diagram(Persistencedf$dgms, diagLim=c(0, 25))
TDA::plot.diagram(Persistencedf$dgms, diagLim=c(0, 25), barcode=TRUE)

# Approximate persistence

Persistenceappx <- persistence(df, path=path, n_perm=250)

TDA::plot.diagram(Persistenceappx$dgms, diagLim=c(0, 25))
TDA::plot.diagram(Persistenceappx$dgms, diagLim=c(0, 25), barcode=TRUE)

# Compute persistent homology from distance matrix or object

G <- VRgraPersistence(df, eps=3.5) # connects all points x, y, if ||x-y|| < eps
EG <- get_edges2D(df, G)

ggplot(df[names(V(G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alPersistencea=0.15) +
  geom_point(size=2) +
  theme_bw() +
  coord_fixed()

D <- distances(G)
PersistenceG <- persistence(D, path=path, distance_matrix=TRUE)

TDA::plot.diagram(PersistenceG$dgms, diagLim=c(0, 25))
TDA::plot.diagram(PersistenceG$dgms, diagLim=c(0, 25), barcode=TRUE)

# View birth and death of the most persisting cocycle

cyclesIndex <- which(PersistenceG$dgms[,1] == 1) # diagram rows corresponding to (1D) cycles
persistenceCycles <- PersistenceG$dgms[cyclesIndex, "Death"] - PersistenceG$dgms[cyclesIndex, "Birth"] # persistence of each cycle
maxPersistingCycle <- which.max(persistenceCycles) # max persisting cycle
PersistenceG$dgms[cyclesIndex[maxPersistingCycle], 2:3] # birth and death of max persisting cycle

# Locate most persisting cocycle on the graPersistence (may look weird)

PersistenceG <- persistence(D, path=path, distance_matrix=TRUE, do_cocycles=TRUE)
cyclesIndex <- which(PersistenceG$dgms[,1] == 1)
persistingCycle <- which.max(PersistenceG$dgms[cyclesIndex, "Death"] - PersistenceG$dgms[cyclesIndex, "Birth"])

Ecocycle <- cbind(V(G)$name[PersistenceG$cocycles[[2]][[persistingCycle]][,1]], V(G)$name[PersistenceG$cocycles[[2]][[persistingCycle]][,2]])
Ecocycle <- rbind(Ecocycle[D[Ecocycle] < mean(PersistenceG$dgms[cyclesIndex[maxPersistingCycle], 2:3]),])
Vcocycle <- unique(as.character(Ecocycle))
Ecocycle <- cbind(df[Ecocycle[,1],], df[Ecocycle[,2],])
colnames(Ecocycle) <- c("x1", "y1", "x2", "y2")
ggplot(df[names(V(G)),], aes(x=x, y=y)) +
  geom_point(size=2) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alPersistencea=0.15) +
  geom_segment(data=Ecocycle, aes(x=x1, y=y1, xend=x2, yend=y2), color='red', size=1) +
  geom_point(data=df[Vcocycle,], size=3, col="red") +
  theme_bw() +
  coord_fixed()

Dt <- PersistenceG$dperm2all
Dt[sapply(Dt, is.null)] <- Inf
Dt <- as.matrix(Dt)
str(Dt)
