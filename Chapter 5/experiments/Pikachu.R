# Load libraries

devtools::load_all() # load BCB
library("ggplot2") # plotting

# Load and plot data

df <- read.table("Data/Toy/Pikachu.csv", header=FALSE, sep=",")
colnames(df) <- c("x", "y")
ggplot(df, aes(x = x, y = y)) +
  geom_point(size = 3) +
  theme_bw() +
  coord_fixed()

# Conduct backbone pipeline

BCB <- backbone(df, eps=3.5, type="Rips")

# View boundary coefficients, pine, and backbone

EG <- get_edges2D(df, BCB$G)
Epine <- get_edges2D(df, BCB$pine)
EBCB <- get_edges2D(df, BCB$B)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.15) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color='black') +
  geom_point(size=2, aes(col=BCB$f)) +
  geom_point(data=df[V(BCB$B)$name,], size=3, col="red") +
  geom_segment(data=EBCB, aes(x=x1, y=y1, xend=x2, yend=y2), color='red', size=2) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="BC") +
  theme_bw() +
  coord_fixed()

# View cost according to the number of leaves

ggplot(BCB$cost, aes(x=leaves, y=cost)) +
  geom_line(aes(group=component, col=component)) +
  geom_point(aes(group=component, col=component)) +
  xlab("number of leaves") +
  ylab("relative cost") +
  theme_bw()

# Increase the total number of leaves, account for global cost only

BCB <- get_new_leaves(BCB, leaves=16)
EBCB <- get_edges2D(df, BCB$B)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.15) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color='black') +
  geom_point(size=2, aes(col=BCB$f)) +
  geom_point(data=df[V(BCB$B)$name,], size=3, col="red") +
  geom_segment(data=EBCB, aes(x=x1, y=y1, xend=x2, yend=y2), color='red', size=2) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="BC") +
  theme_bw() +
  coord_fixed()

# Increase the total number of leaves, standardize cost for each component

BCB <- get_new_leaves(BCB, leaves=16, stdize=TRUE)
EBCB <- get_edges2D(df, BCB$B)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.15) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color='black') +
  geom_point(size=2, aes(col=BCB$f)) +
  geom_point(data=df[V(BCB$B)$name,], size=3, col="red") +
  geom_segment(data=EBCB, aes(x=x1, y=y1, xend=x2, yend=y2), color='red', size=2) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="BC") +
  theme_bw() +
  coord_fixed()

# Increase number of leaves component wise

BCB <- get_new_leaves(BCB, leaves=c(7, 4, 2, 2, 1))
EBCB <- get_edges2D(df, BCB$B)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.15) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color='black') +
  geom_point(size=2, aes(col=BCB$f)) +
  geom_point(data=df[V(BCB$B)$name,], size=3, col='red') +
  geom_segment(data=EBCB, aes(x=x1, y=y1, xend=x2, yend=y2), color='red', size=2) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="BC") +
  theme_bw() +
  coord_fixed()

# Investigate whether there are cycles missing from the backbone

BCB <- check_for_cycles(BCB)
TDA::plot.diagram(BCB$diagram, diagLim=c(0, 20))
legend(17.5, 5, legend=c("H0", "H1"), col=c("black", "red"), pch=c(19, 2), pt.lwd=2, box.lty=0)

Ecocycles <- do.call("rbind", BCB$repCycles)
Ecocycles <- cbind(df[Ecocycles[,1],], df[Ecocycles[,2],])
colnames(Ecocycles) <- c("x1", "y1", "x2", "y2")

ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.15) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color='black') +
  geom_point(size=1, aes(col=BCB$f)) +
  geom_point(data=df[V(BCB$B)$name,], size=2, col='red') +
  geom_segment(data=EBCB, aes(x=x1, y=y1, xend=x2, yend=y2), color='red', size=1) +
  geom_segment(data=Ecocycles, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', size=2, alpha=0.75) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="BC") +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size=20),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Compare persistent homology with that of the original data
# The 'persistence' function allows a more efficient computation in Python when using arbitrary distance matrices
# The 'persistence' function does not allow one to locate representative cycles

PersistenceDF <- persistencePython(X=df) # original point cloud data
TDA::plot.diagram(PersistenceDF$dgms, diagLim=c(0, 20))
legend(17.5, 5, legend=c("H0", "H1"), col=c("black", "red"), pch=c(19, 2), pt.lwd=2, box.lty=0)

PersistenceG <- persistencePython(X=distances(BCB$G), distance_matrix=TRUE) # metric space induced by graph
TDA::plot.diagram(PersistenceG$dgms, diagLim=c(0, 20))
legend(17.5, 5, legend=c("H0", "H1"), col=c("black", "red"), pch=c(19, 2), pt.lwd=2, box.lty=0)

PersistenceSG <- persistencePython(X=distances(BCB$G), distance_matrix=TRUE, n_perm=length(V(BCB$B))) # sample of metric space induced by graph
TDA::plot.diagram(PersistenceSG$dgms, diagLim=c(0, 20))
legend(17.5, 5, legend=c("H0", "H1"), col=c("black", "red"), pch=c(19, 2), pt.lwd=2, box.lty=0)

