# Load libraries
devtools::load_all() # load BCB
library("ggplot2") # plotting

# Sample and plot data

set.seed(42)
rho <- sqrt(runif(250))
theta <- runif(250, 0, 2 * pi)
df.D <- data.frame(x=rho * cos(theta), y=rho * sin(theta))
ggplot(df.D, aes(x=x, y=y)) +
  geom_point(size=3) +
  coord_fixed() +
  theme_bw()

# Construct and plot Vietoris-Rips graph

vr.D <- VRgraph(df.D, 10)
Evr.D <- get_edges2D(df.D, vr.D)
ggplot(df.D, aes(x=x, y=y)) +
  geom_segment(data = Evr.D, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.01) +
  geom_point(size=3) +
  coord_fixed() +
  theme_bw()

# Compute boundary coefficient

BC <- boundary_coefficient(vr.D)
ggplot(df.D[names(V(vr.D)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.D, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.01) +
  geom_point(size=2, aes(col = BC)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# View the BC pine

BCB <- backbone(vr.D)
Epine <- get_edges2D(df.D, BCB$pine)
ggplot(df.D[names(V(vr.D)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.D, aes(x = x1, y = y1, xend = x2, yend = y2), color='black', alpha=0.01) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color="red", size=1) +
  geom_point(size=1, aes(col=BC)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  theme_bw() +
  coord_fixed() +
  theme(text = element_text(size=20),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Compare to ordinary local cluster coefficient

LCC <- transitivity(vr.D, type="local")
ggplot(df.D[names(V(vr.D)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.D, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.01) +
  geom_point(size=3, aes(col=LCC)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  coord_fixed() +
  theme_bw()

# Compare to local efficiency

Eff <- brainwaver::local.efficiency(as_adjacency_matrix(vr.D), as_adjacency_matrix(vr.D, attr="weight", sparse=FALSE))[[2]]
Eff <- round(Eff, 9)
ggplot(df.D[names(V(vr.D)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.D, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.01) +
  geom_point(size=2, aes(col=Eff)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Compare to Barrat's generalized local cluster coefficient

net.D <- apply(as_edgelist(vr.D), 2, as.numeric)
net.D <- cbind(net.D, E(vr.D)$weight)
colnames(net.D) <- c("i", "j", "w")
net.D <- rbind(net.D, net.D[,c(2, 1, 3)])
B_LCC <- tnet::clustering_local_w(net.D)[,2]
ggplot(df.D[names(V(vr.D)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.D, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.01) +
  geom_point(size=3, aes(col=B_LCC)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  coord_fixed() +
  theme_bw()

# Compare to Onnela's generalized local cluster coefficient

O_LCC <-  DirectedClustering::DirectedClustering::ClustF(as_adjacency_matrix(vr.D, attr = "weight", sparse = FALSE))$LocalCC
ggplot(df.D[names(V(vr.D)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.D, aes(x = x1, y = y1, xend = x2, yend = y2), color='black', alpha=0.01) +
  geom_point(size=3, aes(col=O_LCC)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  coord_fixed() +
  theme_bw()

# Onnela's local cluster coefficient performs very similar when outliers are absent
# However, in the presence of outliers, it is outperformed by the boundary coefficient

df.C <- read.table("Data/Toy/curvC.csv", header=FALSE, sep=",")
colnames(df.C) <- c("x", "y")
ggplot(df.C, aes(x=x, y=y)) +
  geom_point(size=3) +
  coord_fixed() +
  theme_bw()

vr.C <- VRgraph(df.C, 10)
Evr.C <- get_edges2D(df.C, vr.C)

BC <- boundary_coefficient(vr.C)
ggplot(df.C[names(V(vr.C)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.C, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.2) +
  geom_point(size=3, aes(col=BC)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  coord_fixed() +
  theme_bw()

O_LCC <-  DirectedClustering::ClustF(as_adjacency_matrix(vr.C, attr = "weight", sparse = FALSE))$LocalCC
ggplot(df.C[names(V(vr.C)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.C, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.2) +
  geom_point(size=2, aes(col=O_LCC)) +
  scale_colour_gradientn(colours = topo.colors(7)) +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Compare to Zhangs's generalized local cluster coefficient

Z_LCC <- qgraph::clustZhang(vr.D)[,1]
ggplot(df.D[names(V(vr.D)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.D, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.01) +
  geom_point(size=3, aes(col=Z_LCC)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  coord_fixed() +
  theme_bw()

# Compare to eccentricity

Ecc <- eccentricity(vr.D)
ggplot(df.D[names(V(vr.D)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.D, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.01) +
  geom_point(size=3, aes(col=Ecc)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  coord_fixed() +
  theme_bw()

# Compare to betweenness

Betw <- betweenness(vr.D)
ggplot(df.D[names(V(vr.D)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.D, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.01) +
  geom_point(size=3, aes(col=Betw)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  coord_fixed() +
  theme_bw()

# Let's see how betweenness performs with curving structures

df.Y <- distinct(read.table("Data/Toy/curvY.csv", header=FALSE, sep=",")) # we cannot compute betweenness in the presence of duplicate rows
colnames(df.Y) <- c("x", "y")
ggplot(df.Y, aes(x=x, y=y)) +
  geom_point(size=3) +
  coord_fixed() +
  theme_bw()

vr.Y <- VRgraph(df.Y, 10)
Evr.Y <- get_edges2D(df.Y, vr.Y)

BC <- boundary_coefficient(vr.Y)
ggplot(df.Y[names(V(vr.Y)),], aes(x=x, y=y)) +
  geom_segment(data=Evr.Y, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.2) +
  geom_point(size=3, aes(col=BC)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  coord_fixed() +
  theme_bw()

Betw <- betweenness(vr.Y)
ggplot(df.Y[names(V(vr.Y))[order(Betw)],], aes(x=x, y=y)) +
  geom_segment(data=Evr.Y, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.2) +
  geom_point(size=2, aes(col=Betw[order(Betw)])) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  coord_fixed() +
  labs(color="Betw") +
  theme_bw() +
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
