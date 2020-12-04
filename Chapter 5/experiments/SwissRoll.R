# Load libraries
devtools::load_all() # load BCB
library("ggplot2") # 2d plotting
library("plot3D") # 3d plotting

# Generate and plot data

n <- 1000
set.seed(420)
df <- IDmining::SwissRoll(n)
df[,"z"] <- (df[,"z"] - min(df["z"])) / (max(df[,"z"]) - min(df[,"z"])) - 1 / 2 # narrow to graph structured backbone
G <- VRgraph(df, eps=0.75)
fit <- data.frame(cmdscale(distances(G)))
colnames(fit) <- c("x", "y")
scatter3D(x=df[V(G)$name, "x"], y=df[V(G)$name, "z"], z=df[V(G)$name, "y"], colvar=fit[,"x"],
          bty="b2", theta=20, phi=10, ylab="z", zlab="y", ticktype="detailed",
          colkey=FALSE, pch=19)

# Conduct backbone pipeline

BCB <- backbone(G)

# View boundary coefficients and backbone

nodes1 <- setdiff(names(V(G)), names(V(BCB$B)))
nodes2 <- names(V(BCB$B)[degree(BCB$B) == 1])
nodes2 <- names(shortest_paths(BCB$B, nodes2[1], nodes2[2])$vpath[[1]])
scatter3D(x=df[nodes1, "x"], y=df[nodes1, "z"], z=df[nodes1, "y"], colvar=fit[nodes1, 1], colkey=FALSE,
          bty="b2", theta=20, phi=10, ylab="z", zlab="y", ticktype="detailed", pch=19, alpha=0.35)
scatter3D(x=df[nodes2, "x"], y=df[nodes2, "z"], z=df[nodes2, "y"], col="black", pch=19, cex=2, add=TRUE)
scatter3D(x=df[nodes2, "x"], y=df[nodes2, "z"], z=df[nodes2, "y"], col="black", add=TRUE, type="l", lwd=5)

EG <- get_edges2D(fit, G)
EBCB <- get_edges2D(fit, BCB$B)
ggplot(fit, aes(x=x, y=y)) +
  geom_point(size=2, aes(col=BCB$f), alpha=0.75) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), size=2, alpha=0.005) +
  geom_segment(data=EBCB, aes(x=x1, y=y1, xend=x2, yend=y2), size=1) +
  geom_point(data=fit[names(V(BCB$B)),], size=3) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  xlab("Isomap1") +
  ylab("Isomap2") +
  labs(col="BC") +
  theme_bw() +
  theme(text = element_text(size=20))

# Compare with MST

MSTB <- MST_backbone(G, 2)
nodes3 <- setdiff(names(V(G)), names(V(MSTB)))
nodes4 <- names(V(MSTB)[degree(MSTB) == 1])
nodes4 <- names(shortest_paths(MSTB, nodes4[1], nodes4[2])$vpath[[1]])
scatter3D(x=df[nodes3, "x"], y=df[nodes3, "z"], z=df[nodes3, "y"], colvar=fit[nodes3, 1], colkey=FALSE,
          bty="b2", theta=20, phi=10, ylab="z", zlab="y", ticktype="detailed", pch=19, alpha=0.35)
scatter3D(x=df[nodes4, "x"], y=df[nodes4, "z"], z=df[nodes4, "y"], col="black", pch=19, cex=1, add=TRUE)
scatter3D(x=df[nodes4, "x"], y=df[nodes4, "z"], z=df[nodes4, "y"], col="black", add=TRUE, type="l", lwd=3)

EMSTB <- get_edges2D(fit, MSTB)
ggplot(fit, aes(x=x, y=y)) +
  geom_point(size=2, aes(col=fit[,1]), alpha=0.75) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), size=2, alpha=0.005) +
  geom_segment(data=EMSTB, aes(x=x1, y=y1, xend=x2, yend=y2), size=1) +
  geom_point(data=fit[names(V(MSTB)),], size=3) +
  scale_colour_gradientn(colours=jet.col(100)) +
  xlab("Isomap1") +
  ylab("Isomap2") +
  labs(col="Isomap1") +
  theme_bw() +
  theme(text = element_text(size=20))

# Compare with longest shortest path

FPSB <- furthest_point_backbone(G, 2)
nodes5 <- setdiff(names(V(G)), names(V(FPSB)))
nodes6 <- names(V(FPSB)[degree(FPSB) == 1])
nodes6 <- names(shortest_paths(FPSB, nodes6[1], nodes6[2])$vpath[[1]])
scatter3D(x=df[nodes5, "x"], y=df[nodes5, "z"], z=df[nodes5, "y"], colvar=fit[nodes5, 1], colkey=FALSE,
          bty="b2", theta=20, phi=10, ylab="z", zlab="y", ticktype="detailed", pch=19, alpha=0.25)
scatter3D(x=df[nodes6, "x"], y=df[nodes6, "z"], z=df[nodes6, "y"], col="black", pch=19, cex=2, add=TRUE)
scatter3D(x=df[nodes6, "x"], y=df[nodes6, "z"], z=df[nodes6, "y"], col="black", add=TRUE, type="l", lwd=5)

EFPSB <- get_edges2D(fit, FPSB)
ggplot(fit, aes(x=x, y=y)) +
  geom_point(size=2, aes(col=fit[,1]), alpha=0.75) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), size=2, alpha=0.005) +
  geom_segment(data=EFPSB, aes(x=x1, y=y1, xend=x2, yend=y2), size=1) +
  geom_point(data=fit[names(V(FPSB)),], size=3) +
  scale_colour_gradientn(colours=jet.col(100)) +
  xlab("Isomap1") +
  ylab("Isomap2") +
  labs(col="Isomap1") +
  theme_bw() +
  theme(text = element_text(size=20))

# Compare with medoids Steiner backbone

MSB <- medoids_steiner_backbone(G, 2)
nodes7 <- setdiff(names(V(G)), names(V(MSB)))
nodes8 <- names(V(MSB)[degree(MSB) == 1])
nodes8 <- names(shortest_paths(MSB, nodes8[1], nodes8[2])$vpath[[1]])
scatter3D(x=df[nodes7, "x"], y=df[nodes7, "z"], z=df[nodes7, "y"], colvar=fit[nodes7, 1], colkey=FALSE,
          bty="b2", theta=20, phi=10, ylab="z", zlab="y", ticktype="detailed", pch=19, alpha=0.25)
scatter3D(x=df[nodes8, "x"], y=df[nodes8, "z"], z=df[nodes8, "y"], col="black", pch=19, cex=2, add=TRUE)
scatter3D(x=df[nodes8, "x"], y=df[nodes8, "z"], z=df[nodes8, "y"], col="black", add=TRUE, type="l", lwd=5)

EMSB <- get_edges2D(fit, MSB)
ggplot(fit, aes(x=x, y=y)) +
  geom_point(size=2, aes(col=fit[,1]), alpha=0.75) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), size=2, alpha=0.005) +
  geom_segment(data=EMSB, aes(x=x1, y=y1, xend=x2, yend=y2), size=1) +
  geom_point(data=fit[names(V(MSB)),], size=3) +
  scale_colour_gradientn(colours=jet.col(100)) +
  xlab("Isomap1") +
  ylab("Isomap2") +
  labs(col="Isomap1") +
  theme_bw() +
  theme(text = element_text(size=20))

# Quantitative measures backbone

backbones_evaluate(list(MSB, FPSB, MSTB, BCB$B), G)
