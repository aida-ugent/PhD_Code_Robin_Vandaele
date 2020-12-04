# Load libraries
devtools::load_all() # load BCB
library("igraphdata") # karate network

# Load and plot data

data(karate)
G <- upgrade_graph(karate)
E(G)$weight <- 1 / E(G)$weight
V(G)$label.color <- rep("black", length(V(G)))
V(G)$label.color[c(1, length(V(G)))] <- "white"
set.seed(5)
layout <- layout.auto(G)
plot(G, layout=layout)

# Conduct backbone pipeline

BCB <- backbone(G)

# View BC, BC-pine, and backbone

color_breaks <- c("blue", "green", "yellow")
V(G)$color <- colorRampPalette(color_breaks)(length(BCB$f))[rank(BCB$f)]
legend_colors <- colorRampPalette(color_breaks)(length(BCB$f))
legend_image <- rev(as.raster(matrix(legend_colors, ncol = 1)))
E(G)$color <- rep("gray", length(E(G)))
E(G)$color[get.edge.ids(G, t(ends(BCB$pine, E(BCB$pine))))] <- "black"
E(G)$color[get.edge.ids(G, t(ends(BCB$B, E(BCB$B))))] <- "red"
edge.width <- rep(1, length(E(G)))
edge.width[get.edge.ids(G, t(ends(BCB$pine, E(BCB$pine))))] <- 2
edge.width[get.edge.ids(G, t(ends(BCB$B, E(BCB$B))))] <- 4
layout(matrix(1:2, ncol=2), width=c(3, 1), height=c(1, 1))
plot(G, layout=layout, edge.width=edge.width, vertex.size=15, vertex.label.cex=1.25)
plot(c(0,2), c(0,1), type="n", axes=FALSE, xlab="", ylab="", main="BC")
text(x=1.5, y = seq(0, 1, l=5), labels=round(seq(min(BCB$f), max(BCB$f), l=5), 2), cex=1.5)
rasterImage(legend_image, 0, 0, 1, 1)
dev.off()

# View backbone

plot(BCB$B, vertex.size=25, vertex.label.cex=1.75, edge.width=10, edge.color='black')

# Quantitative measures backbone

MSB <- medoids_steiner_backbone(G, 2)
FPSB <- furthest_point_backbone(G, 2)
MSTB <- MST_backbone(G, 2)

backbones_evaluate(list(MSB, FPSB, MSTB, BCB$B), G)
