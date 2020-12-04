# Load libraries
devtools::load_all() # load BCB
library("ggplot2") # plotting

# Load and plot data

cells_info <- readRDS("Data/Cell/Real/cellbench-SC1_luyitian")
fit <- cbind(data.frame(dyndimred::dimred(cells_info$expression, method="dm_diffusionmap", ndim=5)), cells_info$grouping)
colnames(fit)[1:2] <- c("x", "y")
colnames(fit)[ncol(fit)] <- "group"
fit$group <- factor(fit$group)
fit[,1] <- -fit[,1]
cols <- c("H1975"=rgb(250, 127, 113, maxColorValue = 255),
          "HCC827"=rgb(178, 221, 104, maxColorValue = 255),
          "H2228"=rgb(254, 254, 178, maxColorValue = 255),
          "H1975,H2228,HCC827"=rgb(252, 179, 97, maxColorValue = 255))
G <- graph_from_data_frame(cells_info$milestone_network)
V(G)$color <- cols
op <- par(mar = c(0, 0, 0, 0))
plot(G, vertex.shape="rectangle", vertex.label.cex=1, vertex.size=c(75, 75, 75, 150),
     edge.arrow.size=0.25, edge.color="black")
par(op)
ggplot(fit, aes(x=x, y=y, fill=group)) +
  geom_point(size=3, col="black", pch=21) +
  scale_colour_manual(values=cols, aesthetics="fill") +
  coord_fixed() +
  theme_bw() +
  theme(legend.title=element_text(size=15, face="bold"), legend.text=element_text(size = 12))

# Conduct backbone pipeline on original data

BCB <- backbone(cells_info$expression)

# View backbone

EG <- get_edges2D(fit[,1:2], BCB$G)
Epine <- get_edges2D(fit[,1:2], BCB$pine)
ggplot(fit[V(BCB$G)$name,], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color="red", alpha=.25, size=1) +
  geom_point(size=3, aes(fill=group), pch=21) +
  geom_point(data=fit[V(BCB$B)$name,], fill="blue", size=5, pch=21) +
  scale_colour_manual(values=cols, aesthetics="fill") +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Conduct backbone pipeline on a dimensionality reduction of the data

BCB <- backbone(fit[,1:3])
BCB2 <- get_new_leaves(BCB, sum(degree(G) == 1))

# View backbone

EG <- get_edges2D(fit[,1:2], BCB$G)
Epine <- get_edges2D(fit[,1:2], BCB$pine)
EBCB <- get_edges2D(fit[,1:2], BCB$B)
EBCB2 <- get_edges2D(fit[,1:2], BCB2$B)
ggplot(fit[V(BCB$G)$name,], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color="red", alpha=.75, size=1) +
  geom_point(size=3, pch=21, aes(fill=group)) +
  geom_segment(data=EBCB, aes(x= x1, y=y1, xend=x2, yend=y2), col="black", size=2.75) +
  geom_segment(data=EBCB, aes(x=x1, y=y1, xend=x2, yend=y2), col="blue", size=2, alpha=0.75) +
  geom_point(data=fit[V(BCB$B)$name,], fill="blue", size=5, pch=21) +
  scale_colour_manual(values=cols, aesthetics="fill") +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Visualize BC-pine, and compare with regular MST

set.seed(20)
layout_pine <- layout_with_fr(BCB$pine)
edge.color <- rep(rgb(0.5, 0.5, 0.5), length(E(BCB$pine)))
edge.color[get.edge.ids(BCB$pine, as.vector(t(ends(BCB$B, E(BCB$B)))))] <- "blue"
edge.width <- rep(1, length(E(BCB$pine)))
edge.width[get.edge.ids(BCB$pine, as.vector(t(ends(BCB$B, E(BCB$B)))))] <- 2
plot(BCB$pine, vertex.size=2, vertex.label=NA, vertex.color=cols[fit[V(BCB$pine)$name, 6]], vertex.frame.color="black",
     edge.color=edge.color, edge.width=edge.width, layout=layout_pine, asp=0.5)

MSTcells <- mst(BCB$G)
set.seed(97)
layout_MST <- layout_with_fr(MSTcells)
plot(MSTcells, vertex.size=2, vertex.label=NA, vertex.color=cols[fit[V(MSTcells)$name, 6]], vertex.frame.color="black",
     edge.color=rgb(0.5, 0.5, 0.5), edge.width=1, layout=layout_MST, asp=0.5)

# Quantitative measures backbone

MSB <- medoids_steiner_backbone(BCB$G, 3) # k-medoids Steiner Tree
FPSB <- furthest_point_backbone(BCB$G, 3) # backbone through furthest point sample
MSTB <- MST_backbone(BCB$G, 3) # backbone through regular MST

backbones_evaluate(list(MSB, FPSB, MSTB, BCB$B), BCB$G)
