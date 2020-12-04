# Load libraries
devtools::load_all() # load BCB
library("ggplot2") # plotting

# Load and plot data

cells_info <- readRDS("Data/Cell/Synthetic/Cell_Data_Synth")
fit <- cbind(data.frame(dyndimred::dimred(cells_info$expression, method="mds", ndim=2)),
             cells_info$prior_information$groups_id$group_id[gtools::mixedorder(cells_info$prior_information$groups_id$cell_id)])
colnames(fit)[1:2] <- c("x", "y")
colnames(fit)[ncol(fit)] <- "group"
fit$group <- factor(fit$group)
cols <- c("1"=rgb(250, 127, 113, maxColorValue=255),
          "2"=rgb(252, 179, 97, maxColorValue=255))
G <- graph_from_data_frame(cells_info$milestone_network)
V(G)$color <- cols
plot(G, layout=layout_as_tree, edge.arrow.size=0.25, edge.color="black")
ggplot(fit, aes(x = x, y = y, fill = group)) +
  geom_point(size=3, col="black", pch=21) +
  scale_colour_manual(values=cols, aesthetics="fill") +
  coord_fixed() +
  theme_bw()

# Conduct backbone pipeline

BCB <- backbone(cells_info$expression)

# View backbone in original data

EG <- get_edges2D(fit[,1:2], BCB$G)
Epine <- get_edges2D(fit[,1:2], BCB$pine)
EBCB <- get_edges2D(fit[,1:2], BCB$B)
ggplot(fit[V(BCB$G)$name,], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.75, size=1) +
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

set.seed(12)
layout_pine <- layout_with_fr(BCB$pine)
edge.color <- rep(rgb(0.5, 0.5, 0.5), length(E(BCB$pine)))
edge.color[get.edge.ids(BCB$pine, as.vector(t(ends(BCB$B, E(BCB$B)))))] <- "blue"
edge.width <- rep(1, length(E(BCB$pine)))
edge.width[get.edge.ids(BCB$pine, as.vector(t(ends(BCB$B, E(BCB$B)))))] <- 2
plot(BCB$pine, vertex.size=2, vertex.label=NA, vertex.color=cols[fit[V(BCB$pine)$name, 3]], vertex.frame.color="black",
     edge.color=edge.color, edge.width=edge.width, layout=layout_pine, asp=0.5)

MSTcells <- mst(BCB$G)
set.seed(4)
layout_MST <- layout_with_fr(MSTcells)
plot(MSTcells, vertex.size=2, vertex.label=NA, vertex.color=cols[fit[V(MSTcells)$name, 3]], vertex.frame.color="black",
     edge.color=rgb(0.5, 0.5, 0.5), edge.width=1, layout=layout_MST, asp=0.5)

# Quantitative measures backbone

MSB <- medoids_steiner_backbone(BCB$G, 2) # k-medoids Steiner Tree
FPSB <- furthest_point_backbone(BCB$G, 2) # backbone through furthest point sample
MSTB <- MST_backbone(BCB$G, 2) # backbone through regular MST

backbones_evaluate(list(MSB, FPSB, MSTB, BCB$B), BCB$G)

# Investigate other backbones

EMSB <- get_edges2D(fit[,1:2], MSB) # k-medoids Steiner Tree
ggplot(fit[V(BCB$G)$name,], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_point(size=3, pch=21, aes(fill=group)) +
  geom_segment(data=EMSB, aes(x= x1, y=y1, xend=x2, yend=y2), col="black", size=2.75) +
  geom_segment(data=EMSB, aes(x=x1, y=y1, xend=x2, yend=y2), col="blue", size=2, alpha=0.75) +
  geom_point(data=fit[V(MSB)$name,], fill="blue", size=5, pch=21) +
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

EFPSB <- get_edges2D(fit[,1:2], FPSB) # backbone through furthest point sample
ggplot(fit[V(BCB$G)$name,], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_point(size=3, pch=21, aes(fill=group)) +
  geom_segment(data=EFPSB, aes(x= x1, y=y1, xend=x2, yend=y2), col="black", size=2.75) +
  geom_segment(data=EFPSB, aes(x=x1, y=y1, xend=x2, yend=y2), col="blue", size=2, alpha=0.75) +
  geom_point(data=fit[V(FPSB)$name,], fill="blue", size=5, pch=21) +
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

EMSTB <- get_edges2D(fit[,1:2], MSTB) # backbone through regular MST
ggplot(fit[V(BCB$G)$name,], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_point(size=3, pch=21, aes(fill=group)) +
  geom_segment(data=EMSTB, aes(x= x1, y=y1, xend=x2, yend=y2), col="black", size=2.75) +
  geom_segment(data=EMSTB, aes(x=x1, y=y1, xend=x2, yend=y2), col="blue", size=2, alpha=0.75) +
  geom_point(data=fit[V(MSTB)$name,], fill="blue", size=5, pch=21) +
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
