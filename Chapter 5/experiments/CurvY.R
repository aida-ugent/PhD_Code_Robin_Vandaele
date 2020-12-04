# Load libraries
devtools::load_all() # load BCB
library("ggplot2") # plotting

# Load and plot data

df <- read.table("Data/Toy/curvY.csv", header=FALSE, sep=",")
colnames(df) <- c("x", "y")
ggplot(df, aes(x=x, y=y)) +
  geom_point(size=3) +
  theme_bw() +
  coord_fixed()

# Conduct backbone pipeline

BCB <- backbone(df, eps=10, type="Rips", assign=TRUE)

# View boundary coefficients

EG <- get_edges2D(df, BCB$G)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.15) +
  geom_point(size=2, aes(col=BCB$f)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="BC") +
  theme_bw() +
  coord_fixed() +
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# View the BC pine

Epine <- get_edges2D(df, BCB$pine)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x = x1, y = y1, xend = x2, yend = y2), color='black', alpha=0.15) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color="red", size=1) +
  geom_point(size=1, aes(col=BCB$f)) +
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

# View backbone and branch assignments

Epine <- get_edges2D(df, BCB$pine)
EBCB <- get_edges2D(df, BCB$B)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.5, size=1) +
  geom_point(size=3, fill=BCB$col, alpha=0.5, pch=21) +
  geom_segment(data=EBCB, aes(x= x1, y=y1, xend=x2, yend=y2), col="black", size=2.75) +
  geom_segment(data=cbind(EBCB, as.factor(BCB$branch)),
               aes(x=x1, y=y1, xend=x2, yend=y2, col=as.factor(BCB$branch)), size=2) +
  geom_point(data=df[V(BCB$B)$name,], fill=BCB$col[V(BCB$B)$name], size=5, pch=21) +
  scale_color_manual(values=BCB$palette) +
  guides(col = guide_legend(title="Path")) +
  coord_fixed() +
  theme_bw()

# View cost according to the number of leaves

ggplot(BCB$cost, aes(x=leaves, y=cost)) +
  geom_line(aes(group=component, col=component)) +
  geom_point(aes(group=component, col=component)) +
  xlab("number of leaves") +
  ylab("relative cost") +
  theme_bw()

# Increase number of leaves

BCB <- get_new_leaves(BCB, leaves=3)
EBCB <- get_edges2D(df, BCB$B)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.5, size=1) +
  geom_point(size=3, fill=BCB$col, alpha=0.5, pch=21) +
  geom_segment(data=EBCB, aes(x= x1, y=y1, xend=x2, yend=y2), col="black", size=2.75) +
  geom_segment(data=cbind(EBCB, as.factor(BCB$branch)),
               aes(x=x1, y=y1, xend=x2, yend=y2, col=as.factor(BCB$branch)), size=2) +
  geom_point(data=df[V(BCB$B)$name,], fill=BCB$col[V(BCB$B)$name], size=5, pch=21) +
  scale_color_manual(values=BCB$palette) +
  guides(col = guide_legend(title="Path")) +
  coord_fixed() +
  theme_bw()

# View the backbone and betweenness of the pine

Betw <- betweenness(BCB$pine)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", size=1) +
  geom_point(size=1, aes(col=Betw)) +
  geom_segment(data=EBCB, aes(x= x1, y=y1, xend=x2, yend=y2), col="red", size=2.75) +
  geom_point(data=df[names(V(BCB$B)),], size=3, aes(col=Betw[names(V(BCB$B))])) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="Betw") +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Qualitatively compare various backbones

ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) + # backbone mined from BC-pine
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.5, size=1) +
  geom_point(size=2, aes(col=BCB$f)) +
  geom_segment(data=EBCB, aes(x=x1, y=y1, xend=x2, yend=y2), size=1.5, col="red") +
  geom_point(data=df[V(BCB$B)$name,], col="red", size=4) +
  geom_point(data=df[V(BCB$B)$name[degree(BCB$B)==1],], size=6) +
  geom_point(data=df[V(BCB$B)$name[degree(BCB$B)==1],], col="orange", size=4) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="BC") +
  coord_fixed() +
  theme_bw()

MSB <- medoids_steiner_backbone(BCB$G, 3) # k-median steiner tree
EMSB <- get_edges2D(df, MSB)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_point(size=1) +
  geom_segment(data=EMSB, aes(x=x1, y=y1, xend=x2, yend=y2), size=1.5, col="red") +
  geom_point(data=df[V(MSB)$name,], col="red", size=4) +
  geom_point(data=df[MSB$terminals,], size=6) +
  geom_point(data=df[MSB$terminals,], col="orange", size=4) +
  coord_fixed() +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

FPSB <- furthest_point_backbone(BCB$G, 3) # furthest point sample
EFPSB <- get_edges2D(df, FPSB)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_point(size=1) +
  geom_segment(data=EFPSB, aes(x=x1, y=y1, xend=x2, yend=y2), size=1.5, col="red") +
  geom_point(data=df[V(FPSB)$name,], col="red", size=4) +
  geom_point(data=df[V(FPSB)$name[degree(FPSB)==1],], size=6) +
  geom_point(data=df[V(FPSB)$name[degree(FPSB)==1],], col="orange", size=4) +
  coord_fixed() +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

MSTB <- MST_backbone(BCB$G, 3) # backbone mined through MST
EMSTB <- get_edges2D(df, MSTB)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", alpha=0.1) +
  geom_point(size=1) +
  geom_segment(data=EMSTB, aes(x=x1, y=y1, xend=x2, yend=y2), size=1, col="red") +
  geom_point(data=df[V(MSTB)$name,], col="red", size=2) +
  geom_point(data=df[V(MSTB)$name[degree(MSTB)==1],], size=6) +
  geom_point(data=df[V(MSTB)$name[degree(MSTB)==1],], col="orange", size=4) +
  coord_fixed() +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
