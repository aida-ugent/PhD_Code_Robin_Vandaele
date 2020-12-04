# Load libraries
devtools::load_all() # load BCB
library("ggplot2") # plotting
library("randomcoloR") # branch assignment

# Load and plot data

df <- read.table("Data/Toy/curvX.csv", header=FALSE, sep=",")
colnames(df) <- c("x", "y")
ggplot(df, aes(x=x, y=y)) +
  geom_point(size=3) +
  theme_bw() +
  coord_fixed()

# Conduct backbone pipeline

BCB <- backbone(df, eps=8, type="Rips", assign=TRUE)

# View boundary coefficients

EG <- get_edges2D(df, BCB$G)
ggplot(df[names(V(BCB$G)),], aes(x=x, y=y)) +
  geom_segment(data=EG, aes(x = x1, y = y1, xend = x2, yend = y2), color='black', alpha=0.15) +
  geom_point(size=2, aes(col=BCB$f)) +
  scale_colour_gradientn(colours=topo.colors(7)) +
  labs(col="BC") +
  theme_bw() +
  coord_fixed()

# View cost according to the number of leaves

ggplot(BCB$cost, aes(x=leaves, y=cost)) +
  geom_line(aes(group=component, col=component), size=1.5) +
  geom_vline(xintercept=length(which(degree(BCB$B) == 1)), linetype="dashed", size=1.5) +
  geom_point(aes(group=component, col=component, size=1.5)) +
  xlab("number of leaves") +
  ylab("relative cost") +
  theme_bw() +
  theme(text = element_text(size=20), legend.position="none")

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
  theme_bw() +
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15)) +
  theme(text = element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
