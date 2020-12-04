# Load libraries
devtools::load_all() # load BCB
library("ggplot2") # plotting

# Load and plot data

df.EQ <- rbind(read.table("Data/EQ/query1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE),
               read.table("Data/EQ/query2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE),
               read.table("Data/EQ/query3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE),
               read.table("Data/EQ/query4.csv", header=TRUE, sep=",", stringsAsFactors=FALSE),
               read.table("Data/EQ/query5.csv", header=TRUE, sep=",", stringsAsFactors=FALSE),
               read.table("Data/EQ/query6.csv", header=TRUE, sep=",", stringsAsFactors=FALSE),
               read.table("Data/EQ/query7.csv", header=TRUE, sep=",", stringsAsFactors=FALSE))
ggplot(df.EQ, aes(x=longitude, y=latitude)) +
  geom_point() +
  theme_bw() +
  coord_fixed()

# sample data

set.seed(63)
I <- sample(seq(nrow(df.EQ)), 5000, replace=FALSE)
df <- df.EQ[I, c("longitude", "latitude")]
ggplot(df, aes(x=longitude, y=latitude)) +
  geom_point() +
  theme_bw() +
  coord_fixed()

# Conduct backbone pipeline with custom great circle distances

X <- as.dist(fields::rdist.earth(df, df))
attr(X, "Labels") <- rownames(df)
BCB <- backbone(X, max_leaves=50)

# Move the apparant isolated triangular island from the left to the right for illustration

isle <- which((df.EQ[,"longitude"] < -165) & (abs(df.EQ[,"latitude"] + 25) < 25))
df[I %in% isle, "longitude"] <- df[I %in% isle, "longitude"] + 360
df.EQ[isle, "longitude"] <- df.EQ[isle, "longitude"] + 360

# View backbone

EG <- get_edges2D(df, BCB$G)
Epine <- get_edges2D(df, BCB$pine)
EBCB <- get_edges2D(df, BCB$B)
ggplot(df, aes(x=longitude, y=latitude)) +
  borders("world", col=rgb(0.5, 0.5, 0.5, 0.25), fill=rgb(0.5, 0.5, 0.5, 0.25)) +
  geom_point(data=df.EQ, alpha=0.05) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), col="darkgrey") +
  geom_segment(data=Epine, aes(x=x1, y=y1, xend=x2, yend=y2), col="black") +
  geom_point(size=3, aes(col=BCB$f), alpha=0.3) +
  geom_segment(data=EBCB, aes(x=x1, y=y1, xend=x2, yend=y2), col="red") +
  geom_point(data=df[V(BCB$B)$name,], col="red", size=1.5) +
  scale_colour_gradientn(colors=topo.colors(7)) +
  theme_bw() +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=30, face="bold"),
        legend.title=element_text(size=30), legend.text=element_text(size=25),
        legend.key.height=unit(3, "line")) +
  labs(col="BC") +
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

# Increase number of leaves

BCB <- get_new_leaves(BCB, leaves=35)
EBCB <- get_edges2D(df, BCB$B)

P1 <- ggplot(df, aes(x=longitude, y=latitude)) +
  borders("world", col=rgb(0.5, 0.5, 0.5, 0.25), fill=rgb(0.5, 0.5, 0.5, 0.25)) +
  geom_point(data=df.EQ, alpha=0.05) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), col="darkgrey") +
  geom_point(size=3, col="blue", alpha=0.3) +
  geom_segment(data=EBCB, aes(x=x1, y=y1, xend=x2, yend=y2), col="red") +
  geom_point(data=df[V(BCB$B)$name,], col="red", size=1.5) +
  theme_bw() +
  geom_text(label="Our Method", x=180, y=80, size=5) +
  theme(text=element_text(size=20), axis.text=element_text(size=10), axis.title=element_text(size=15, face="bold")) +
  coord_fixed()

# Compute with k-medoids Steiner tree

MSB <- medoids_steiner_backbone(BCB$G, 35)
EMSB <- get_edges2D(df, MSB)

P2 <- ggplot(df, aes(x=longitude, y=latitude)) +
  borders("world", col=rgb(0.5, 0.5, 0.5, 0.25), fill=rgb(0.5, 0.5, 0.5, 0.25)) +
  geom_point(data=df.EQ, alpha=0.05) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), col="darkgrey") +
  geom_point(size=3, col="blue", alpha=0.3) +
  geom_segment(data=EMSB, aes(x=x1, y=y1, xend=x2, yend=y2), col="red") +
  geom_point(data=df[V(MSB)$name,], col="red", size=1.5) +
  theme_bw() +
  geom_text(label="FacilitySteiner", x=170, y=80, size=5) +
  theme(text=element_text(size=20), axis.text=element_text(size=10), axis.title=element_text(size=15, face="bold")) +
  coord_fixed()

# Compute furthest point sample backbone

FPSB <- furthest_point_backbone(BCB$G, 35)
EFPSB <- get_edges2D(df, FPSB)

P3 <- ggplot(df, aes(x=longitude, y=latitude)) +
  borders("world", col=rgb(0.5, 0.5, 0.5, 0.25), fill=rgb(0.5, 0.5, 0.5, 0.25)) +
  geom_point(data=df.EQ, alpha=0.05) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), col="darkgrey") +
  geom_point(size=3, col="blue", alpha=0.3) +
  geom_segment(data=EFPSB, aes(x=x1, y=y1, xend=x2, yend=y2), col="red") +
  geom_point(data=df[V(FPSB)$name,], col="red", size=1.5) +
  theme_bw() +
  geom_text(label="FarthestPoint", x=170, y=80, size=5) +
  theme(text=element_text(size=20), axis.text=element_text(size=10), axis.title=element_text(size=15, face="bold")) +
  coord_fixed()

# Compute backbone mined from MST

MSTB <- MST_backbone(BCB$G, 35)
EMSTB <- get_edges2D(df, MSTB)

P4 <- ggplot(df, aes(x=longitude, y=latitude)) +
  borders("world", col=rgb(0.5, 0.5, 0.5, 0.25), fill=rgb(0.5, 0.5, 0.5, 0.25)) +
  geom_point(data=df.EQ, alpha=0.05) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), col="darkgrey") +
  geom_point(size=3, col="blue", alpha=0.3) +
  geom_segment(data=EMSTB, aes(x=x1, y=y1, xend=x2, yend=y2), col="red") +
  geom_point(data=df[V(MSTB)$name,], col="red", size=1.5) +
  theme_bw() +
  geom_text(label="CLOFinMSF", x=170, y=80, size=5) +
  theme(text=element_text(size=20), axis.text=element_text(size=10), axis.title=element_text(size=15, face="bold")) +
  coord_fixed()

# Qualitatively compare backbones

ggpubr::ggarrange(P2, P3, P4, P1, nrow=2, ncol=2)

# Quantitative measures backbone

backbones_evaluate(list(MSB, FPSB, MSTB, BCB$B), BCB$G)

# Check for cycles in the backbone

write.table(distances(BCB$G, v=V(BCB$B)$name, to=V(BCB$B)$name), "Server/X.csv") # Save distance matrix for persistence on server
BCB$diagram <- read.table("Server/diagram.csv") # load diagram
op <- par(mar = c(3.25, 3.25, 1, 1))
TDA::plot.diagram(BCB$diagram, diagLim=c(0, 10000))
legend(8000, 5000, legend=c("H0", "H1"), col=c("black", "red"), pch=c(19, 2), pt.lwd=2, box.lty=0)
par(op)

BCB$repCycles <- lapply(jsonlite::fromJSON("Server/cycleLocation.json")[BCB$diagram[,"dimension"]==1], function(CL){
  return(cbind(V(BCB$B)$name[CL[,1]], V(BCB$B)$name[CL[,2]]))
})
maxPersistingCycles <- order((BCB$diagram[,"Death"] - BCB$diagram[,"Birth"])[BCB$diagram[,"dimension"]==1], decreasing=TRUE)[1]
Ecocycles <- unique(do.call("rbind", BCB$repCycles[maxPersistingCycles]))
cocyclesG <- graph_from_data_frame(data.frame(from=Ecocycles[,1], to=Ecocycles[,2], weight=distances(BCB$G)[Ecocycles]), directed=FALSE)
cycles  <- components(cocyclesG)
maxDiamCycle <- which.max(sapply(1:cycles$no, function(c) diameter(induced.subgraph(cocyclesG, which(cycles$membership==c)))))
Ecocycles <- Ecocycles[Ecocycles[,1] %in% V(cocyclesG)$name[cycles$membership==maxDiamCycle],] # representative cycle with maximal diamater
Ecocycles <- cbind(df[Ecocycles[,1],], df[Ecocycles[,2],])
colnames(Ecocycles) <- c("x1", "y1", "x2", "y2")

P1 <- ggplot(df, aes(x=longitude, y=latitude)) +
  borders("world", col=rgb(0.5, 0.5, 0.5, 0.25), fill=rgb(0.5, 0.5, 0.5, 0.25)) +
  geom_point(data=df.EQ, alpha=0.05) +
  geom_segment(data=EG, aes(x=x1, y=y1, xend=x2, yend=y2), col="darkgrey") +
  geom_point(size=3, col="blue", alpha=0.3) +
  geom_segment(data=EBCB, aes(x=x1, y=y1, xend=x2, yend=y2), col="red") +
  geom_point(data=df[V(BCB$B)$name,], col="red", size=1.5) +
  geom_segment(data=Ecocycles, aes(x=x1, y=y1, xend=x2, yend=y2), color="orange", size=1.5, alpha=0.8) +
  theme_bw() +
  geom_text(label="Our Method", x=170, y=80, size=5) +
  theme(text=element_text(size=20), axis.text=element_text(size=10), axis.title=element_text(size=15, face="bold")) +
  coord_fixed()
