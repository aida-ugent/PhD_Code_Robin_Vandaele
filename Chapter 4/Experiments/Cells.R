### Set-up ### 

# Set working directory (change accordingly)

setwd("~/Documents/PhD_code/Chapter 4")

# Load libraries

library("ggplot2")
library("igraph")
library("TDAmapper")

# Source functions

for(file in list.files("R")){
  source(file.path("R", file))
}

# Load and plot data

df.stems <- read.table("Data/Cells/stemcells.csv", header = TRUE, sep = ";")
df.stems.pca <- prcomp(df.stems[,-6], center = TRUE, scale. = TRUE)$x
df.stems.pca <- cbind(data.frame(df.stems.pca), df.stems$label)
colnames(df.stems.pca)[6] <- "label"
ggplot(df.stems.pca, aes(x = PC1, y = PC2, col = label)) + 
  geom_point(size = 5) + coord_fixed() + theme_bw(base_size = 15) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        legend.position = c(0.117, 0.15), legend.key.size = unit(1.5, 'lines')) +
  guides(colour = guide_legend(override.aes = list(size = 5), title = "Cell type")) +
  scale_color_manual(breaks = c("LT-HSC", "ST-HSC", "CMP", "CLP"), 
                     values = c("red", "yellow", "blue", "cyan"))

# Construct Vietoris-Rips graph from data

eps <- 2
D.stems <- dist(df.stems[,-6])
vr.stems <- VRgraph(D.stems, eps)

# Classify and visualize local-global topologies

start_time <- Sys.time()
lg.stems <- classLG(vr.stems, 2)
end_time <- Sys.time()
message(paste("time: ", toString(end_time - start_time), "s", sep = ""))
lg.stems.text <- paste("(", lg.stems$degree, ", ", lg.stems$cycles, ")", sep = "")
ggplot(cbind(df.stems.pca, lg.stems.text), aes(x = PC1, y = PC2, colour = lg.stems.text)) + 
  geom_point(size = 5) + theme_bw(base_size = 15) + 
  coord_fixed() + guides(colour = guide_legend("(Degree, Cycles)")) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        legend.position = c(0.2, 0.1), legend.key.size = unit(1.5, 'lines'))

# Construct underlying graph topology

start_time <- Sys.time()
GT.stems <- constructGT(df.stems, vr.stems, lg.stems, D.stems, 4, method = "mcquitty")
end_time <- Sys.time()
message(paste("time: ", toString(end_time - start_time), "s", sep = ""))
E <- sapply(E(GT.stems), function(e) as.vector(get.edges(GT.stems, e)))
E <- data.frame(matrix(unlist(apply(E, 2, function(e) 
  cbind(df.stems.pca[GT.stems$centers[e[1]], 1:2], 
        df.stems.pca[GT.stems$centers[e[2]], 1:2]))), 
  ncol = 4, byrow = TRUE))
ggplot(cbind(df.stems.pca, lg.stems.text, GT.stems$membership), aes(x = PC1, y = PC2)) + 
  geom_point(size = 5, aes(shape = lg.stems.text, colour = GT.stems$membership)) +
  geom_point(data = df.stems.pca[GT.stems$centers,], col = 'black', size = 4) + 
  theme_bw(base_size = 15) + coord_fixed() + 
  guides(shape = guide_legend("(Degree, Cycles)"), colour = guide_legend("Cluster")) +
  geom_segment(data = E, aes(x = X1, y = X2, xend = X3, yend = X4), size = 2) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        legend.position = "right", legend.key.size = unit(1.5, 'lines'))

# Construct and visualize Mapper graph (filter ~ PCA1)

start_time <- Sys.time()
stems.mapper <- mapper1D(
  distance_matrix = dist(df.stems[,-6]),
  filter_values = df.stems.pca[,1],
  num_intervals = 10,
  percent_overlap = 50,
  num_bins_when_clustering = 10)
end_time <- Sys.time()
message(paste("time: ", toString(end_time - start_time), "s", sep = ""))
plot(graph.adjacency(stems.mapper$adjacency, mode="undirected")

# Vizualize groups on PCA plot

I <- unlist(stems.mapper$points_in_level)
groups <- factor(unlist(sapply(1:length(stems.mapper$points_in_level), function(i){
  rep(i, length(stems.mapper$points_in_level[[i]]))
})))
ggplot(df.stems.pca[I,], aes(x = PC1, y = PC2, col = groups)) + 
  geom_point() + coord_fixed()

# Vizualize branches and bifurcation area on PCA plot

I1 <- unique(unlist(stems.mapper$points_in_vertex[1:8]))
I2 <- stems.mapper$points_in_vertex[[9]]
I3 <- stems.mapper$points_in_vertex[[10]]
I4 <- stems.mapper$points_in_vertex[[11]]
groups <- factor(c(rep(1, length(I1)), rep(2, length(I2)), 
                   rep(3, length(I3)), rep(4, length(I4))))
ggplot(df.stems.pca[c(I1, I2, I3, I4),], aes(x = PC1, y = PC2, col = groups)) + 
  geom_point(size = 5) + coord_fixed() + theme_bw(base_size = 15) +   
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 20),
        legend.position = "right", legend.key.size = unit(1.5, 'lines'))
