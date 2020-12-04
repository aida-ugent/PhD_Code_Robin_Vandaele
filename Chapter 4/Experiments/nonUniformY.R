### Set-up ### 

# Set working directory (change accordingly)

setwd("~/Documents/PhD_code/Chapter 4")

# Load libraries

library("ggplot2")
library("igraph")

# Source functions

for(file in list.files("R")){
  source(file.path("R", file))
}

# Load and plot data

df.Y <- read.table("Data/Toy/ToyY.csv", header = FALSE, sep = ",")
colnames(df.Y) <- c("x", "y")
ggplot(df.Y, aes(x = x, y = y)) + geom_point() + coord_fixed()

# Construct and visualize Vietoris-Rips graph from data

D.Y <- dist(df.Y)
eps <- 15
vr.Y <- VRgraph(D.Y, eps)
E <- sapply(E(vr.Y), function(e) as.vector(get.edges(vr.Y, e)))
E <- data.frame(matrix(unlist(apply(E, 2, function(e) 
  cbind(df.Y[e[1],], df.Y[e[2],]))), ncol = 4, byrow = TRUE))
ggplot(df.Y, aes(x = x, y = y)) + coord_fixed() +
  geom_segment(data = E, aes(x = X1, y = X2, xend = X3, yend = X4), 
               color = 'red') + geom_point(size = 3, alpha = 0.25)

# Classify and visualize local-global topologies

start_time <- Sys.time()
lg.Y <- classLG(vr.Y, 3)
end_time <- Sys.time()
message(paste("time: ", toString(end_time - start_time), "s", sep = ""))
lg.Y.text <- paste("(", lg.Y$degree, ", ", lg.Y$cycles, ")", sep = "")
ggplot(cbind(df.Y, lg.Y.text), aes(x = x, y = y, colour = lg.Y.text)) + 
  geom_point(size = 3) + theme_bw(base_size = 15) + 
  coord_fixed() + guides(colour = guide_legend("(Degree, Cycles)")) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Construct underlying graph topology

start_time <- Sys.time()
GT.Y <- constructGT(df.Y, vr.Y, lg.Y, D.Y, 4, method = "complete")
end_time <- Sys.time()
message(paste("time: ", toString(end_time - start_time), "s", sep = ""))
E <- sapply(E(GT.Y), function(e) as.vector(get.edges(GT.Y, e)))
E <- data.frame(matrix(unlist(apply(E, 2, function(e) 
  cbind(df.Y[GT.Y$centers[e[1]],], df.Y[GT.Y$centers[e[2]],]))), 
  ncol = 4, byrow = TRUE))
ggplot(cbind(df.Y, lg.Y.text, GT.Y$membership), aes(x = x, y = y)) + 
  geom_point(size = 3, aes(shape = lg.Y.text, colour = GT.Y$membership)) +
  geom_point(data = df.Y[GT.Y$centers,], size = 5) + 
  theme_bw(base_size = 15) + coord_fixed() + 
  guides(shape = guide_legend("(Degree, Cycles)"), colour = guide_legend("Cluster")) +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  geom_segment(data = E, aes(x = X1, y = X2, xend = X3, yend = X4), size = 3, alpha = 1)
