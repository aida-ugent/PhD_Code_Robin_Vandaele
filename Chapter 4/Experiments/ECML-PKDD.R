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

df.A <- read.table("Data/Toy/ECMLPKDD.csv", sep = ",")
ggplot(df.A, aes(x = V1, y = V2)) + geom_point() + coord_fixed()

# Construct and visualize Vietoris-Rips graph from data

D.A <- dist(df.A)
eps <- 3.5
vr.A <- VRgraph(D.A, eps)
E <- sapply(E(vr.A), function(e) as.vector(get.edges(vr.A, e)))
E <- data.frame(matrix(unlist(apply(E, 2, function(e) 
  cbind(df.A[e[1],], df.A[e[2],]))), ncol = 4, byrow = TRUE))
ggplot(df.A, aes(x = V1, y = V2)) + coord_fixed() +
  geom_segment(data = E, aes(x = X1, y = X2, xend = X3, yend = X4), 
               color = 'red') + geom_point(size = 3, alpha = 0.25)

# Classify and visualize local-global topologies

start_time <- Sys.time()
lg.A <- classLG(vr.A, 3)
end_time <- Sys.time()
message(paste("time: ", toString(end_time - start_time), "s", sep = ""))
lg.A.text <- paste("(", lg.A$degree, ", ", lg.A$cycles, ")", sep = "")
ggplot(cbind(df.A, lg.A.text), aes(x = V1, y = V2, colour = lg.A.text)) +
  geom_segment(data = E, aes(x = X1, y = X2, xend = X3, yend = X4), color = 'black') +
  geom_point(size = 3, alpha = 0.25) + 
  geom_point(size = 3) + theme_bw(base_size = 15) + 
  coord_fixed() + guides(colour = guide_legend("(Degree, Cycles)")) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Construct underlying graph topology

start_time <- Sys.time()
GT.A <- constructGT(df.A, vr.A, lg.A, D.A, 4, method = "complete")
end_time <- Sys.time()
message(paste("time: ", toString(end_time - start_time), "s", sep = ""))
E <- sapply(E(GT.A), function(e) as.vector(get.edges(GT.A, e)))
E <- data.frame(matrix(unlist(apply(E, 2, function(e) 
  cbind(df.A[GT.A$centers[e[1]],], df.A[GT.A$centers[e[2]],]))), 
  ncol = 4, byrow = TRUE))
ggplot(cbind(df.A, lg.A.text, GT.A$membership), aes(x = V1, y = V2)) + 
  geom_point(size = 3, aes(shape = lg.A.text, colour = GT.A$membership)) +
  scale_shape_manual(values = seq(0, 7)) +
  geom_point(data = df.A[GT.A$centers,], col = 'black', size = 3) + 
  theme_bw(base_size = 15) + coord_fixed() + 
  guides(shape = guide_legend("(Degree, Cycles)"), 
         colour = guide_legend("Cluster")) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  geom_segment(data = E, aes(x = X1, y = X2, xend = X3, yend = X4))
