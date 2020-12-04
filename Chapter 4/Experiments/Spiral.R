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

df.S <- read.table("Data/Toy/spiral.csv", sep = ",")
ggplot(df.S, aes(x = V1, y = V2)) + geom_point() + coord_fixed()

# Construct and visualize Vietoris-Rips graph from data

D.S <- dist(df.S)
eps <- 3
vr.S <- VRgraph(D.S, eps)
E <- sapply(E(vr.S), function(e) as.vector(get.edges(vr.S, e)))
E <- data.frame(matrix(unlist(apply(E, 2, function(e) 
  cbind(df.S[e[1],], df.S[e[2],]))), ncol = 4, byrow = TRUE))
ggplot(df.S, aes(x = V1, y = V2)) + coord_fixed() +
  geom_segment(data = E, aes(x = X1, y = X2, xend = X3, yend = X4), 
               color = 'red') + geom_point(size = 3, alpha = 0.25)

# Classify and visualize local-global topologies

start_time <- Sys.time()
lg.S <- classLG(vr.S, 3)
end_time <- Sys.time()
message(paste("time: ", toString(end_time - start_time), "s", sep = ""))
lg.S.text <- paste("(", lg.S$degree, ", ", lg.S$cycles, ")", sep = "")
ggplot(cbind(df.S, lg.S.text), aes(x = V1, y = V2, colour = lg.S.text)) + 
  geom_point(size = 3) + theme_bw(base_size = 15) + 
  coord_fixed() + guides(colour = guide_legend("(Degree, Cycles)")) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Construct underlying graph topology

start_time <- Sys.time()
GT.S <- constructGT(df.S, vr.S, lg.S, D.S, 3, method = "complete")
end_time <- Sys.time()
message(paste("time: ", toString(end_time - start_time), "s", sep = ""))
E <- sapply(E(GT.S), function(e) as.vector(get.edges(GT.S, e)))
E <- data.frame(matrix(unlist(apply(E, 2, function(e) 
  cbind(df.S[GT.S$centers[e[1]],], df.S[GT.S$centers[e[2]],]))), 
  ncol = 4, byrow = TRUE))
ggplot(cbind(df.S, lg.S.text, GT.S$membership), aes(x = V1, y = V2)) + 
  geom_point(size = 3, aes(shape = lg.S.text, colour = GT.S$membership)) +
  geom_point(data = df.S[GT.S$centers,], col = 'black', size = 3) + 
  theme_bw(base_size = 15) + coord_fixed() + 
  guides(shape = guide_legend("(Degree, Cycles)"), colour = guide_legend("Cluster")) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  geom_segment(data = E, aes(x = X1, y = X2, xend = X3, yend = X4))

# Divide long branch into consecutive components

GT.S <- splitBranch(GT.S, vr.S, D.S, 3, 30)
E <- sapply(E(GT.S), function(e) as.vector(get.edges(GT.S, e)))
E <- data.frame(matrix(unlist(apply(E, 2, function(e) 
  cbind(df.S[GT.S$centers[e[1]], 1:2], df.S[GT.S$centers[e[2]], 1:2]))), 
  ncol = 4, byrow = TRUE))
ggplot(cbind(df.S, lg.S.text, GT.S$membership), aes(x = V1, y = V2)) + 
  geom_point(size = 3, aes(shape = lg.S.text, colour = GT.S$membership)) +
  geom_point(data = df.S[GT.S$centers,], size = 4) + 
  geom_segment(data = E, aes(x = X1, y = X2, xend = X3, yend = X4), size = 2, alpha = 0.5) +
  theme_bw(base_size = 15) + coord_fixed() +
  guides(shape = guide_legend("(Degree, Cycles)"), colour = guide_legend("Cluster")) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
