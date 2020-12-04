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

df.E <- read.table("Data/Toy/df_Ellipse.csv", header = TRUE, sep = ",")
df.E[,1] <- NULL
ggplot(df.E, aes(x = x, y = y)) + geom_point() + coord_fixed()

# Construct and visualize Vietoris-Rips graph from data

D.E <- dist(df.E)
eps <- 0.3
vr.E <- VRgraph(D.E, eps)
E <- sapply(E(vr.E), function(e) as.vector(get.edges(vr.E, e)))
E <- data.frame(matrix(unlist(apply(E, 2, function(e) 
  cbind(df.E[e[1],], df.E[e[2],]))), ncol = 4, byrow = TRUE))
ggplot(df.E, aes(x = x, y = y)) + coord_fixed() +
  geom_segment(data = E, aes(x = X1, y = X2, xend = X3, yend = X4), 
               color = 'red') + geom_point(size = 3, alpha = 0.25)

# Classify and visualize local-global topologies

start_time <- Sys.time()
lg.E <- classLG(vr.E, 2)
end_time <- Sys.time()
message(paste("time: ", toString(end_time - start_time), "s", sep = ""))
lg.E.text <- paste("(", lg.E$degree, ", ", lg.E$cycles, ")", sep = "")
ggplot(cbind(df.E, lg.E.text), aes(x = x, y = y, colour = lg.E.text)) + 
  geom_point(size = 3) + theme_bw(base_size = 15) +
  coord_fixed() + guides(colour = guide_legend("(Degree, Cycles)")) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Construct underlying graph topology

start_time <- Sys.time()
GT.E <- constructGT(df.E, vr.E, lg.E, D.E, 6, method = "complete")
end_time <- Sys.time()
message(paste("time: ", toString(end_time - start_time), "s", sep = ""))
E <- sapply(E(GT.E), function(e) as.vector(get.edges(GT.E, e)))
E <- data.frame(matrix(unlist(apply(E, 2, function(e) 
  cbind(df.E[GT.E$centers[e[1]],], df.E[GT.E$centers[e[2]],]))), 
  ncol = 4, byrow = TRUE))
ggplot(cbind(df.E, lg.E.text, GT.E$membership), aes(x = x, y = y)) + 
  geom_point(size = 3, aes(shape = lg.E.text, colour = GT.E$membership)) +
  geom_point(data = df.E[GT.E$centers,], col = 'black', size = 3) + 
  theme_bw(base_size = 15) + coord_fixed() + 
  guides(shape = guide_legend("(Degree, Cycles)"), colour = guide_legend("Cluster")) +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  geom_segment(data = E, aes(x = X1, y = X2, xend = X3, yend = X4))
