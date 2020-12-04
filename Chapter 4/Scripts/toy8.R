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

# Extract and plot 8

df.8 <- df.A[which(df.A[,1] > 87.5 & df.A[,2] < 60),]
ggplot(df.8, aes(x = V1, y = V2)) + geom_point() + coord_fixed()

# Construct and visualize Vietoris-Rips graph from data

D.8 <- dist(df.8)
eps <- 3.5
vr.8 <- VRgraph(D.8, eps)
E <- sapply(E(vr.8), function(e) as.vector(get.edges(vr.8, e)))
E <- data.frame(matrix(unlist(apply(E, 2, function(e) 
  cbind(df.8[e[1],], df.8[e[2],]))), ncol = 4, byrow = TRUE))
ggplot(df.8, aes(x = V1, y = V2)) + coord_fixed() + theme_bw(base_size = 15) +
  geom_segment(data = E, size = 1, aes(x = X1, y = X2, xend = X3, yend = X4), 
               color = 'red') + geom_point(size = 3, alpha = 0.5)

# Visualize the idea behind LTDA

circleFun <- function(center = c(0,0), r = 1, npoints = 100){
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

c <- 64
I1 <- which(as.matrix(D.8)[c,] <= 5)
I2 <- which(as.matrix(D.8)[c,] > 10)
I3 <- setdiff(seq(nrow(df.8)), c(I1, I2))
ggplot() + coord_fixed() + theme_bw(base_size = 15) +
  geom_segment(data = E, size = 1, aes(x = X1, y = X2, xend = X3, yend = X4), color = 'black') + 
  geom_point(data = df.8[I1,], aes(x = V1, y = V2), size = 5, color = 'red', alpha = 0.75) +
  geom_point(data = df.8[I2,], aes(x = V1, y = V2), size = 5, color = 'blue', alpha = 0.75) +
  geom_point(data = df.8[I3,], aes(x = V1, y = V2), size = 5, color = 'green', alpha = 0.75) +
  geom_point(data = df.8[c,], aes(x = V1, y = V2), col = 'black', size = 7) + 
  geom_path(data = circleFun(as.numeric(df.8[c,]), 5, npoints = 100), aes(x = x, y = y),
            size = 3, alpha = 0.75) +
  geom_path(data = circleFun(as.numeric(df.8[c,]), 10, npoints = 100), aes(x = x, y = y),
            size = 3, alpha = 0.75)
