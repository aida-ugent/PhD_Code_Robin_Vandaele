# Load libraries

library("ggplot2") # plotting

# Sample and plot data

set.seed(42)
n <- 50
theta <- runif(2 * n, 0, 2 * pi)
ycenters <- c(0, 3)
df <- data.frame(x=cos(theta) + rnorm(2 * n, sd=0.1), y=rep(ycenters, each=n) + sin(theta) + rnorm(2 * n, sd=0.1))
ggplot(df, aes(x=x, y=y)) +
  geom_point(size=1) +
  coord_fixed() +
  theme_bw()

# Persistent homology for connected components and cycles

filtration <- TDA::ripsFiltration(df, maxdimension=1, maxscale=max(dist(df)) / 2 + 0.1)
diag <- TDA::filtrationDiag(filtration, maxdimension=1, library="Dionysus", printProgress=TRUE)
op <- par(mar = c(3.25, 1, 1, 1))
TDA::plot.diagram(diag[["diagram"]], diagLim=c(0,  2.5), barcode=TRUE)
legend(2, 100, legend=c("H0", "H1"), col=c("black", "red"), lty=1, lwd=2, box.lty=0)
par(op)
op <- par(mar = c(3.25, 3.25, 1, 1))
TDA::plot.diagram(diag[["diagram"]], diagLim=c(0, 2.5))
legend(2, 1, legend=c("H0", "H1"), col=c("black", "red"), pch=c(19, 2), pt.lwd=2, box.lty=0)
par(op)

# View simplicial complex for different distance parameters

epsilons <- c(0.01, 0.5, 0.75, 1, 1.6, 2.5)
simpPlots <- list()
for(idx in 1:length(epsilons)){
  if(idx == 1){
    m <- -Inf
    nodes <- data.frame(x=numeric(0), y=numeric(0))
    edges <- data.frame(x1=numeric(0), y1=numeric(0), x2=numeric(0), y2=numeric(0))
    triangles  <- data.frame(id=integer(0), x=numeric(0), y=numeric(0))
    tId <- 0
  }
  for(simp_idx in 1:length(filtration$cmplx)){
    if(filtration$values[simp_idx] <= epsilons[idx] & filtration$values[simp_idx] > m){
      simp <- filtration$cmplx[[simp_idx]]
      if(length(simp)==1) nodes[nrow(nodes) + 1,] <- df[simp,]
      else if(length(simp)==2) edges[nrow(edges) + 1,] <- as.numeric(c(df[simp[1],], df[simp[2],]))
      else if(length(simp==3)){
        tId <- tId + 1
        triangles[(nrow(triangles) + 1):(nrow(triangles) + 3),] <- cbind(tId, df[simp,])
      }
    }
  }
  m <- epsilons[idx]
  simpPlots[[length(simpPlots) + 1]] <- ggplot(nodes, aes(x=x, y=y)) +
    geom_segment(data=edges, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", size=0.5, alpha=0.5) +
    geom_point(size=1.5, alpha=0.75) +
    geom_polygon(data=triangles, aes(group=id), fill="green", alpha=0.25) +
    geom_segment(data=edges, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", size=0.5, alpha=0.05) +
    geom_point(size=1.5, alpha=0.25) +
    coord_fixed() +
    theme_bw() +
    ggtitle(latex2exp::TeX(sprintf("$\\epsilon = %g$", epsilons[idx]))) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    theme(plot.title = element_text(hjust = 0.5, size=22), text = element_text(size=20))
}

gridExtra::grid.arrange(grobs=simpPlots, ncol=length(epsilons))
