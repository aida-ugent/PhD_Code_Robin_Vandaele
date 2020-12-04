# Load libraries

library("ggplot2") # plotting

# Sample and plot data

set.seed(42)
n <- 50
theta <- runif(n, 0, 2 * pi)
df <- data.frame(x=cos(theta) + rnorm(n, sd=0.1), y=sin(theta) + rnorm(n, sd=0.1))
ggplot(df, aes(x=x, y=y)) +
  geom_point(size=1) +
  coord_fixed() +
  theme_bw()

# Persistent homology for connected components and cycles

filtration <- TDA::ripsFiltration(df, maxdimension=1, maxscale=2)
diag <- TDA::filtrationDiag(filtration, maxdimension=1, library="Dionysus", printProgress=TRUE)
TDA::plot.diagram(diag[["diagram"]], diagLim=c(0,  2.5))
legend(2, 1, legend=c("H0", "H1"), col=c("black", "red"), pch=c(19, 2), pt.lwd=2, box.lty=0)

# View simplicial complex for different distance parameters

epsilons <- c(0.01, 0.5, 0.75, 2)
simpPlots <- list()
for(eps in epsilons){
  simpI <- which(filtration$values < eps)
  nodes <- data.frame(x=numeric(0), y=numeric(0))
  edges <- data.frame(x1=numeric(0), y1=numeric(0), x2=numeric(0), y2=numeric(0))
  triangles  <- data.frame(id=integer(0), x=numeric(0), y=numeric(0))
  tId <- 0
  for(simp in filtration$cmplx[simpI]){
    if(length(simp)==1) nodes[nrow(nodes) + 1,] <- df[simp,]
    else if(length(simp)==2) edges[nrow(edges) + 1,] <- as.numeric(c(df[simp[1],], df[simp[2],]))
    else if(length(simp==3)){
      tId <- tId + 1
      triangles[(nrow(triangles) + 1):(nrow(triangles) + 3),] <- cbind(tId, df[simp,])
    }
  }
  simpPlots[[length(simpPlots) + 1]] <- ggplot(nodes, aes(x=x, y=y)) +
    geom_segment(data=edges, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", size=0.5, alpha=0.5) +
    geom_point(size=1.5, alpha=0.75) +
    geom_polygon(data=triangles, aes(group=id), fill="green", alpha=0.25) +
    geom_segment(data=edges, aes(x=x1, y=y1, xend=x2, yend=y2), color="black", size=0.5, alpha=0.05) +
    geom_point(size=1.5, alpha=0.25) +
    coord_fixed() +
    theme_bw() +
    ggtitle(latex2exp::TeX(sprintf("$\\epsilon = %g$", eps))) +
    theme(plot.title = element_text(hjust = 0.5, size=22))
}

gridExtra::grid.arrange(grobs=simpPlots, ncol=2)
