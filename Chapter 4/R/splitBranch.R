#' Divide a branch into a given number of consecutive clusters according to a
#' breadth first search traversion
#'
#' @param GT constructed underlying graph topology of the data represented by g
#' @param g  proximity graph of the data set
#' @param D  pairwise distance object of the full data frame
#' @param k  number of consecutive clusters 
#' @param i  index of the cluster in GT to be divided
#' 
#' @return   modified underlying graph topology with given branch divided into k components

splitBranch <- function(GT, g, D, i, k){
 I <- which(GT$membership == i)
 subD <- distances(g, I, I)
 root <- which(subD == max(subD), arr.ind = TRUE)[1]
 bfs <- bfs(induced_subgraph(g, I), root = root, dist = TRUE)
 m <- max(bfs$dist)
 newM <- lapply(1:(k - 1), function(t){
   which(bfs$dist %in% (floor((t - 1) * m / k) + 1):(floor(t * m / k)))
 })
 newM[[1]] <- c(root, newM[[1]])
 newM[[k]] <- setdiff(1:length(I), unlist(newM))
 s <- max(V(GT))
 levels(GT$membership) <- c(levels(GT$membership), (s + 1):(s + k - 1))
 for(t in 1:(k - 1)){
   GT$membership[I][newM[[t]]] <- s + t
 }
 N <- neighbors(GT, i)
 GT <- delete_edges(GT, incident(GT, i))
 oldM <- lapply(N, function(t) which(GT$membership == t))
 GT <- add_vertices(GT, k - 1)
 GT <- add_edges(GT, c(sapply((s + 1):(s + k - 2), function(t) c(t, t + 1))))
 GT <- add_edges(GT, c(i, s + k - 1))
 GT$centers[c((s + 1):(s + k - 1), i)] <- 
   sapply(1:k, function(j) center(I[newM[[j]]], D))
 for(n in N){
   subD <- mean(as.matrix(D)[as.matrix(expand.grid(newM[[1]], oldM[[1]]))])
 }
 newN <- sapply(seq(N), function(j){
   c((s+ 1):(s + k - 1), i)[which.min(sapply(seq(newM), 
                                             function(t) mean(as.matrix(D)[
                                               as.matrix(expand.grid(I[newM[[t]]], 
                                                                     oldM[[j]]))])))]
 })
 GT <- add_edges(GT, c(rbind(N, newN)))
 for(j in seq(newN)){
   n <- intersect(neighbors(GT, newN[j]), c((s + 1):(s + k - 1), i))
   GT$centers[newN[j]] <- oldM[[j]][which.max(as.matrix(D)[
     as.matrix(expand.grid(oldM[[j]], GT$centers[n]))])]
 }
 return(GT)
}