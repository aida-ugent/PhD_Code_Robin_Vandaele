#' Construct a graph topology based on a local-global classification
#'
#' @param df     data frame of point cloud data with an underlying 
#'               graph-structured topology
#' @param g      proximity graph of the data set
#' @param lg     data frame containg local-global classification of 
#'               considered data set
#' @param D      pairwise distance object of full data frame
#' @param r      distance parameter to force separation of emerging 
#'               branches
#' @param method linking method to be used to force separation of 
#'               emerging branches by means of hierarchical clustering,
#'               standard set to "complete" 
#'           
#' @return      an igraph object representing the underlying topology 
#'              of df

constructGT <- function(df, g, lg, D, r, method = "complete"){
  M <- list(); N <- integer(0); centers <- integer(0); adj <- list()
  I <- which(lg$degree > 2)
  groups <- 2^lg[I, 1] * 3^lg[I, 2]
  names <- unique(groups)
  clusts <- list("membership" = integer(length(I)), "no" = 0)
  for(name in names){
    J <- which(groups == name)
    Cs <- clusters(induced_subgraph(g, I[J]))
    clusts$membership[J] <- clusts$no + Cs$membership
    clusts$no <- clusts$no + Cs$no
  }
  if(clusts$no > 0){
    N <- 1:clusts$no
    M <- lapply(N, function(i) I[which(clusts$membership == i)])
    centers <- sapply(N, function(i) center(M[[i]], D))
    cum_nclusts <- c(0, cumsum(lg[centers, "degree"]))
    adj <- rep(list(integer(0)), tail(N, 1) + tail(cum_nclusts, 1))
    adj[N] <- lapply(N, function(i) {
      (length(N) + cum_nclusts[i] + 1):(length(N) + cum_nclusts[i + 1])
    })
    adj_bifs <- integer(length(N))
    for(i in N[-length(N)]){
      for(j in (i + 1):length(N)){
        if(isConnected(g, M[[i]], M[[j]])){
          adj_bifs[i] <- adj_bifs[i] + 1
          adj_bifs[j] <- adj_bifs[j] + 1
          Mi <- adj[[i]][adj_bifs[i]]
          Mj <- adj[[j]][adj_bifs[j]]
          adj[[i]][adj_bifs[i]] <- j
          adj[[j]][adj_bifs[j]] <- i
          adj <- lapply(adj[1:(length(adj) - 2)], function(x){
            newx <- x
            newx[x > Mi] <- x[x > Mi] - 1
            newx[x > Mj] <- x[x > Mj] - 2
            return(newx)
          })
        }
      }
    }
    I <- which((lg$degree == 1) | (lg$degree == 2))
    for(i in N){
      M[(length(M) + 1):(length(M) + lg[centers[i], "degree"] - adj_bifs[i])] <-
        separate(centers[i], I, r, g, lg[centers[i], "degree"] - adj_bifs[i], D, method)
    }
    if(length(N) > 1){
      OV <- matrix(apply(expand.grid(1:length(M), 1:length(M)), 1, 
                         function(x){
                          sum(M[[x[1]]] %in% M[[x[2]]]) / length(M[[x[1]]])
                           }), nrow = length(M))
      OV <- which(OV > 0, arr.ind = TRUE)
      OV <- OV[OV[,1] < OV[,2],]
      while(nrow(OV) > 0){
        for(i in seq(nrow(OV))){
          M[c(OV[i, 1], OV[i, 2])] <- 
            splitsOverlapping(centers, M, OV[i, 1], OV[i, 2], adj, lg, D)
        }
        OV <- matrix(apply(expand.grid(1:length(M), 1:length(M)), 1, 
                           function(x){
                             sum(M[[x[1]]] %in% M[[x[2]]]) / length(M[[x[1]]])
                           }), nrow = length(M))
        OV <- matrix(which(OV > 0, arr.ind = TRUE), ncol = 2)
        OV <- matrix(OV[OV[,1] < OV[,2],],  ncol = 2)
      }
    }
  }
  I <- setdiff(1:nrow(df), unlist(M))
  if(length(I) > 0){
    groups <- 2^lg[I, 1] * 3^lg[I, 2]
    names <- unique(groups)
    groups <- lapply(names, function(i) I[which(groups == i)])
    clusts <- lapply(groups, function(x) clusters(induced_subgraph(g, x)))
    i <- which(names == 12)
    if(length(i) == 1){
      NO <- clusts[[i]]$no
      for(c in 1:NO){
        W <- which(clusts[[i]]$membership == c)
        C <- groups[[i]][W]
        if(!isConnected(g, unlist(M), C)){
          newM <- separateCycle(g, D, C, method)
          J <- newM == 4
          newM[J] <- c
          newM[!J] <- newM[!J] + clusts[[i]]$no
          clusts[[i]]$membership[W] <- newM
          clusts[[i]]$no <- clusts[[i]]$no + 3
        }
      }
    }
    cum_nclusts <- c(0, cumsum(sapply(clusts, function(l) l$no)))
    n_clusts <- tail(cum_nclusts, 1)
    M[(length(M) + 1):(length(M) + n_clusts)] <- rep(list(NA), n_clusts)
    for(i in seq(length(groups))){
      M[(length(M) - n_clusts + cum_nclusts[i] + 1):
          (length(M) - n_clusts + cum_nclusts[i] + clusts[[i]]$no)] <- 
        lapply(1:clusts[[i]]$no, 
               function(j) groups[[i]][which(clusts[[i]]$membership == j)])
      }
    centers[(length(N) + 1):length(M)] <- 
      sapply((length(N) + 1):length(M), function(i) center(M[[i]], D))
    adj[(length(M) - n_clusts + 1):length(M)] <- rep(list(integer(0)), n_clusts)
    for(i in (length(M) - n_clusts + 1):length(M)){
      for(j in (length(N) + 1):length(M)){
        if(i != j && !(j %in% adj[[i]])){
          if(isConnected(g, M[[i]], M[[j]])){
            adj[[i]] <- c(adj[[i]], j)
          }
        }
      }
    }
    for(u in seq(length(adj))){
      for(v in adj[[u]]){
        adj[[v]] <- c(adj[[v]], u)
      } 
    }
    for(i in N){
      for(u in adj[[i]]){
        if(length(adj[[u]]) < lg[centers[u], "degree"]){
          for(v in unlist(adj[N[-(1:i)]])){
            if(length(adj[[v]]) < lg[centers[v], "degree"]){
              if(isConnected(g, M[[u]], M[[v]])){
                adj[[u]] <- c(adj[[u]], v)
                adj[[v]] <- c(adj[[v]], u)
              }
            }
            if(length(adj[[u]]) == lg[centers[u], "degree"]) break
          }
        }
      }
    }
  }
  GT <- as.undirected(graph_from_adj_list(adj))
  I <- which(degree(GT) == 1)
  for(i in I){
    centers[i] <- M[[i]][which.max(as.matrix(D)[
      as.matrix(expand.grid(
        M[[i]], centers[as.integer(neighbors(GT, i))]))])]
  }
  I <- order(unlist(M))
  M <- sapply(M, function(x) length(x))
  M <- factor(unlist(sapply(1:length(M), function(i) rep(i, M[i]))))[I]
  GT$membership <- M
  GT$centers <- centers
  return(GT)
}

#' Compute a distance object D_I as a distance object D,
#' restricted to the set I X I
#'
#' @param D pairwise distance object of full data frame
#' @param I indices for which all pairwise distances are to be 
#'          extracted
#' 
#' @return  distance object D restricted to I

subD <- function(D, I){
  D_I <- expand.grid(I, I)
  D_I <- D_I[which(D_I[,2] < D_I[,1]),]
  D_I <- attr(D, "Size") * (D_I[,2] - 1) - 
    D_I[,2] * (D_I[,2] - 1) / 2 + D_I[,1] - D_I[,2]
  D_I <- D[D_I]
  attr(D_I, "Size") <- length(I)
  return(D_I)
}

#' Separate branches near a bifurcation point
#'
#' @param c      center for which emerging branches are to be seperated
#' @param I      indices of unassigned data points
#' @param r      distance parameter to force separation of emerging 
#'               branches
#' @param g      proximity graph of the data set
#' @param no     number of clusters to separate
#' @param D      pairwise distance object of full data frame
#' @param method linking method to be used to force separation of 
#'               emerging branches by means of hierarchical clustering
#'           
#' @return       list of membership to detected branches

separate <- function(c, I, r, g, no, D, method){
  adj_to_c <- intersect(I, unlist(ego(g, r, c)))
  clusters <- hclust(subD(D, adj_to_c), method = method)
  cut <- as.factor(cutree(clusters, no))
  L <- lapply(1:no, function(i) adj_to_c[which(cut == i)])
  return(L)
}

#' Split two overlapping membership sets according to distance to
#' adjacent bifurcations nodes in two approximately equal sized sets
#'
#' @param centers indices of bifurcation centers
#' @param M       list of indices representing membership to 
#'                different clusters
#' @param i       first index in M of overlapping indices
#' @param j       second index in M of overlapping indices
#' @param adj     current constructed adjacency list
#' @param D       pairwise distance object of full data frame
#' @param lg     data frame containg local-global classification of 
#'               considered data set
#'           
#' @return       list of two splitted overlapping components

splitsOverlapping <- function(centers, M, i, j, adj, lg, D){
  U <- union(M[[i]], M[[j]])
  bi <- which.max(i - length(centers) <= cumsum(lg[centers, "degree"]))
  Ord <- order(as.matrix(D)[as.matrix(expand.grid(centers[bi], U))])
  return(list(U[Ord[1:floor(length(U) / 2)]],
                U[Ord[-(1:floor(length(U) / 2))]]))
}

#' Separate an isolated cycle in four connected parts
#'
#' @param g       proximity graph of the data set
#' @param D       pairwise distance object of full data frame
#' @param I       indices of the circle in the full data frame
#' @param method  linking method to be used to force separation
#'                of two non adjacent components in the circle 
#'           
#' @return        new vector of membership to one of four connected parts
#'                of the circle

separateCycle <- function(g, D, I, method){
  BFS <- bfs(induced_subgraph(g, I), root = 1, dist = TRUE, 
             unreachable = FALSE)
  m <- max(BFS$dist)
  I2 <- which(BFS$dist %in% (floor(m / 3) + 1):(2 * floor(m / 3)))
  I3 <- which(BFS$dist %in% (2 * floor(m / 3) + 1):m)
  cut2 <- cutree(hclust(subD(D, I[I2]), method = method), 2)
  I4 <- I2[which(cut2 == 1)]
  I2 <- I2[which(cut2 == 2)]
  memb <- rep(1, length(I))
  memb[I2] <- 2
  memb[I3] <- 3
  memb[I4] <- 4
  return(memb)
}

#' Determine if two sets of vertices are connected in a given graph
#'
#' @param g  proximity graph of the data set
#' @param I1 first set of vertex indices
#' @param I2 second set of vertex indices
#' 
#' @return   boolean value indicating if there is an edge between I1 and I2

isConnected <- function(g, I1, I2){
  neighboring <- unlist(ego(g, 1, nodes = I1))
  isCon <- FALSE
  for(v in I2){if(v %in% neighboring){isCon <- TRUE; break}}
  return(isCon)
}
