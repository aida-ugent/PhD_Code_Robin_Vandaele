#' Evaluate each backbone in a list according to various metrics
#'
#' @param BL list of backbone graphs
#' @param G  the orignal graph
#'
#' @return  a data frame with length(L) rows, where each row represents the evaluation
#'          of the corresponding backbone according to various metrics

backbones_evaluate <- function(BL, G){
  Comp <- components(G)
  dG <- distances(G)
  uv <- which(dG == max(dG[dG != Inf]), arr.ind=TRUE)[1,]
  memberNodes <- lapply(1:Comp$no, function(C) V(G)$name[which(Comp$membership == C)])
  rm(Comp)
  Center <- sapply(memberNodes, function(N) N[which.min(apply(dG[N, N], 1, max, na.rm=TRUE))])
  totalVar <- sum(apply(as.matrix(dG[,Center], ncol=length(Center)), 1, min, na.rm=TRUE))
  print("Computing original average commute times...")
  OmegaG <- matrix(rep(NaN, length(V(G))^2), nrow=length(V(G)))
  rownames(OmegaG) <- names(V(G))
  colnames(OmegaG) <- names(V(G))
  for(N in memberNodes) OmegaG[N, N] <- linkprediction::proxfun(induced_subgraph(G, N), N, N, method="act")
  diag(OmegaG) <- 0
  metrics <- do.call("rbind", lapply(1:length(BL), function(i){
    print(paste0("Computing projection on backbone ", i, "..."))
    GprojB <- project_on_backbone(G, BL[[i]], dG)
    print(paste0("Computing average commute times of projection on backbone ", i, "..."))
    OmegaGprojB <- matrix(rep(NaN, length(V(GprojB))^2), nrow=length(V(G)))
    rownames(OmegaGprojB) <- names(V(G))
    colnames(OmegaGprojB) <- names(V(G))
    for(N in memberNodes) OmegaGprojB[N, N] <- linkprediction::proxfun(induced_subgraph(GprojB, N), N, N, method="act")
    diag(OmegaGprojB) <- 0
    return(data.frame(V = length(V(BL[[i]])) / length(V(G)),
                      E = length(E(BL[[i]])) / length(E(G)),
                      R = 1 - sum(dG[cbind(GprojB$addedV, GprojB$connectTo)]) / totalVar,
                      sigma =  dG[uv[1], uv[2]] / distances(GprojB, names(V(G))[uv[1]], names(V(G))[uv[2]])[1, 1],
                      cor_act = cor(as.numeric(OmegaG), as.numeric(OmegaGprojB), use="complete.obs"),
                      leaves = sum(degree(BL[[i]]) == 1)))
  }))
  return(metrics)
}
