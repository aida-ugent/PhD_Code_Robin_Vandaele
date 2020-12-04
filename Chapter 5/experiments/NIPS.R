# Load libraries
devtools::load_all() # load BCB
library("ggplot2") # plotting

# Load, prepare, and plot data

authors_meta <- rjson::fromJSON(file="Data/NIPS/author_meta.json")
G <- simplify(graph_from_data_frame(read.table("Data/NIPS/author_adjacency.csv", header=FALSE, sep=",",
                                               col.names = c("author1", "author2", "weight")),
                                    directed=FALSE))
comp <- components(G)
members <- which(comp$membership == which.max(comp$csize))
G <- induced_subgraph(G, members)
V(G)$name <- as.character(seq(length(V(G))))
V(G)$author <- unlist(rjson::fromJSON(file="Data/NIPS/id_author.json")[members])
E(G)$weight <- 2 / E(G)$weight
set.seed(420)
layout <- layout_with_lgl(G)
plot(G, vertex.size=1, vertex.label=NA, vertex.frame.color=rgb(0.5, 0.5, 0.5, 0.25),
     edge.color=rgb(0.5, 0.5, 0.5, 0.1), layout=layout)

# Conduct backbone pipeline

BCB <- backbone(G, max_leaves=10, assign=TRUE)

# View boundary coefficients

color_breaks <- c("blue", "green", "yellow")
colors <- colorRampPalette(color_breaks)(length(BCB$f))[rank(BCB$f)]
legend_colors <- colorRampPalette(rev(color_breaks))(length(BCB$f))
legend_image <- as.raster(matrix(legend_colors, ncol=1))
layout(matrix(1:2, ncol=2), width=c(2, 1), height=c(1, 1))
plot(G, vertex.size=1, vertex.label=NA, vertex.frame.color=rgb(0.5, 0.5, 0.5, 0.25),
     vertex.color=colors, edge.color=rgb(0.5, 0.5, 0.5, 0.1), layout=layout)
plot(c(0,2), c(0,1), type="n", axes=F, xlab="", ylab="", main="BC")
text(x=1.5, y=seq(0, 1, l=5), labels=round(seq(min(BCB$f), max(BCB$f), l=5), 2))
rasterImage(legend_image, 0, 0, 1, 1)
dev.off()

# View backbone

plot(BCB$B, vertex.label.color="black", vertex.size=5, vertex.label.cex=0.75, palette=categorical_pal(8),
     vertex.label=V(BCB$B)$author, vertex.label.degree=90, vertex.color=BCB$col[as.integer(V(BCB$B)$name)],
     vertex.frame.color=rgb(0.5, 0.5, 0.5, 0.25), edge.color=BCB$palette[BCB$branch], edge.width = 10,
     layout=layout_as_tree(BCB$B))

# View cost according to the number of leaves

ggplot(BCB$cost, aes(x=leaves, y=cost)) +
  geom_line(aes(group=component, col=component), size=1.5) +
  geom_vline(xintercept=length(which(degree(BCB$B) == 1)), linetype="dashed", size=1.5) +
  geom_point(aes(group=component, col=component, size=1.5)) +
  xlab("number of leaves") +
  ylab("relative cost") +
  theme_bw() +
  theme(text = element_text(size=20), legend.position="none")

# Increase the number of leaves

BCB <- get_new_leaves(BCB, 5, G, assign=TRUE)
plot(BCB$B, vertex.label.color="black", vertex.size=5, vertex.label.cex=1.35, palette=categorical_pal(8),
     vertex.label=V(BCB$B)$author, vertex.label.degree=90, vertex.color=BCB$col[as.integer(V(BCB$B)$name)],
     vertex.frame.color=rgb(0.5, 0.5, 0.5, 0.25), edge.color=BCB$palette[BCB$branch], edge.width = 10,
     layout=layout_as_tree(BCB$B), vertex.label.font=2, main="NeurIPS", label.degree = pi/4)

# View assignment induced by the backbone

plot(G, vertex.size=1, vertex.label=NA, vertex.frame.color=rgb(0.5, 0.5, 0.5, 0.25),
     vertex.color=BCB$col, edge.color = rgb(0.5, 0.5, 0.5, 0.1), layout=layout)

# Quantitative measures backbone

MSB <- medoids_steiner_backbone(G, 5) # k-medoids Steiner Tree
FPSB <- furthest_point_backbone(G, 5) # backbone through furthest point sample
MSTB <- MST_backbone(G, 5) # backbone through regular MST

backbones_evaluate(list(MSB, FPSB, MSTB, BCB$B), G)

# Check for correlation between year since active/citations and backbone
# compare with mst and full built tree trough furthest point sampling

n_citations <- sapply(authors_meta[V(G)$author], function(l) l$n_citation)
year_since_active <- sapply(authors_meta[V(G)$author], function(l) l$year)

current_pruned_mst <- mst(G)
current_pruned_fps <- furthest_point_backbone(G, length(V(G)))
current_pruned_BCB <- BCB$pine
original_size <- length(V(G))
pruning <- data.frame(times = 0, MST_size = 1, FPS_size = 1, BCB_size = 1,
                      MST_avg_citations = mean(n_citations), FPS_avg_citations = mean(n_citations), BCB_avg_citations = mean(n_citations),
                      MST_avg_year = mean(year_since_active), FPS_avg_year = mean(year_since_active), BCB_avg_year = mean(year_since_active))

for(i in seq(length(V(G)))){
  pruning[i + 1, 1] <- i
  if(length(V(current_pruned_mst)) > 3) current_pruned_mst <- hard_prune(current_pruned_mst, times=1)
  if(length(V(current_pruned_BCB)) > 3) current_pruned_fps <- hard_prune(current_pruned_fps, times=1)
  if(length(V(current_pruned_BCB)) > 3) current_pruned_BCB <- hard_prune(current_pruned_BCB, times=1)
  pruning[i + 1, 2] <- length(V(current_pruned_mst)) / original_size
  pruning[i + 1, 3] <- length(V(current_pruned_fps)) / original_size
  pruning[i + 1, 4] <- length(V(current_pruned_BCB)) / original_size
  pruning[i + 1, 5] <- mean(n_citations[V(current_pruned_mst)$author])
  pruning[i + 1, 6] <- mean(n_citations[V(current_pruned_fps)$author])
  pruning[i + 1, 7] <- mean(n_citations[V(current_pruned_BCB)$author])
  pruning[i + 1, 8] <- mean(year_since_active[V(current_pruned_mst)$author])
  pruning[i + 1, 9] <- mean(year_since_active[V(current_pruned_fps)$author])
  pruning[i + 1, 10] <- mean(year_since_active[V(current_pruned_BCB)$author])
  if(max(length(V(current_pruned_mst)), length(V(current_pruned_fps)), length(V(current_pruned_BCB))) <= 3) break
}

P1 <- ggplot(reshape2::melt(pruning[,c(1, 2, 3, 4)], id="times"), aes(x=times, y=value, col=variable)) +
  ylab("Fraction of original graph size") +
  xlab("Times pruned") +
  scale_color_manual(name="", labels = c("regular MST", "FPS tree", "BC-pine"), values = c("red", "green", "blue")) +
  geom_line() +
  theme_bw() +
  theme(legend.position="bottom", legend.text=element_text(size=25), text=element_text(size=20))

P2 <- ggplot(reshape2::melt(pruning[,c(1, 5, 6, 7)], id="times"), aes(x=times, y=value, col=variable)) +
  ylab("Mean number of citations") +
  xlab("Times pruned") +
  scale_color_manual(name="", labels = c("regular MST", "FPS tree", "BC-pine"), values = c("red", "green", "blue")) +
  geom_line() +
  theme_bw() +
  theme(legend.position="bottom", legend.text=element_text(size=25), text=element_text(size=20))

P3 <- ggplot(reshape2::melt(pruning[,c(1, 8, 9, 10)], id="times"), aes(x=times, y=value, col=variable)) +
  ylab("Mean year since active") +
  xlab("Times pruned") +
  scale_color_manual(name="", labels = c("regular MST", "FPS tree", "BC-pine"), values = c("red", "green", "blue")) +
  geom_line() +
  theme_bw() +
  theme(legend.position="bottom", legend.text=element_text(size=25), text=element_text(size=20))

ggpubr::ggarrange(P1, P2, P3, ncol=3, legend="bottom", common.legend=TRUE)
