# Load libraries
devtools::load_all() # load BCB

# Load and plot data

load("Data/GoT/union_edges.RData")
load("Data/GoT/union_characters.RData")
G <- graph_from_data_frame(union_edges, vertices=union_characters, directed=FALSE)
color_vertices <- union_characters %>% group_by(house, color) %>% summarise(n = n()) %>% filter(!is.na(color))
colors <- as.factor(color_vertices$color)
newColors <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9')
color_vertices$color <- newColors[sapply(color_vertices$color, function(col) which(levels(colors) == col))]
G <- simplify(G)
V(G)$name[V(G)$name=="Walda Frey daughter of Merrett"] <- "Walda Frey"
V(G)$color[!is.na(V(G)$color)] <-  sapply(V(G)$color[!is.na(V(G)$color)], function(col) newColors[which(levels(colors) == col)])
V(G)$color[is.na(V(G)$color)] <- "white"
set.seed(69)
layout <- layout_with_fr(G)
plot(G, layout=layout, vertex.label=gsub(" ", "\n", V(G)$name), vertex.shape=V(G)$shape,
     vertex.color=V(G)$color, vertex.size=(V(G)$popularity + 0.5) * 5, vertex.frame.color="gray",
     vertex.label.color="black", vertex.label.cex=0.8, edge.arrow.size=0.5, edge.color=E(G)$color,
     edge.lty=E(G)$lty)
legend(x=-2, y=1.5, legend=c(NA, NA, as.character(color_vertices$house)), pch=19, pt.cex=5, cex=1.75,
       col=c(NA, NA, color_vertices$color), bty="n", ncol=1, y.intersp=0.3, x.intersp=0.3)

# Conduct backbone pipeline with ordinary LCC

old <- Sys.time()
f <- transitivity(G, type="local")
new <- Sys.time() - old
print(paste("Local cluster coefficients obtained in", round(new, 3), attr(new, "units")))
f[degree(G) == 1] <- 1
LCCB <- backbone(G, f=f)

# View LCC and LCC-pine

color_breaks <- c("blue", "green", "yellow")
colors <- colorRampPalette(color_breaks)(length(f))[rank(f)]
legend_colors <- colorRampPalette(rev(color_breaks))(length(f))
legend_image <- as.raster(matrix(legend_colors, ncol = 1))
E(G)$color <- rep("gray", length(E(G)))
E(G)$color[get.edge.ids(G, t(ends(LCCB$pine, E(LCCB$pine))))] <- "red"
layout(matrix(1:2, ncol=2), width=c(2, 1), height=c(1, 1))
plot(G, layout=layout, vertex.label=gsub(" ", "\n", V(G)$name), vertex.shape=V(G)$shape,
     vertex.color=colors, vertex.size=(V(G)$popularity + 0.5) * 5, vertex.frame.color="gray",
     vertex.label.color="black", vertex.label.cex=0.8, edge.arrow.size=0.5, edge.lty=E(G)$lty)
plot(c(0,2), c(0,1), type="n", axes=FALSE, xlab="", ylab="", main="LCC")
text(x=1.5, y = seq(min(f), max(f), l=5), labels=seq(min(f), max(f), l=5))
rasterImage(legend_image, 0, 0, 1, 1)
dev.off()

# Increase number of leaves and view backbone

LCCB <- get_new_leaves(LCCB, 10, G, stdize=TRUE)
E(G)$color <- adjustcolor(rep("gray", length(E(G))), 0.5)
E(G)$color[get.edge.ids(G, t(ends(LCCB$B, E(LCCB$B))))] <- adjustcolor("red", 0.75)
edge.width <- rep(0.75, length(E(G)))
edge.width[get.edge.ids(G, t(ends(LCCB$pine, E(LCCB$pine))))] <- 3
backboneI <- which(V(G)$name %in% V(LCCB$B)$name)
V(G)$color[-backboneI] <- adjustcolor(V(G)$color[-backboneI], 0.75)
vertex.size <- rep(5, length(V(G)))
vertex.size[-backboneI] <- 1
vertex.label <- V(G)$name
vertex.label[-backboneI] <- ""
vertex.label <- gsub(" ", "\n", vertex.label)
vertex.label[vertex.label=="Asha\n(Yara)\nGreyjoy"] <- "Asha (Yara)\nGreyjoy"
vertex.label[vertex.label=="Damon\nLannister\nson\nof\nJason"] <- "Damon Lannister\nson of Jason"
vertex.label[vertex.label=="Brandon\nStark\nson\nof\nCregan"] <- "Brandon Stark\nson of Cregan"
vertex.label[vertex.label=="Princess\nof\nDorne"] <- "Princess of\nDorne"
vertex.label[vertex.label=="Aegon\nV\nTargaryen"] <- "Aegon V Targaryen"
vertex.label[vertex.label=="Unknown\nmother\nTyrell" ] <- "Unknown\nmother\nTyrell"
vertex.label[vertex.label=="Aerys\nII\nTargaryen"] <- "Aerys II\nTargaryen"
vertex.label[vertex.label=="Luthor\nTyrell\nson\nof\nMoryn"] <- "Luthor\nTyrell\n son of Moryn"
vertex.label[vertex.label=="Elia\nMartell"] <- "Elia Martell"
vertex.label[vertex.label=="Tywin\nLannister"] <- "Tywin Lannister"
vertex.label[vertex.label=="Tommen\nBaratheon"] <- "Tommen Baratheon"
vertex.label[vertex.label=="Betha\nBlackwood"] <- "Betha Blackwood"
vertex.label[vertex.label=="Duncan\nTargaryen"] <- "Duncan Targaryen"
vertex.label[vertex.label=="Shaera\nTargaryen"] <- "Shaera Targaryen"
vertex.label[vertex.label=="Rhaelle\nTargaryen"] <- "Rhaelle Targaryen"
vertex.label[vertex.label=="Jon\nSnow"] <- "Jon Snow"
vertex.label[vertex.label=="Edwyle\nStark"] <- "Edwyle Stark"
vertex.label[vertex.label=="Rickard\nStark"] <- "Rickard Stark"
vertex.label[vertex.label=="Eddard\nStark"] <- "Eddard Stark"
vertex.label[vertex.label=="Margaery\nTyrell"] <- "Margaery Tyrell"
layout[V(G)$name=="Aegon V Targaryen", 2] <- layout[V(G)$name=="Aegon V Targaryen", 2] + 0.5
layout[V(G)$name=="Jon Snow", 2] <- layout[V(G)$name=="Jon Snow", 2] + 0.5
layout[V(G)$name=="Tommen Baratheon", 2] <- layout[V(G)$name=="Tommen Baratheon", 2] - 0.5
layout[V(G)$name=="Luthor Tyrell", 1] <- layout[V(G)$name=="Luthor Tyrell", 1] + 0.25
layout[V(G)$name=="Damion Lannister", 2] <- layout[V(G)$name=="Damion Lannister", 2] - 0.5
par(mar=rep(0, 4))
plot(G, layout=layout, vertex.label=vertex.label, vertex.size=vertex.size, vertex.frame.color="gray",
     vertex.shape="circle", vertex.label.color="black", vertex.label.font=2, vertex.label.cex=1, edge.width=edge.width, asp=1.25)
legend(x=-2, y=1, legend=c(NA, NA, as.character(color_vertices$house)), pch=19, pt.cex=3, cex=1.25,
       col=c(NA, NA, color_vertices$color), bty="n", ncol=1, y.intersp=1.5, x.intersp=0.75)
dev.off()

# Quantitative measures backbone

MSB <- medoids_steiner_backbone(G, c(8, 2)) # k-medoids Steiner Tree
FPSB <- furthest_point_backbone(G, c(8, 2)) # backbone through furthest point sample
MSTB <- MST_backbone(G, c(8, 2)) # backbone through regular MST

backbones_evaluate(list(MSB, FPSB, MSTB, LCCB$B), G)

# Check for cycles in the backbone

LCCB <- check_for_cycles(LCCB)
TDA::plot.diagram(LCCB$diagram, diagLim=c(0, 7))
legend(6, 3, legend=c("H0", "H1"), col=c("black", "red"), pch=c(19, 2), pt.lwd=2, box.lty=0)

maxPersistingCycles <- order((LCCB$diagram[,"Death"] - LCCB$diagram[,"Birth"])[LCCB$diagram[,"dimension"]==1], decreasing=TRUE)[1:2]
Ecocycles <- unique(do.call("rbind", LCCB$repCycles[maxPersistingCycles]))
GwithRepCycle <- add_edges(G, as.character(t(Ecocycles)))
E(GwithRepCycle)$color[is.na(E(GwithRepCycle)$color)] <- adjustcolor("orange")
new.edge.width <- c(edge.width, rep(3, nrow(Ecocycles)))

par(mar=rep(0, 4))
plot(GwithRepCycle, layout=layout, vertex.label=vertex.label, vertex.size=vertex.size, vertex.frame.color="gray",
     vertex.shape="circle", vertex.label.color="black", vertex.label.font=2, vertex.label.cex=1, edge.width=new.edge.width, asp=1.25)
legend(x=-2, y=1, legend=c(NA, NA, as.character(color_vertices$house)), pch=19, pt.cex=3, cex=1.25,
       col=c(NA, NA, color_vertices$color), bty="n", ncol=1, y.intersp=.45, x.intersp=0.75)
dev.off()
