# Load libraries
devtools::load_all() # load BCB

# Load and plot data

G <- simplify(graph_from_data_frame(read.table("Data/HarryPotter/relations.csv", header=TRUE, sep=",",
                                               col.names = c("potter1", "potter2", "type")),
                                    directed=FALSE))
G <- read.table("Data/HarryPotter/relations.csv", header=TRUE, sep=",", col.names = c("from", "to", "weight"))
G <- G[G[,"weight"] != "-", 1:2]
G <- simplify(graph_from_data_frame(G, directed=FALSE))
V(G)$name <- as.character(read.table("Data/HarryPotter/characters.csv", quote="", header=TRUE, sep=",")[as.integer(V(G)$name) + 1, 2])
V(G)$name[V(G)$name == "\"Bartemius \"\"Barty\"\" Crouch Jr.\""] <- "Bartemius Crouch Jr."
V(G)$name[V(G)$name == "\"Bartemius \"\"Barty\"\" Crouch Sr.\""] <- "Bartemius Crouch Sr."
set.seed(17)
layout <- layout_with_graphopt(G)
small_component1 <- which(V(G)$name %in% c("Petunia Dursley", "Vernon Dursley", "Dudley Dursley"))
small_component2 <- which(V(G)$name %in% c("Tom Riddle Sr.", "Mary Riddle"))
layout[V(G)$name=="Albus Dumbledore", 1] <- layout[V(G)$name=="Albus Dumbledore", 1] - 5
layout[small_component1, 1] <- layout[small_component1, 1] - 250
layout[small_component1, 2] <- layout[small_component1, 2] - 75
layout[small_component2, 1] <- layout[small_component2, 1] - 300
layout[small_component2, 2] <- layout[small_component2, 2] - 25
plot(G, layout=layout, vertex.label=gsub(" ", "\n", V(G)$name), vertex.shape=V(G)$shape,
     vertex.frame.color="gray", vertex.label.color="black", vertex.label.cex=0.8, edge.arrow.size=0.5,
     edge.color=E(G)$color, edge.lty=E(G)$lty)


# Conduct backbone pipeline with ordinary LCC

old <- Sys.time()
f <- transitivity(G, type="local")
new <- Sys.time() - old
print(paste("Local cluster coefficients obtained in", round(new, 3), attr(new, "units")))
f[degree(G) == 1] <- 1
LCCB <- backbone(G, f=f)

# View LCC and LCC-pine

color_breaks <- c("blue", "green", "yellow")
V(G)$color <- adjustcolor(colorRampPalette(color_breaks)(length(f))[rank(f, ties.method="max")], 0.5)
frame.color <- rep("black", length(V(G)))
frame.color[V(G)$name %in% V(LCCB$B)$name] <- "red"
legend_colors <- adjustcolor(colorRampPalette(color_breaks)(length(f)), 0.5)
legend_image <- rev(as.raster(matrix(legend_colors, ncol = 1)))
E(G)$color <- rep("gray", length(E(G)))
E(G)$color[get.edge.ids(G, t(ends(LCCB$pine, E(LCCB$pine))))] <- adjustcolor("blue", 0.35)
E(G)$color[get.edge.ids(G, t(ends(LCCB$B, E(LCCB$B))))] <- adjustcolor("red", 0.5)
edge.width <- rep(0.75, length(E(G)))
edge.width[get.edge.ids(G, t(ends(LCCB$pine, E(LCCB$pine))))] <- 1
edge.width[get.edge.ids(G, t(ends(LCCB$B, E(LCCB$B))))] <- 4
labels <- c("Harry Potter", "Ron Weasley", "Hermione Granger", "Severu Snape", "Albus Dumbledore", "James Potter", "Lily Potter",
            "Neville Longbottom", "Bartemius Crouch Jr.", "Lord Voldemort", "Lucius Malfoy", "Draco Malfoy",
            "Bartemius Crouch Sr.", "Igor Karkaroff", "Olympe Maxime", "Fluffy", "Aragog", "Mary Riddle", "Tom Riddle Sr.",
            "Rubeus Hagrid", "Dobby", "Cedric Diggory", "Sirius Black", "George Weasley", "Luna Lovegood", "Vernon Dursley", "Dudley Dursley",
            "Petunia Dursley", "Hedwig", "Moaning Myrtle", "Severus Snape", "Vincent Crabbe Sr.")
labelsI <- which(V(G)$name %in% labels)
charactersInPlot <- V(G)$name
charactersInPlot[-labelsI] <- ""
vertex.size <- rep(7, length(V(G)))
vertex.size[-labelsI] <- 1
layout(matrix(1:2, ncol=2), width=c(3, 1), height=c(1, 1))
plot(G, layout=layout, dge.width=edge.width, vertex.size=vertex.size, vertex.label.cex=1.25, vertex.label.color="black",
     vertex.frame.color=frame.color, vertex.label.font=2, asp=0, vertex.label=charactersInPlot)
plot(c(0,2), c(0,1), type="n", axes=FALSE, xlab="", ylab="", main="LCC")
text(x=1.5, y = seq(0, 1, l=5), labels=round(seq(min(f), max(f), l=5), 2), cex=1.5)
rasterImage(legend_image, 0, 0, 1, 1)
dev.off()

# Quantitative measures backbone

MSB <- medoids_steiner_backbone(G, c(2, 1, 2))
FPSB <- furthest_point_backbone(G, c(2, 1, 2))
MSTB <- MST_backbone(G, c(2, 1, 2))

backbones_evaluate(list(MSB, FPSB, MSTB, LCCB$B), G)

# Check for cycles in the backbone

LCCB <- check_for_cycles(LCCB, G=G)
op <- par(mar = c(3.25, 3.25, 1, 1))
TDA::plot.diagram(LCCB$diagram, diagLim=c(0, 3))
legend(2, 1.5, legend=c("H0", "H1"), col=c("black", "red"), pch=c(19, 2), pt.lwd=2, box.lty=0)
par(op)

Ecocycles <- do.call("rbind", LCCB$repCycles)
GwithRepCycle <- add_edges(G, as.character(t(Ecocycles)))
E(GwithRepCycle)$color[is.na(E(GwithRepCycle)$color)] <- adjustcolor("orange")
new.edge.width <- c(edge.width, rep(3, nrow(Ecocycles)))

layout(matrix(1:2, ncol=2), width=c(3, 1), height=c(1, 1))
plot(GwithRepCycle, layout=layout, edge.width=new.edge.width, vertex.size=vertex.size, vertex.label.cex=1.25, vertex.label.color="black",
     vertex.frame.color=frame.color, vertex.label.font=2, asp=0, vertex.label=charactersInPlot)
plot(c(0,2), c(0,1), type="n", axes=FALSE, xlab="", ylab="", main="LCC")
text(x=1.5, y = seq(0, 1, l=5), labels=round(seq(min(f), max(f), l=5), 2), cex=1.5)
rasterImage(legend_image, 0, 0, 1, 1)
dev.off()

### summary image PhD ###

I <- which(components(G)$membership == 1)
G2 <- induced_subgraph(G, V(G)$name[I])
B2 <- induced_subgraph(LCCB$B, intersect(V(G2)$name[I], V(LCCB$B)$name))
B2 <- add_edges(B2, c("Harry Potter", "Albus Dumbledore"))
B2 <- delete_vertices(B2, "Sirius Black")
E(G2)$color <- rep("darkgrey", length(E(G2)))
E(G2)$color[get.edge.ids(G2, t(ends(B2, E(B2))))] <- adjustcolor("red", 0.5)
edge.width2 <- rep(0.75, length(E(G2)))
edge.width2[get.edge.ids(G2, t(ends(B2, E(B2))))] <- 7
layout2 <- layout[I,]
layout2[V(G2)$name=="Moaning Myrtle", 1] <- layout2[V(G2)$name=="Moaning Myrtle", 1] - 10
labelsI2 <- which(V(G2)$name %in% labels)
charactersInPlot2 <- V(G2)$name
charactersInPlot2[-labelsI2] <- ""
vertex.size2 <- rep(7, length(V(G2)))
vertex.size2[-labelsI2] <- 1
V(G2)$color <- rep(adjustcolor("lightblue", 0.75), length(V(G2)))
V(G2)$color[V(G2)$name %in% V(B2)$name] <- adjustcolor("red", 0.75)
p <- par(mar=c(0, 0, 0, 0) + 0.1)
plot(G2, layout=layout2, edge.width=edge.width2, vertex.size=vertex.size2, vertex.label.cex=1.25, vertex.label.color="black",
     vertex.label.font=2, asp=0, vertex.label=charactersInPlot2); p
