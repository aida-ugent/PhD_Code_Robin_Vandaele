trim_deg2 <- function(g) {
  get_deg2 <- function(x) {
    dd <- degree(x)
    trim <- V(x)[names(dd[dd==2])]
  }
  ng <- g
  trim <- get_deg2(ng)
  while(length(trim)) {
    tv <- trim[1]
    touch <- adjacent_vertices(ng, tv)[[1]]
    ng <- delete_edges(ng, E(ng)[tv %--% touch])
    ng <- add_edges(ng, touch$name)
    ng <- delete_vertices(ng, V(ng)[tv])
    trim <- get_deg2(ng)
  }
  ng
}
