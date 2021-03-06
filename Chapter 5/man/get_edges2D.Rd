% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_edges2D.R
\name{get_edges2D}
\alias{get_edges2D}
\title{Given a two-dimensional dataframe df and a graph G = (V, E), where E is a list of edges
between nodes of df marked by their indexes, return a four-dimensional dataframe containing
the coordinates of the endpoints of each segment defined by df and E, suitable to plot with
geom_segment (ggplot2)}
\usage{
get_edges2D(df, G)
}
\arguments{
\item{df}{two-dimensional dataframe containing the points that must be connected through
segments}

\item{G}{Either a graph containing the edges (u, v) marking a segments between the u-th and
v-th point of df, or a list of edges i}
}
\value{
a four-dimensional dataframe containing the coordinates of the endpoints of the
          edges to plot
}
\description{
Given a two-dimensional dataframe df and a graph G = (V, E), where E is a list of edges
between nodes of df marked by their indexes, return a four-dimensional dataframe containing
the coordinates of the endpoints of each segment defined by df and E, suitable to plot with
geom_segment (ggplot2)
}
