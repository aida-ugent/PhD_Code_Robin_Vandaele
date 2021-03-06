% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kNNgraph.R
\name{kNNgraph}
\alias{kNNgraph}
\title{Construct a undirected kNN graph from a distance object by connecting each node
to its k closest neighbors}
\usage{
kNNgraph(X, k = 10)
}
\arguments{
\item{X}{either a data.frame/matrix containing observations,
or a pairwise distance object to be used for constructing the complex}

\item{k}{number of neighbors in the kNN graph, standard is 10}
}
\value{
undirected k-nearest neighbor igraph object
}
\description{
Construct a undirected kNN graph from a distance object by connecting each node
to its k closest neighbors
}
