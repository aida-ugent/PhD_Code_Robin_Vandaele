% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hard_prune.R
\name{hard_prune}
\alias{hard_prune}
\title{Prune all leaves in a graph either a given amount of times or until there are no more than a specified number of leaves}
\usage{
hard_prune(G, times = 1, leaves = NA)
}
\arguments{
\item{G}{graph to be pruned}

\item{times}{the number of times all leaves must be pruned, standard set to 1}

\item{leaves}{(maximum) number of leaves to be kept, 'times' is ignored if provided}
}
\value{
the pruned graph
}
\description{
Prune all leaves in a graph either a given amount of times or until there are no more than a specified number of leaves
}
