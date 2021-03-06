% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prune.R
\name{prune}
\alias{prune}
\title{Prune a graph a given number of times by choosing the number of times/leaves and used method}
\usage{
prune(G, times = 1, leaves = NA, method = "hard", weight = NULL,
  vcost = NULL)
}
\arguments{
\item{G}{graph to be pruned}

\item{times}{times the graph must be pruned, only used if 'method' is 'hard' or 'soft', standard set to 1
may also be a vector with size the number of components of 'the (constructed) graph 'G'}

\item{leaves}{(maximum) number of leaves to be included in the subtree, standard set to NA, 'times' is ignored if provided
may also be a vector with size the number of components of 'the (constructed) graph 'G'}

\item{method}{used method to perform the pruning, standard set to 'hard':
- 'hard': iteratively prune all leaves in the graph
- 'soft': iteratively prune the leaves representing the fewest nodes
- 'search': search for a set of leaves resulting in a high cost tree optimal if G itself is a tree)}

\item{weight}{optional reweighting vector of the edges of G (used only with search prune)}

\item{vcost}{optional vertex cost for which the sum is to be optimized (used only with search prune, ignores weight if provided)}
}
\value{
A list containing the following items
               - B: the pruned graph derived from 'G'
               - membership: the membership of the 'X' points/nodes according to the components of the (constructed) graph
               - (if method is 'search') cost: a data.frame corresponding to the obtained cost according to the number of leaves and component
               - (if method is 'search') full_cost: a vector of the full cost of each component
               - (if method is 'search') includedE: a list of lists (one for each connected component) containing the added edges of 'G' at each iteration
}
\description{
Prune a graph a given number of times by choosing the number of times/leaves and used method
}
