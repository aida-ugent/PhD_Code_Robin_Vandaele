#' Compute persistent homology in Python using the ripser library (https://pypi.org/project/ripser/)
#'
#' @param X               A matrix, data frame, or distance object.
#' @param path            The path to the folder containing the python script for computing
#'                        persistent homology.
#' @param maxdim          Maximum homology dimension computed. Will compute all dimensions lower
#'                        than and equal to this value. For 1, H_0 and H_1 will be computed.
#' @param thresh          Maximum distances considered when constructing filtration.
#'                        If infinity, compute the entire filtration.
#' @param coeff           Compute homology with coefficients in the prime field Z/pZ for p=coeff.
#' @param distance_matrix Indicator that X is a distance matrix, if not we compute
#'                        distances in X using the chosen metric.
#' @param do_cocycles     Indicator of whether to compute cocycles, if so, we compute and store
#'                        cocycles in the 'cocycles' member variable
#' @param metric          The metric to use when calculating distance between instances in a
#'                        feature array. If metric is a string, it must be one of the options
#'                        specified in pairwise_distances, including "euclidean", "manhattan",
#'                        or "cosine". Alternatively, if metric is a callable function, it is
#'                        called on each pair of instances (rows) and the resulting value
#'                        recorded. The callable should take two arrays from X as input and
#'                        return a value indicating the distance between them.
#' @param n_perm          The number of points to subsample in a "greedy permutation,"
#'                        or a furthest point sampling of the points.  These points
#'                        will be used in lieu of the full point cloud for a faster
#'                        computation, at the expense of some accuracy, which can
#'                        be bounded as a maximum bottleneck distance to all diagrams
#'                        on the original point set.
#'
#' @return A persistent homology object as computed by the ripser.ripser function in Python
#'         diagrams are post-processed as to be able to plot them using the TDA library in R

# Load libraries

persistencePython <- function(X, path=file.path(getwd(), "Python"), maxdim=1, thresh=Inf, coeff=2, distance_matrix=FALSE,
                        do_cocycles=FALSE, metric="euclidean", n_perm=NULL) {

  # Save arguments to a JSON file
  write(suppressWarnings(RJSONIO::toJSON(list(X=as.matrix(X), path=path, maxdim=1, thresh=thresh, coeff=coeff,
                                     distance_matrix=distance_matrix, do_cocycles=do_cocycles, metric=metric, n_perm=n_perm))),
        file=file.path(path, "arguments.json"))

  # Compute persistent homology in Python
  system("mprof clean")
  script <- file.path(path, "Persistence.py")
  system(paste("mprof run --include-children", script, path), ignore.stdout=TRUE)

  # Read results
  results <- RJSONIO::fromJSON(file.path(path, "results.json"))
  file.remove(file.path(path, "arguments.json"))
  file.remove(file.path(path, "results.json"))
  print(paste("Persistent cohomology computation performed in", round(results$elapsed_time, 3), "seconds"))

  # Read memory consumption in MB (per 0.1 second)
  mprof <- system("mprof list", intern=TRUE)
  results$mprof <- read.delim(substr(mprof, stringr::str_locate(mprof, "mprofile")[,"start"], # stored in current working directory
                                     stringr::str_locate(mprof, ".dat")[,"end"]), skip=1, sep=' ')[,2]
  system("mprof clean")

  # Post-process diagrams
  results$dgms <- do.call("rbind", lapply(1:length(results$dgms), function(i) {
    cbind(i - 1, do.call("rbind", lapply(results$dgms[[i]], function(bar){
      if(is.list(bar)){
        return(c(bar[[1]], Inf))
      }
      else return(bar)
    })))
  }))
  colnames(results$dgms) <- c("dimension", "Birth", "Death")
  class(results$dgms) <- "diagram"

  # Post-process cocycles
  if(!is.null(results$cocycles)){
    results$cocycles <- lapply(results$cocycles, function(ccs) {
      lapply(ccs, function(cc) do.call("rbind", cc) + 1)
    })
  }
  colnames(results$dgms) <- c("dimension", "Birth", "Death")
  class(results$dgms) <- "diagram"

  # Post-process distance matrix
  results$dperm2all <- sapply(results$dperm2all, function(row){sapply(row, function(d) if(is.null(d)) return(Inf) else return(d))})

  # Post-process sample index
  results$idx_perm <- results$idx_perm + 1

  return(results)
}
