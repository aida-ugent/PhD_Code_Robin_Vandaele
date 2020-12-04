# Code to wrap the backbone_CG function using the dynverse library
# More information can be found on https://dynverse.org/

definition <- dynwrap::definition(
  method = def_method(
    id = "bcb"
  ),
  parameters = def_parameters(
    dynparam::integer_parameter(
      id = "dim",
      default = 20,
      distribution = dynparam::uniform_distribution(2, 100),
      description = "Maximal number of dimensions the data will be reduced to"
    ),
    dynparam::integer_parameter(
      id = "k",
      default = 10,
      distribution = dynparam::uniform_distribution(5, 100),
      description = "The number of neighbors to build a kNN graph"
    ),
    dynparam::character_parameter(
      id = "dimred",
      default = "diffusionmap",
      values = c("pca", "diffusionmap"),
      description = "The type of dimensionality reduction"
    )
  ),
  wrapper = def_wrapper(
    input_required = "expression",
    input_optional = c("groups_network")
  )
)

run_fun <- function(expression, priors, parameters, seed, verbose){

  # perform a dimensionality reduction on the expression matrix (Diffusion Map)
  if(parameters$dimred == "pca"){
    dimred <- irlba::prcomp_irlba(expression, n=min(parameters$dim, nrow(expression) - 1, ncol(expression) - 1))
    ndim <- dim(dimred$x)[2]
    curve <- dimred$sdev[1:ndim]
  }
  else if(parameters$dimred == "diffusionmap"){
    dimred <- diffusionMap::diffuse(calculate_distance(expression, method="pearson"), neigen=parameters$dim, delta=10e-5)
    ndim <- dim(dimred$X)[2]
    curve <- dimred$eigenmult[1:ndim]
  }
  else stop("dimred must be 'pca' or 'diffusionmap'")

  # select optimal number of dimensions if ndim is large enough
  if (ndim > 3){
    # this code is adapted from the expermclust() function in TSCAN
    # the only difference is in how PCA is performed
    # (they specify scale. = TRUE and we leave it as FALSE)
    x <- 1:ndim
    optpoint1 <- which.min(sapply(2:10, function(i){
      x2 <- pmax(0, x - i)
      sum(lm(curve ~ x + x2)$residuals^2 * rep(1:2, each=10))
      }))

    # this is a simple method for finding the "elbow" of a curve, from
    # https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
    x <- cbind(1:ndim, curve)
    line <- x[c(1, nrow(x)),]
    proj <- princurve::project_to_curve(x, line)
    optpoint2 <- which.max(proj$dist_ind) - 1

    # we will take more than 3 components only if both methods recommend it
    optpoint <- max(c(min(c(optpoint1, optpoint2)), 3))
  } else optpoint <- ndim

  if(parameters$dimred == "pca") dimred <- dimred$x[,seq_len(optpoint)]
  else if(parameters$dimred == "diffusionmap") dimred <- dimred$X[,seq_len(optpoint)]
  rownames(dimred) <- rownames(expression)

  # check if number of leaves are given, estimate otherwise
  if(!is.null(priors$groups_network)){
    n_leaves <- sum(table(unlist(priors$groups_network)) == 1)
    if(n_leaves < 2) n_leaves <- NA
  } else n_leaves <- NA

  # compute backbone object
  milestone_network <- backbone(dimred, k=parameters$k, leaves=n_leaves, max_leaves=30)$B

  # convert to dynwrap model
  milestone_ids <- V(milestone_network)$name
  milestone_network <- igraph::as_data_frame(milestone_network)
  colnames(milestone_network)[3] <- "length"
  milestone_network[,"directed"] <- FALSE
  dynwrap::wrap_data(cell_ids=rownames(expression)) %>%
    dynwrap::add_dimred_projection(milestone_ids=milestone_ids,
                                   milestone_network=milestone_network,
                                   dimred=dimred,
                                   dimred_milestones=dimred[milestone_ids,])
}

ti_bcb <- create_ti_method_r(definition, run_fun)
