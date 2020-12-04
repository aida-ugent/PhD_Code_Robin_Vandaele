# Load libraries

library("igraph")
library("tidyverse")
library("dynbenchmark")
library("gridExtra")
library("dplyr")
library("gdata")

# Make design for experiments

metrics <- c("correlation", "him", "featureimp_wcor", "F1_branches")
ks <- c(5, 10, 15, 20, 25)
dimreds <- c("pca", "diffusionmap")
methods <- list("slingshot", "paga")
methods_with_priors <- character()
for(k in ks){
  for(dimred in dimreds){
    methods[[length(methods) + 1]] <- BCB:::ti_bcb(k=k, dim=dim, dimred=dimred)
    methods[[length(methods)]]$method$id <- paste("bcb-", "k", k, "-", dimred, sep="")
    methods[[length(methods) + 1]] <- BCB:::ti_bcb(k=k, dim=dim, dimred=dimred)
    methods[[length(methods)]]$method$id <- paste("bcb-", "k", k, "-", dimred, "-prior", sep="")
    methods_with_priors <- c(methods_with_priors, methods[[length(methods)]]$method$id)
  }
}
design <-  dynbenchmark::benchmark_generate_design(
  dynbenchmark::list_datasets()$id,
  methods=methods,
  priors=list(
    list(id="leaves",
         method_id=methods_with_priors,
         set=c("groups_network")),
    list(id="NA",
         method_id=c("slingshot", "paga")),
    list(id="leaves",
         method_id=c("slingshot", "paga"),
         set=c("start_id"))
  ),
  num_repeats = 1
)

# Submit benchmark experiment to cluster

benchmark_submit(
  design = design,
  metrics = metrics,
  qsub_grouping = "{method_id}/{prior_id}",
  qsub_params =
    list(
      timeout = 10800,
      memory = "15G"
    ),
  local_output_folder="Data/Dynverse/results",
  remote_output_folder="/scratch/irc/shared/bcb/"
)

# Download results from cluster locally

benchmark_fetch_results(local_output_folder="Data/Dynverse/results")

# Bind results to view in R

output <- benchmark_bind_results(local_output_folder="Data/Dynverse/results", load_models=FALSE) # set load_models=TRUE for confusion matrices

# Extract subset of experiments for analysis

topologies <- lapply(load_datasets()$prior_information, function(p) graph.data.frame(p$groups_network))

ks <- c(5, 10, 15, 20, 25)
dimreds <- c("diffusionmap")
priors <- c(FALSE, TRUE)
methods <- "slingshot"
for(k in ks){
  for(dimred in dimreds){
    for(prior in priors){
      methods <- c(methods, paste("bcb-", "k", k, "-", dimred, if(prior) "-prior" else "", sep=""))
    }
  }
}
methodsI <- which(output[["method_id"]] %in% methods)
errorI <- which(output[["dataset_id"]] %in% unique(output[["dataset_id"]][which(output[["error_message"]] != "")]))
topologiesI <- which(output[["dataset_id"]] %in% output[["dataset_id"]][
  which(sapply(topologies, function(t) TRUE))
  ])
outputI <- intersect(intersect(methodsI, setdiff(1:length(output[[1]]), errorI)), topologiesI)
dataset_id <- intersect(1:length(unique(output[["dataset_id"]])), which(output[["dataset_id"]] %in% output[["dataset_id"]][outputI]))
methods <- unique(output[["method_id"]][outputI])
metrics <- c("correlation", "him", "featureimp_wcor", "F1_branches", "time_method", "max_mem_gb")

# Check number of real and synthetic datasets

sum(startsWith(output[["dataset_id"]][dataset_id], "real"))
sum(startsWith(output[["dataset_id"]][dataset_id], "synthetic"))

# Metrics - cumulative plots

cumplots <- list()
for(metric in metrics){
  results <- data.frame(dataset_id=dataset_id, metric=output[[metric]][outputI],
                        method_id=output$method_id[outputI])
  if(metric == "time_method") results[is.na(results[,"metric"]), "metric"] <- Inf
  for(id in methods){
    I <- results[,"method_id"] == id
    results[I, "metric"] <- sort(results[I, "metric"])
    results[I, "dataset_id"] <- seq(sum(I))
  }
  with_prior <- endsWith(as.character(results$method_id), "prior")
  if(sum(with_prior) > 0 & sum(with_prior) < nrow(results)){
    results$method_id <- as.character(results$method_id)
    results$method_id[with_prior] <- substr(results$method_id[with_prior], 1,
                                            nchar(results$method_id[with_prior]) - 6)
    results$method_id <- factor(results$method_id)
    results$prior <- factor(with_prior)
    cumplots[[length(cumplots) + 1]] <- ggplot(results, aes(dataset_id, metric)) +
      geom_line(aes(col=method_id, linetype=prior)) +
      ggtitle(metric) +
      theme_bw() +
      xlab("") +
      theme(plot.title=element_text(hjust=0.5), text = element_text(size=20))
  }
  else cumplots[[length(cumplots) + 1]] <- ggplot(results, aes(dataset_id, metric)) +
    geom_line(aes(col=method_id)) +
    ggtitle(metric) +
    theme_bw() +
    xlab("") +
    theme(plot.title=element_text(hjust=0.5), text = element_text(size=20))
}
ggpubr::ggarrange(plotlist=cumplots, ncol=2, nrow=3, common.legend=TRUE)

# Number of cells - distribution

ncells <- sapply(load_datasets()$cell_ids, function(ci) length(ci)) # number of cells per dataset
hist(ncells[dataset_id], xlab="Number of cells")

# Metrics - plots according to number of cells

ncellplots <- list()
for(metric in metrics){
  results <- data.frame(metric=output[[metric]][outputI], method_id=output$method_id[outputI], ncells=ncells[dataset_id])
  ncellplots[[length(ncellplots) + 1]] <- ggplot(results, aes(ncells, metric)) +
    geom_line(aes(col=method_id)) +
    ggtitle(metric) +
    theme_bw() +
    xlab("number of cells") +
    theme(plot.title=element_text(hjust=0.5))
}
do.call("grid.arrange", c(ncellplots, ncol=2))

# Number of features - distribution

nfeatures <- sapply(load_datasets()$feature_info, function(fi) nrow(fi)) # number of genes per dataset
hist(nfeatures[dataset_id], xlab="Number of features")

# Metrics - plots according to number of cells

nfeatureplots <- list()
for(metric in metrics){
  results <- data.frame(metric=output[[metric]][outputI], method_id=output$method_id[outputI], nfeatures=nfeatures[dataset_id])
  nfeatureplots[[length(nfeatureplots) + 1]] <- ggplot(results, aes(nfeatures, metric)) +
    geom_line(aes(col=method_id)) +
    ggtitle(metric) +
    theme_bw() +
    xlab("number of features") +
    theme(plot.title=element_text(hjust=0.5))
}
do.call("grid.arrange", c(nfeatureplots, ncol=2))

# Topology - distribution

true_leaves <- sapply(topologies, function(t) sum(degree(t) == 1))[dataset_id] # leaves in topology
table_true_leaves <- table(true_leaves)
plot(table_true_leaves, xlab="(True) Leaves", ylab="Frequency")

# Topology - confusion matrices

output <- benchmark_bind_results(local_output_folder="Data/Dynverse/results", load_models=TRUE)
topologies <- lapply(load_datasets()$prior_information, function(p) graph.data.frame(p$groups_network))
ks <- c(5, 10, 15, 20, 25)
dimreds <- c("diffusionmap")
priors <- c(FALSE)
methods <- "slingshot"
for(k in ks){
  for(dimred in dimreds){
    for(prior in priors){
      methods <- c(methods, paste("bcb-", "k", k, "-", dimred, if(prior) "-prior" else "", sep=""))
    }
  }
}
methodsI <- which(output[["method_id"]] %in% methods)
errorI <- which(output[["dataset_id"]] %in% unique(output[["dataset_id"]][which(output[["error_message"]] != "")]))
topologiesI <- which(output[["dataset_id"]] %in% output[["dataset_id"]][
  which(sapply(topologies, function(t) TRUE))
  ])
outputI <- intersect(intersect(methodsI, setdiff(1:length(output[[1]]), errorI)), topologiesI)
dataset_id <- intersect(1:length(unique(output[["dataset_id"]])), which(output[["dataset_id"]] %in% output[["dataset_id"]][outputI]))
methods <- unique(output[["method_id"]][outputI])

leaveplots <- list()
for(method in methods){
  I <- output$method_id == method
  pred_leaves <- sapply(output$model[I][dataset_id], function(m) sum(table(unlist(m$milestone_network[,c(1, 2)])) == 1))
  too_many_leaves <- which(pred_leaves > max(true_leaves))
  print(paste(length(too_many_leaves), "leaves were predicted above maximal possible leaves for method", method))
  if(length(too_many_leaves > 0)){
    print("Predicted leaves discared:")
    print(sort(unique(pred_leaves[too_many_leaves])))
    pred_leaves <- pred_leaves[-too_many_leaves]
    this_true_leaves <- true_leaves[-too_many_leaves]
  }
  else this_true_leaves <- true_leaves
  possible_leaves <- sort(union(unique(this_true_leaves) , unique(pred_leaves)))
  data_leaves <- data.frame(id=1:length(pred_leaves), this_true_leaves=this_true_leaves, pred_leaves=pred_leaves)
  data_leaves <- rbind(as.data.frame(aggregate(id~., data_leaves, length)),
                       cbind(setdiff(data.frame(expand.grid(this_true_leaves=possible_leaves, pred_leaves=possible_leaves)),
                                     data_leaves[,2:3]), id=0))
  data_leaves[,3] <- as.numeric(data_leaves[,3] / table(this_true_leaves)[as.character(data_leaves[,1])])
  data_leaves[is.na(data_leaves[,3]), 3] <- 0
  leaveplots[[length(leaveplots) + 1]] <- ggplot(data_leaves, aes(x=factor(pred_leaves), y=factor(this_true_leaves))) +
    geom_tile(aes(fill=id)) +
    geom_text(aes(label=round(id, 1)), size=5) +
    scale_x_discrete(name="Predicted Leaves") +
    scale_y_discrete(name="True Leaves", limits=as.character(rev(possible_leaves))) +
    scale_fill_gradient2(mid="red", high="green") +
    theme(plot.title=element_text(hjust=0.5), legend.title=element_blank(), text=element_text(size=7)) +
    ggtitle(method) +
    coord_fixed() +
    theme(text = element_text(size=20))
}
do.call("grid.arrange", c(leaveplots, ncol=2))
