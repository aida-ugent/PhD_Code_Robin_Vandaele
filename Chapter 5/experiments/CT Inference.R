# Load libraries
devtools::load_all() # load BCB
library("tidyverse") # read RDS files

# Infer trajectory of example data set

dataset1 <- dynwrap::example_dataset
eval1 <- dyneval::evaluate_ti_method(
  dataset = dataset1,
  method = ti_bcb(),
  parameters = list(dimred="pca"),
  metrics = c("correlation", "him", "featureimp_wcor", "F1_branches")
)
model1 <- eval1$models[[1]]
dynplot::plot_graph(model1)
dynplot::plot_dimred(model1, expression_source=dataset1)
eval1$summary$correlation
eval1$summary$him
eval1$summary$featureimp_wcor
eval1$summary$F1_branches

# Simulate bifurcating trajectory data and infer the trajectory

set.seed(99)
dataset2 <- dyntoy::generate_dataset(model=dyntoy::model_bifurcating(), num_cells=1000)
eval2 <- dyneval::evaluate_ti_method(
  dataset = dataset2,
  method = ti_bcb(),
  parameters = list(),
  metrics = c("correlation", "him", "featureimp_wcor", "F1_branches")
)
model2 <- eval2$models[[1]]
dynplot::plot_graph(model2)
dynplot::plot_dimred(model2, expression_source=dataset2)
eval2$summary$correlation
eval2$summary$him
eval2$summary$featureimp_wcor
eval2$summary$F1_branches

# Infer the trajectory of real gene expression data

dataset3 <- readRDS("Data/Cell/Real/Cell_Data_Real")
eval3 <- dyneval::evaluate_ti_method(
  dataset = dataset3,
  method = ti_bcb(),
  parameters = list(),
  metrics = c("correlation", "him", "featureimp_wcor", "F1_branches"),
  give_priors = c("groups_network")
)
model3 <- eval3$models[[1]]
dynplot::plot_graph(model3)
dynplot::plot_dimred(model3, expression_source=dataset3)
eval3$summary$correlation
eval3$summary$him
eval3$summary$featureimp_wcor
eval3$summary$F1_branches

# Infer the trajectory of second real gene expression data

dataset4 <- readRDS("Data/Cell/Real/cellbench-SC1_luyitian")
eval4 <- dyneval::evaluate_ti_method(
  dataset = dataset4,
  method = ti_bcb(),
  parameters = list(),
  metrics = c("correlation", "him", "featureimp_wcor", "F1_branches")
)
model4 <- eval4$models[[1]]
dynplot::plot_graph(model4)
dynplot::plot_dimred(model4, expression_source=dataset4)
eval4$summary$correlation
eval4$summary$him
eval4$summary$featureimp_wcor
eval4$summary$F1_branches
