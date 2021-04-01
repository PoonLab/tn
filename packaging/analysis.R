
#' Runs an AIC analysis on a range of cluster sets
#' The AIC obtained is based on a predictive model of cluster growth
AIC.analysis <- function(cluster.range, model.formula=Growth~Size+Time, predictor.transformations=list("Time"=function(x){mean(x)})) {
  #' @param cluster.range: An inputted range of cluster sets across multiple parameters
  #' @param model.formula: The full model for the prediction of growth. This will be compared to a null Growth~Size model
  #' @param predictor.transformations: A named list of transformation functions for each predictor variable
  #' @return: A data.table of analysis results. Several important summary values such as cluster size summaries are also reported here
  
}