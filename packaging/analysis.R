#' Runs an AIC analysis on a range of cluster sets
#' The AIC obtained is based on a predictive model of cluster growth
#' NOTE: The default additional parameter for this analysis is "Time". This may or may not be a row in inputted cluster.range data
#' The default outcome variable is growth. This generally means that this function is expecting a cluster with annotated growth data.
fit.analysis <- function(cluster.data, mc.cores=1, null.formula=Growth~Size, full.formula=Growth~Size+Time, 
                         predictor.model=function(f, x){glm(formula=f, data = x, family="poisson")},
                         predictor.transformations=list("Time"=function(x){mean(x)})) {
  #' @param cluster.data: Inputted set(s) of cluster data. May or may not be sorted into ranges
  #' @param mc.cores: A parallel option
  #' @param predictor.model: An inputted predictive model function to be used on the data set.
  #' @param full.formula: The full model for the prediction of growth. This will be compared to a null Growth~Size model
  #' @param predictor.transformations: A named list of transformation functions for each predictor variable
  #' @return: A data.table of analysis results. Several important summary values such as null and full AIC are proposed here
  
  #Check inputs
  predictors <- names(predictor.transformations)
  setIDs <- unique(cluster.data[,SetID])
  formula.elements <- unlist(lapply(as.character(full.formula), function(x){strsplit(x, " ")[[1]]}))
  if(!all((formula.elements)%in%c(colnames(cluster.data), "+", "~", "-", "*"))){
    warning("Predictors for formula may not be in the range of cluster data")
  }
  if(!all((predictors)%in%colnames(cluster.data))){
    stop("Predictors referenced in transform step are not in the range of cluster data")
  }
  if(!("RangeID"%in%colnames(cluster.data))){
    warning("No range ID, by default this will be set to 0 for all sets")
    cluster.data[,"RangeID" := 0]
  }
  if(!("Growth"%in%colnames(cluster.data))){
    warning("No Growth information from clusters. By default this will be set to 0 for all sets")
    cluster.data[,"Growth" := 0]
  }
  
  #Transform cluster data for modelling based on inputs
  model.data <- cluster.data[, c("Header", "Size", "Growth", "SetID", "RangeID")]
  model.data[, (predictors) := lapply(predictors, function(x){
    sapply(cluster.data[, get(x)], function(z){(predictor.transformations[[x]])(z)})
  })]
    
  #Obtain fit data for each cluster set
  cluster.analysis <- dplyr::bind_rows(
    parallel::mclapply(setIDs, function(id) {
      print(id)
      DT <- model.data[SetID==id, ]
      suppressWarnings(null.fit <- predictor.model(null.formula, DT))
      suppressWarnings(full.fit <- predictor.model(full.formula, DT))
      
      res <- data.table::data.table("NullFit"=list(null.fit), "FullFit"=list(full.fit), "SetID"=DT[1,SetID], "RangeID"=DT[1,RangeID])
      return(res)
    }, mc.cores=mc.cores))
  
  return(cluster.analysis)
}

