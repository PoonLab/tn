#' Runs a given clustering method over a range of parameters values.
multi.cluster <- function(cluster.method, param.list, mc.cores=1, verbose=T, rangeID=0) {
  #'@param t: The input tree file, annotated with vertex and edge information
  #'@param param.list: A named list of parameter sets. Each must correspond to the clustering method used. 
  #'@param rangeID: If several different parameter ranges are used, the rangeID can identify them.
  #'@param mc.cores: A parallel option
  #'@param verbose: An output monitoring option
  #'@return: A larger data.table with parameter sets noted
  
  #Cluster method loop
  cluster.range <- parallel::mclapply(1:length(param.list), function(i){
    x <- param.list[[i]]
    x$setID <- i
    if(verbose){
      flush.console()
      print(paste0(i, " of ", length(param.list)))
    }
    do.call(cluster.method, x)
  }, mc.cores=mc.cores)
  
  cluster.range <- dplyr::bind_rows(cluster.range)
  suppressWarnings(cluster.range[,"RangeID" := rangeID])
  
  return(cluster.range)
}

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
      DT <- model.data[SetID==id, ]
      suppressWarnings(null.fit <- predictor.model(null.formula, DT))
      suppressWarnings(full.fit <- predictor.model(full.formula, DT))
      
      res <- data.table::data.table("NullFit"=list(null.fit), "FullFit"=list(full.fit), "SetID"=DT[1,SetID], "RangeID"=DT[1,RangeID])
      return(res)
    }, mc.cores=mc.cores))
  
  return(cluster.analysis)
}

#'Plots the difference in AIC across the proposed and null models
#'This will reach a central optima, as extremes produce low AIC differences
#'Greater negative values mean larger improvement relative to null model
plot.aic.diff <- function(res){
  #'@param res: The result of a fit.analysis() run.
  #'@return: A set of AIC differences.
  
  #Check inputs
  if(!all(c("NullFit", "FullFit")%in%colnames(res))){
    stop("NullFit and FullFit are not names in result output. Ensure that fit.analysis() was run to
         obtain the result plotted here.")
  }
  
  #Get AIC info and create plot 
  null.aic <- sapply(res$NullFit, function(x){x$aic})
  full.aic <- sapply(res$FullFit, function(x){x$aic})
  aic.diff <- full.aic-null.aic
  
  par(mfrow=c(2, 1), mar = c(0,4.2,1,2), cex.lab=1.2)
  plot(x=res$SetID, type="n", ylim=c(0, max(c(null.aic,full.aic))),
       xlab="", ylab="Akaike's Information Criterion", xaxt='n')
  
  #Background
  bg <- par('usr')
  rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='blanchedalmond', border=NA)
  abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
  abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
  abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
  abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
  abline(h=0)
  
  #Plot Profile with minimum located
  abline(v=which.min(aic.diff), lty=3, lwd=2)
  polygon(c(0, res$SetID, max(res$SetID)), c(0, full.aic, 0) , col=rgb(1,0,0,0.4))
  polygon(c(0, res$SetID, max(res$SetID)), c(0, null.aic, 0) , col=rgb(0,1,1,0.4))
  legend("topright", bg="white",
         legend=c("Full Model", "Null Model", "Overlap"), 
         fill=c(rgb(1,0,0,0.4), rgb(0,1,1,0.4), rgb(0,0.65,0.65,0.7)))
  
  #Plot aic.diff
  par(mar=c(5,4.2,1,2))
  plot(x=res$SetID, type="n", ylim=c(min(aic.diff), max(aic.diff)),
       xlab="SetID", ylab="Difference")
  
  #Background
  bg <- par('usr')
  rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='blanchedalmond', border=NA)
  abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
  abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
  abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
  abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
  
  #Plot difference with minimum located
  lines(res$SetID, aic.diff, lwd=1.6, col="orangered")
  points(res$SetID, aic.diff)
  abline(h=0)
  abline(v=which.min(aic.diff), lty=3, lwd=2)
  text(x=max(res$SetID), y=min(aic.diff), adj=c(1,0), cex=1.2,
       paste0("Highest Loss: ", round(min(aic.diff)),"\nat SetID: ", which.min(aic.diff)))

  return(aic.diff)
}
