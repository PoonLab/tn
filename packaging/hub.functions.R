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