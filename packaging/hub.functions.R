#'Translates a set of sequence headers into a data.frame object for data input.
#'NOTE: This must contain, at minimum, a set of unique sequence id's (labelled ID) and one other variable
#'NOTE: A certain subset of sequences
pull.headers <- function(seqs, var.names, var.transformations=list(), sep="_") {
  #'@param seqs: An inputted alignment using ape's sequence handling
  #'@param var.names: The names of the variables represented in each header. This must contain "ID".
  #'@param var.transformations: A list of transformation functions (such as as.character()) 
  #'these transform each row into it's proper type. by default, each type is set to character.
  #'@param sep: The separator character that splits upthe headers in the fasta file
  #'@return: A data.table object containing the information associated with each sequence.
  
  #Checking Inputs
  if(length(var.names)!=length(unique(var.names))){
    stop("var.names may not contain repeats")
  }
  if(length(var.transformations)==0){
    var.transformations <- lapply(1:length(var.names), function(x){as.character})
  } else {
    if(length(var.names)!=length(var.transformations)) {
      stop("var.names and var.transformations must be equal lengths")
      return(NULL)
    }
  }
  if(("Header"%in%var.names)){
    warning("'Header' is contained within var.names, this will be overwritten")
  }
  
  
  #Split and transform data from headers
  split.headers <- sapply(names(seqs), function(x) {strsplit(x, sep)[[1]]})
  seq.info <- lapply(1:nrow(split.headers), function(i) {
    x <- unname(split.headers[i,])
    x <- var.transformations[[i]](x)
    data.table::data.table(x)
  })
  
  seq.info <- dplyr::bind_cols(seq.info)
  colnames(seq.info) <- var.names
  
  seq.info[,"Header" := names(seqs)]
  
  return(seq.info)
}

#Annotate a subset of the data as "New". This creates an additional column in seq.info
annotate.new <- function(seq.info, which.new=logical(0)){
  #'@param seq.info: A data frame or data.table object containing the sequence meta data
  #'@param which.new: A set of indices of which sequences were to be labelled "new"
  #'@return: seq.info annotated with a true/false "New" column
  
  #Check inputs
  if("New" %in% colnames(seq.info)){
    warning("A column in seq.info is already labelled as New. This will be overwritten")
  }
  if(length(which.new)==0){
    warning("No new sequences are specified. By default, a variable labelled Time will be used to specify new cases.
            This would label cases within the newest year of the time range as new")
    if("Time"%in%colnames(seq.info)){
      new.year <- max(seq.info$Time)-365
      which.new <- which(seq.info$Time>new.year)
    }else{
      stop("No time variable clarified.")
    }
  }
  
  #Annotate New
  seq.info[,"New" := F]
  seq.info[which.new ,"New" := T]
  
  return(seq.info)
}

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