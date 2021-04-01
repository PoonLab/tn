#' Clusters are defined as a series of tips diverging from A high confidence common ancestor.
#' This divergence must be done Through a series of short branches. 
#' NOTE: Running this method requires a tree object with growth.info and path.info defined
step.cluster <- function(t, branch.thresh, boot.thresh=0) {
  #'@param t: The input tree file, annotated with vertex and edge information
  #'@param step.thresh: The maximum branch length criterion defining clusters
  #'@param boot.thresh: The minimum bootstrap criterion defining clusters
  #'@return: A data table which extends a subset of node.info. This includes growth info

  #Input Checking
  if(!is.numeric(branch.thresh) | !is.numeric(boot.thresh) ){
    stop("Clustering criteria must be numeric values")
  }
  if(branch.thresh<0 | boot.thresh<0 | boot.thresh>1){
    warning("Bootstrap criterion should be from 0-1 and 
            branch length criterion should be positive.")
  }
  if(!("path.info"%in%names(t))){
    stop("path.info must be defined for tree")
  }
  if(!("growth.info"%in%names(t))){
    stop("growth.info must be defined for tree")
  }
  
  #Obtain the stopping point in the path based on branch.thresh
  path.stop <- sapply(t$path.info, function(p) {
    h <- which(p["BranchLength",]>branch.thresh)[1]
    c(p[, h], h)
  })
  rownames(path.stop)[4] <- "Height"
  path.stop["Node", is.na(path.stop["Node",])] <- length(t$tip.label)+1
  
  #Check bootstrap requirements, stepping back down clustered paths until they're met.
  i <- which(path.stop["Boot",]<boot.thresh)
  if(length(i)>0) {
    path.stop[,i]  <- sapply(i, function(j){
      p <- t$path.info[[j]]
      p.boots <- p["Boot", 1:path.stop["Height", j]]
      new.h <- which(p.boots>=boot.thresh)[1]
      return(c(p[, new.h], new.h))
    })
  }
  
  #Assign Clusters
  t$node.info[, "Cluster"] <- path.stop["Node",]
  cluster.info <- t$node.info[unique(t$node.info[,Cluster]),]
  t$growth.info[, "Cluster" := t$node.info[t$growth.info$NeighbourNode, Cluster]] 
  t$growth.info[(TermDistance)<=branch.thresh, Cluster := NA]
  
  #Sum bootstrap values within a given cluster
  growth <- t$growth.info[!is.na(Cluster), sum(Bootstrap), by=.(ID, Cluster)]
  growth <- growth[V1>=boot.thresh, Cluster[which.max(V1)], by=.(ID)]
  growth <- table(growth$V1)
  
  #Attach growth info and a set ID to clusters
  cluster.info[, "Growth" := 0]
  cluster.info[ID%in%names(growth), "Growth" := as.numeric(growth)]
  
  cluster.info[, "BranchThresh" := branch.thresh]
  cluster.info[, "BootThresh" := boot.thresh]

  return(cluster.info)
}


##- TO-DO: SOLVE MONOPHYLETIC CLUSTER GROWTH IN A SIMPLE WAY -##
#' Clusters as a monophyletic clade under a high-confidence common ancestor.
#' The pairwise patristic distances in this clade must all 
#' NOTE: Running this method requires a tree object with growth.info defined
mono.pat.cluster <- function(t, dist.thresh, boot.thresh=0, dist.criterion="max.patristic.dist", verbose=T){
  #'@param t: The input tree file, annotated with vertex and edge information
  #'@param dist.criterion: A particular column in node.info that must be less than a distance threshold
  #'@param dist.thresh: The threshold required for clustering.
  #'@param verbose: An output monitoring option
  #'@param boot.thresh: The minimum bootstrap criterion defining clusters
  #'@return: A data table which extends a subset of node.info. This includes growth info
  
  #Input Checking
  if(!is.numeric(branch.thresh) | !is.numeric(boot.thresh) ){
    stop("Clustering criteria must be numeric values")
  }
  if(branch.thresh<0 | boot.thresh<0 | boot.thresh>1){
    warning("Bootstrap criterion should be from 0-1 and 
            branch length criterion should be positive.")
  }
  if(!("growth.info"%in%names(t))){
    stop("growth.info must be defined for tree")
  }
  
  #Cluster Criterion checking
  t$node.info[, "Clustered" := F]
  t$node.info[(Bootstrap>=boot.thresh)&(get(dist.criterion)<=dist.thresh), "Clustered" := T]
  clustered.des <- unlist(t$node.info[(Clustered), Descendants])
  
  times.clustered <- t$node.info[(Clustered), length(which(clustered.des%in%ID)), by=which(Clustered)]
  sub.clusters <- times.clustered[((which<=nrow(t$seq.info))&(V1>1)) | ((which>nrow(t$seq.info))&(V1>0)), which]
  t$node.info[sub.clusters, "Clustered" := F]  
  
  cluster.info <- t$node.info[(Clustered),]

  #Find parent clusters
  t$node.info[,"Cluster":=0]
  for(i in which(t$node.info[,Clustered])) {
      t$node.info[i, "Cluster" := i]
      t$node.info[t$node.info$Descendants[[i]], "Cluster" := i] 
  }
  
  ##- TO-DO: SOLVE MONOPHYLETIC CLUSTER GROWTH IN A SIMPLE WAY -##
  if(verbose){
    warning("Method unfinished. Growth information for clusters not included")
  }
    
  #Attach growth info and a set ID to clusters
  cluster.info[, "Growth" := NA]
  
  cluster.info[, "BranchThresh" := branch.thresh]
  cluster.info[, "BootThresh" := boot.thresh]
  
  return(cluster.info)
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
    if(verbose){
      flush.console()
      print(paste0(i, " of ", length(param.list)))
    }
    do.call(cluster.method, param.list[[i]])
  }, mc.cores=mc.cores)
  
  cluster.range <- dplyr::bind_rows(cluster.range)
  suppressWarnings(cluster.range[,"rangeID" := rangeID])
  
  return(cluster.range)
}