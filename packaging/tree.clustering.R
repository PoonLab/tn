#' Clusters are defined as a series of tips diverging from A high confidence common ancestor.
#' This divergence must be done Through a series of short branches. 
#' NOTE: Running this method requires a tree object with growth.info and path.info defined
step.cluster <- function(t, branch.thresh=0.007, boot.thresh=0, setID=0) {
  #'@param t: The input tree file, annotated with vertex and edge information
  #'@param step.thresh: The maximum branch length criterion defining clusters
  #'@param boot.thresh: The minimum bootstrap criterion defining clusters
  #'@param setID: If several different parameter ranges are used, the setID can identify them
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
  
  #Assign Clusters and update membership info for each
  t$node.info[, "Cluster"] <- path.stop["Node",]
  cluster.set <- t$node.info[unique(t$node.info[,Cluster]),]
  
  cluster.des <- lapply(cluster.set[,ID], function(x){which(t$node.info[1:nrow(t$seq.info),Cluster]%in%x)})
  cluster.seq.cols <- colnames(cluster.set)[6:(ncol(cluster.set)-4)]
  cluster.set[, "Descendants" := cluster.des]
  
  #Re-Obtain membership 
  for(nm in cluster.seq.cols){
    cluster.set[,(nm) := lapply(cluster.des, function(x){
      t$seq.info[x ,get(nm)]
    })]
  }
  cluster.set[, "Membership" := lapply(cluster.des, function(x){t$seq.info[x,ID]})]
  cluster.set[, "Size" := sapply(cluster.des, function(x){length(x)})]
  
  t$growth.info[, "Cluster" := t$node.info[t$growth.info$NeighbourNode, Cluster]] 
  t$growth.info[(TermDistance)<=branch.thresh, Cluster := NA]
  
  #Sum bootstrap values within a given cluster
  growth <- t$growth.info[!is.na(Cluster), sum(Bootstrap), by=.(ID, Cluster)]
  growth <- growth[V1>=boot.thresh, Cluster[which.max(V1)], by=.(ID)]
  growth <- table(growth$V1)
  
  #Attach growth info and a set ID to clusters
  cluster.set[, "Growth" := 0]
  cluster.set[ID%in%names(growth), "Growth" := as.numeric(growth)]
  
  cluster.set[, "BranchThresh" := branch.thresh]
  cluster.set[, "BootThresh" := boot.thresh]
  cluster.set[,"setID" := setID]

  return(cluster.set)
}


##- TO-DO: SOLVE MONOPHYLETIC CLUSTER GROWTH IN A SIMPLE WAY -##
#' Clusters as a monophyletic clade under a high-confidence common ancestor.
#' The pairwise patristic distances in this clade must all 
#' NOTE: Running this method requires a tree object with growth.info defined
mono.pat.cluster <- function(t, dist.thresh, boot.thresh=0, dist.criterion="max.patristic.dist", verbose=T, setID=0){
  #'@param t: The input tree file, annotated with vertex and edge information
  #'@param dist.criterion: A particular column in node.info that must be less than a distance threshold
  #'@param dist.thresh: The threshold required for clustering.
  #'@param verbose: An output monitoring option
  #'@param setID: If several different parameter ranges are used, the setID can identify them.
  #'@param boot.thresh: The minimum bootstrap criterion defining clusters
  #'@return: A data table which extends a subset of node.info. This includes growth info
  
  #Input Checking
  if(!is.numeric(dist.thresh) | !is.numeric(boot.thresh) ){
    stop("Clustering criteria must be numeric values")
  }
  if(dist.thresh<0 | boot.thresh<0 | boot.thresh>1){
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
  
  cluster.set <- t$node.info[(Clustered),]

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
  cluster.set[, "Growth" := NA]
  
  cluster.set[, "DistThresh" := dist.thresh]
  cluster.set[, "BootThresh" := boot.thresh]
  cluster.set[,"setID" := setID]
  
  return(cluster.set)
}