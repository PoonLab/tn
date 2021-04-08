#'Create clusters based on the components of a graph
component.cluster <- function(g, dist.thresh=0.007, setID=0) {
  #'@param g: The input graph, annotated with vertex and edge information
  #'@param dist.thresh: The maximum distance defining which edges are filtered
  #'@param setID: If several different parameter ranges are used, the setID can identify them
  #'@return: A data table which represents cluster information. This includes growth info
  #'Because data.tables are being used, this prevents original values being reassigned via pointer
  
  #Filter edges above the distance threshold and prepare for component finding algorithm
  #All edges from a new sequence are filtered except for their "growth-resolved" edge
  filtered.edges <- g$edge.info<=dist.thresh
  sum(filtered.edges)
  filtered.edges[which(g$seq.info$New),] <- F
  sum(filtered.edges)
  filtered.edges[g$growth.resolved$NewID, g$growth.resolved$Neighbour] <- 
    g$edge.info[g$growth.resolved$NewID, g$growth.resolved$Neighbour]<=dist.thresh
  sum(filtered.edges)
  
  #Run homogenization algorithm to label sequences with their cluster
  seq.cols <- colnames(g$seq.info)
  previous.cluster <- rep(0, nrow(seq.info))
  g$seq.info[, "Cluster" := 1:nrow(seq.info)]
  
  while(all(g$seq.info$Cluster!=previous.cluster)){
    previous.cluster <- g$seq.info$Cluster
    g$seq.info[, Cluster := sapply(1:nrow(g$seq.info), function(i) {
      x <- g$seq.info[which(filtered.edges[i,-i]), Cluster]
      if(length(x)==0){
        return(i)
      }else{
        return(min(x))
      }
    })]
  }
  
  cluster.set <- g$seq.info[!(New),lapply(seq.cols, function(nm){list(get(nm))}),by=Cluster]
  cluster.set[,"Size":=length(V1[[1]]), by=1:nrow(cluster.set)]
  colnames(cluster.set) <- c("ClusterID", seq.cols, "Size")
  cluster.set$New <- NULL
  
  #Attach growth info and set ID
  growth <- table(g$seq.info[(New)&(Cluster%in%cluster.set$ClusterID), Cluster])
  cluster.set[, "Growth" := 0]
  cluster.set[ClusterID%in%as.numeric(names(growth)), Growth := as.numeric(growth)]
  
  cluster.set[,"DistThresh" := dist.thresh]
  cluster.set[,"SetID" := setID]
  
  return(cluster.set)
}