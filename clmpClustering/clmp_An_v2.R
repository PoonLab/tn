require("ape")
require("phangorn")

#Import Tree Data and annotate with sequence ID and Time
#Sequences must be dated with the date separated from the id by '_'. 
##TO-DO: Currently only accepts year dates. Work to allow more specific dates. 
impTree <-function(iFile){
  #@param iFile: The name/path of the input file (expecting a newick file)
  #@preturn: An ape tree object with associated lists of sequence ID and Time
  
  #Creating an ape tree object from the newick file
  t <- read.tree(iFile)
  
  #Obtain lists of sequence ID and Time
  temp <- sapply(t$tip.label, function(x) strsplit(x, '_')[[1]])
  ids <- temp[1,]
  times <- as.numeric(temp[2,])
  
  #Assign those lists to the tree
  t$Time <- times
  t$ID <- ids
  
  #Summarize internal branch length information
  nodes <- unique(t$edge[,1])
  dist <- sapply(nodes, function(x) {
    if(x%in%t$edge[,2]){
      t$edge.length[which(t$edge[,2]%in%x)]
    }
    else {NA}
  })
    
  t$dSum <- data.frame(NodeID=nodes, Dist=dist)
  
  return(t)
}

#A simple function, removing edges that sit above a maximum reporting distance (@param:maxD).
dFilt <- function(iT, maxD) {
  iT$dSum <- subset(iT$dSum, Dist<maxD)
  return(iT)
}

#Create clusters based on component clustering by some measure of genetic distance
clmpClu <- function(iT) {
  #@param iG: The inputted graph. Expecting all vertices, but some edges filtered by distance.
  #@return: The inputted graph, annotated with a cluster size summary and case membership in the vertices section


  #Obtain the descendants of each node
  decs <- lapply(nodes, function(x){
    l <- Descendants(iT,x,"all")
    l <- setdiff(l, 1:length(iT$ID))
    return(l)
  })
  
  #Obtain only decsendant lists who are all within the list of clusters
  clu <- decs[sapply(decs, function(x){
    (length(x)==length(which(x%in%iT$dSum$NodeID)))&&length(x>1)
  })]
  
  #Collapse all vectors that are subsets of other vectors in list
  clu <- clu[!sapply(seq_along(clu), function(i) max(sapply(clu[-i],function(L) all(clu[[i]] %in% L))))]
  
  clu <- lapply(clu, function(x) {
    chd <- iT$edge[which(iT$edge[,1]%in%x),2]
    chd <- chd[chd%in%1:length(iT$ID)]
  })
  
  return(clu)
}

###############Testing

cutoffs <- seq(0.001,0.02,0.0005)
res <- lapply(cutoffs, function(x){
  subT <- dFilt(t, x)
  c <- clmpClu(subT)
  return(c)
})

cnum <- sapply(res, function(c) {
  clustered <- sapply(c, function(x){length(x)})
  sing <- rep(1,length(t$ID)-sum(clustered))
  cnum <- (length(clustered)+length(sing))
  return(cnum)
})


csize <- sapply(res, function(c) {
  clustered <- sapply(c, function(x){length(x)})
  sing <- rep(1,length(t$ID)-sum(clustered))
  csize <- mean(c(sing,clustered))
  return(csize)
})

par(mfrow=c(1,2))
plot(cutoffs, cnum, xlab="Threshold", ylab="Number of Clusters")
plot(cutoffs, csize, xlab="Threshold", ylab="Mean Cluster Size")

#Import Data
treeFile <- "~/Data/Seattle/analysis/FTStsubB.nwk"
iFile <- treeFile
t <- impTree(treeFile)