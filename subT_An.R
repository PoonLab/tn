require("ape")
require("phangorn")

#A simple function, removing sequences that sit above a maximum time point (@param: maxT)
tFilt <- function(iFile, keepT, oFile) {
  
  #Obtain IDs and times of Sequences
  seqs <- read.dna(iFile, format = "fasta", as.character = T)
  IDs <- names(seqs[,1])
  times <- as.numeric(sapply(IDs, function(x) (strsplit(x,'_')[[1]])[[2]]))
  
  #For an even number of sequences at each year
  #t <- as.numeric(names(table(times)))
  #keepInd <- unlist(lapply(t, function(x){head(which(times==x), 100)}))
  
  #keep only those specified
  keepInd <- which(times%in%keepT) 
  
  #Write only those sequences with ambiguity below 1.5% and sequence length above 85%
  write.dna(seqs[keepInd,], oFile, "fasta")
}

#A function to create a temporary tree (currently uses FastTree)
#TODO: Finalize with tempfile, raxmL and potential skipping (probably skipped at the imp level)
makeTree <- function(iFile, MLfunLoc="FastTree",
                     opt="-nt -gtr -log log.txt ", relation=">", oFile="~/tree.tre") {
  
  #Runtext is input
  runtext <- paste("bash -c '", MLfunLoc, opt, iFile, relation, oFile,"'")
  system(runtext)
  
}

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
    
  #A summary of cluster-relevant node information
  t$nSum <- data.frame(NodeID=nodes, Dist=dist, Boot=as.numeric(t$node.label))
  
  return(t)
}

##TO-DO: Fuse two filt functions?

#A simple function, removing edges that sit above a maximum ancestral edge length (@param:maxD).
dFilt <- function(iT, maxD) {
  iT$nSum <- subset(iT$nSum, Dist<maxD)
  return(iT)
}

#A simple function, removing edges that sit under a minimum bootstrap (@param:maxD).
bFilt <- function(iT, maxB) {
  iT$nSum <- subset(iT$nSum, Dist<maxD)
  return(iT)
}

##TO-DO: Slow
#Create clusters based on component clustering by some measure of genetic distance
STClu <- function(iT) {
  #@param iG: The inputted graph. Expecting all vertices, but some edges filtered by distance.
  #@return: The inputted graph, annotated with a cluster size summary and case membership in the vertices section

  nodes <- (length(iT$ID)+1):(length(iT$ID)+iT$Nnode) 

  #Obtain the descendants of each node
  decs <- lapply(nodes, function(x){
    l <- Descendants(iT,x,"all")
    l <- setdiff(l, 1:length(iT$ID))
    return(l)
  })
  
  #Obtain only descendant lists who are all within the list of clusters
  clu <- decs[sapply(decs, function(x){
    (length(x)==length(which(x%in%iT$nSum$NodeID)))&&length(x>1)
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


#Make Option
makeT <- F
if(makeT) {
  
  #TreeFile Tests (FastTree)
  makeTree("~/Data/Seattle/SeattleB_PRO.fas", oFile="~/cT.nwk")
  #tf <- tempfile("cT", fileext=".nwk")
  #makeTree("~/Data/Seattle/SeattleB_PRO.fas", oFile=tf)
  
  cT <- impTree("~/cT.nwk")
  keepT <- head(as.numeric(names(table(cT$Time))), -1)
  nT <- head(as.numeric(names(table(cT$Time))), 1)

  #tf1 <- tempfile("oT", fileext=".fas")
  tFilt("~/Data/Seattle/SeattleB_PRO.fas", keepT, "~/oT.fas")
  tFilt("~/Data/Seattle/SeattleB_PRO.fas", nT, "~/nT.fas")
  
  makeTree("oT.fas", oFile="~/oT.nwk")
  #tf2 <- tempfile("oT", fileext=".nwk")
  #makeTree(tf1, oFile=tf2)
}

#Import trees from files
oT <- impTree("oT.nwk")
cT <- impTree("cT.nwk")

#Plot Test
plotT <- F
if(plotT){
  
  #Clustering Test
  cutoffs <- seq(0.001,0.02,0.0005)
  res <- lapply(cutoffs, function(x){
    subT <- dFilt(oT, x)
    c <- STclu(subT)
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
}
