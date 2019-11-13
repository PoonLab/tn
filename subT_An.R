require("ape")
require("phangorn")

#Obtain times from a properly formatted .fasta file
getT <- function(iFile) {
  
  #Obtain IDs and times of Sequences
  seqs <- read.dna(iFile, format = "fasta", as.character = T)
  IDs <- names(seqs[,1])
  return(as.numeric(sapply(IDs, function(x) (strsplit(x,'_')[[1]])[[2]])))
}

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
  write.dna(seqs[keepInd,], colsep="", oFile, "fasta")
}

#A function to create a temporary tree through FastTree
FastTreeCall <- function(iFile, opt="-nt -gtr -log log.txt ", oFile) {
  
  #Runtext is input
  runText <- paste0("bash -c 'FastTree ", opt, iFile, " > ", oFile,"'")
  system(runText)
}

#A function to create a temporary tree through RaxML
RaxMLCall <- function(iFile, opt="-m GTRCAT -n Tree -p 123 -T 4 -s ", oFile) {
  
  #Runtext is input
  runText <- paste0("bash -c 'raxmlHPC ", opt, iFile, " -w ", oFile,"'")
  system(runText)
}

#Simulate the growth of trees by placing recent sequences as tips on a fixed ML tree
growthSim <- function(oTFile, sFile) {
  
  #Create a temporary reference package location
  unlink(tempdir(), recursive = T)
  tempdir(check=TRUE)
  ref <- tempfile(fileext = ".refpkg")

  #Create reference package using Taxtastic
  runText <- paste0("taxit create -l pol_SubB -P ",  ref, 
                    " --aln-fasta ", sFile,
                    " --tree-stats ", oTFile, "/RAxML_info.Tree",
                    " --tree-file ", oTFile, "/RAxML_result.Tree")
  system(runText)
  
  #Run pplacer on reference package, calling the inputted sequence alignment
  jplace <- tempfile(fileext=".jplace")
  runText <- paste0("bash -c 'pplacer -o ", jplace, " -c ", ref, " ", sFile, "'")
  system(runText)
  
  #Run guppy on pplacer alignment, creating a temporary file of 100 trees (each with 1 of the new cases bound to it)
  oFile <- tempfile(fileext=".tre")
  runText <- paste0("bash -c 'guppy sing --point-mass ", jplace, " -o ", oFile, "'")
  system(runText)
  trees <- read.tree(oFile)
  
  return(trees)
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

#Build a temporary tree from old sequences. Currently uses RaxML for compatibility with pplacer. 
tempTree <- function(sFile="~/Data/Seattle/SeattleB_PRO.fasta", oSFile, oTFile) {
  
  #Obtain Times from sequence file in order to specify the old partition
  times <- getT(sFile)
  keepT <- head(as.numeric(names(table(times))), -1)
  
  tFilt(sFile, keepT, oSFile)
  
  #Create trees, outputted to the temporary tree directory
  RaxMLCall(oSFile, oTFile)
}

###############Testing

#Set up paths to dependencies
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH=paste0(old_path,":~/Desktop/pplacer:~/Desktop/RAxML"))

#For test data
sFile <- "~/subT_An_Files/testS.fasta"
oTFile <- "~/subT_An_Files/TestML/partial"
oSFile <- "~/subT_An_Files/testSF.fasta"

#100 test trees 
trees <- growthSim(oTFile, sFile)

#Make temporary trees
makeTemp <- F
if (makeTemp){
  sFile = "test.fasta"
  oSFile <- tempfile("oS", fileext=".fasta")
  oTFile <- tempdir()
  tempTree(sFile, oSFile, oTFile)
}

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
