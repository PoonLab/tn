require("ape")
require("phangorn")

#Obtain times from a properly formatted .fasta file
getT <- function(iFile) {
  
  #Obtain IDs and times of Sequences
  seqs <- read.dna(iFile, format = "fasta", as.character = T)
  IDs <- names(seqs[,1])
  return(as.numeric(sapply(IDs, function(x) (strsplit(x,'_')[[1]])[[2]])))
}

#Obtain times from a properly formatted .fasta file
getIds <- function(iFile) {
  
  #Obtain IDs and times of Sequences
  seqs <- read.dna(iFile, format = "fasta", as.character = T)
  return(names(seqs[,1]))
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
##TO-DO: Currently Unnused (possibly shouldn't ever be used?)
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
  
  oT <- impTree(paste0(oTFile,"/RAxML_result.Tree"))
  oIds <- oT$tip.label
  nIds <- setdiff(getIds(sFile), oIds)
  
  #The nodes associated with new cases. If these nodes are in clusters, we presume that the new cases are in that cluster as well.
  node.assoc <- sapply(trees, function(t){
    nTip <- which(t$tip.label%in%nIds)
    adj <- Siblings(t, nTip)
    
    #If the sibling is a tip, returns the sibling's parent in complete tree
    if(adj<length(t$tip.label)) {
      adj <- which(oT$tip.label%in%t$tip.label[adj])
      parent <- oT$edge[which(oT$edge[,2]==adj),1]
      return(parent) 
    }
    else{
      return(adj)
    }
  })
  
  #Labelling new nodes consistantly and attatch this as growth info
  ntip.label <- sapply(trees, function(t){t$tip.label[which(t$tip.label%in%nIds)]})
  ntip.number <- length(oT$tip.label)+seq(1,length(nIds),1)
  growth <-  data.frame(ntip.label=nIds, ntip.number=ntip.number, node.assoc=as.numeric(node.assoc), stringsAsFactors = F)
  oT$growth <- growth
  
  return(oT)
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
  ids <- unname(temp[1,])
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
  t$nSum <- data.frame(NodeID=nodes, Dist=dist) 
  
  ##TO-DO: Add bootstrap here. Requires bootstrap to be added in tree creation.
  #bootstraps <- as.numeric(t$node.label)
  #t$nSum$Boot <- bootstraps
  
  return(t)
}


##TO-DO: Fuse two filt functions?

#A simple function, removing edges that sit above a maximum ancestral edge length (@param:maxD).
dFilt <- function(iT, maxD) {
  iT$nSum <- subset(iT$nSum, Dist<=maxD)
  return(iT)
}

#A simple function, removing edges that sit under a minimum bootstrap (@param:maxD).
##Curruntly not used
bFilt <- function(iT, maxB) {
  iT$nSum <- subset(iT$nSum, Dist<=maxD)
  return(iT)
}

####


#Cluster using a moidify subtree-based method
STClu <- function(iT) {

  #Obtain a list of node numbers (matches indexes in the tree)
  nodes <- (length(iT$ID)+1):(length(iT$ID)+iT$Nnode) 

  #Obtain the descendants of each node
  decs <- lapply(nodes, function(x){Descendants(iT,x,"all")})
  cluCon <- sapply(decs, function(x){
    intNodes <- setdiff(x, 1:length(iT$ID))
    (length(intNodes)==length(which(intNodes%in%iT$nSum$NodeID)))
  })
  
  #Obtain only descendant lists who are all within the list of clusters
  clu <- decs[cluCon]
  
  #Check the clustering conditions for all possible subtrees and ensure no cherry clusters are added
  filtCon <- sapply(seq_along(clu), function(i){
    subsetCon <- sum(sapply(clu[-i], function(L) {
      all(clu[[i]] %in% L)
    }))>0
    
    #A cherry alone cannot be a cluster
    cherryCon <- length(clu[[i]])==2
    
    return(!(subsetCon||cherryCon))
  })
  
  clu <- clu[filtCon]

  clu <- lapply(clu, function(x){
    
    #Add new cases if they're associated with the internal node
    gNodes <- which(x%in%iT$growth$node.assoc)
    intNodes <- which(x>length(iT$ID))
                      
    growth <- unlist(sapply(x[gNodes], function(i){
      nTip <- which(iT$growth$node.assoc%in%i)
      iT$growth$ntip.number[nTip]
    }))

    x <- x[-intNodes]
    x <- c(x,growth)
    
    return(x)
  })

  iT$clu <- clu
  
  return(iT)
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

iFile <-"~/subT_An_Files/TestML/partial/RAxML_result.Tree" 

oT <- impTree(iFile) 
oT <- growthSim(oTFile, sFile)

#Make temporary trees
makeTemp <- F
if (makeTemp){
  sFile = "test.fasta"
  oSFile <- tempfile("oS", fileext=".fasta")
  oTFile <- tempdir()
  tempTree(sFile, oSFile, oTFile)
}

#Plot Test
plotT <- T
if(plotT){
  
  #Clustering Test
  dists <- oT$nSum$Dist[-1]
  cutoffs <- seq(min(dists),max(dists),(max(dists)-min(dists))/50)
  
  res <- lapply(cutoffs, function(x){
    subT <- dFilt(oT, x)
    subT <- STClu(subT)
    return(subT$clu)
  })
  
  cnum <- sapply(res, function(c) {
    clustered <- sapply(c, function(x){length(x[x<=length(oT$ID)])})
    sing <- rep(1,length(oT$ID)-sum(clustered))
    cnum <- (length(clustered)+length(sing))
    return(cnum)
  })
  
  csize <- sapply(res, function(c) {
    clustered <- sapply(c, function(x){length(x[x<=length(oT$ID)])})
    sing <- rep(1,length(oT$ID)-sum(clustered))
    csize <- mean(c(sing,clustered))
    return(csize)
  })
  
  gmean <- sapply(res, function(c) {
    mean(sapply(c, function(x){length(which(x>length(oT$ID)))}))
  })
  
  gcover <- sapply(res, function(c) {
    sum(sapply(c, function(x){length(which(x>length(oT$ID)))}))/nrow(oT$growth)
  })
  
  par(mfrow=c(2,2))
  plot(cutoffs, cnum, xlab="Threshold", ylab="Number of Clusters")
  lines(cutoffs, cnum)
  plot(cutoffs[-51], csize[-51], xlab="Threshold", ylab="Mean Cluster Size")
  lines(cutoffs[-51], csize[-51])
  plot(cutoffs[-51], gmean[-51], xlab="Threshold", ylab="Mean Cluster Growth")
  lines(cutoffs[-51], gmean[-51])
  plot(cutoffs, gcover, xlab="Threshold", ylab="Total Growth Coverage")
  lines(cutoffs, gcover)
}
