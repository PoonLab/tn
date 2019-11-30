require("ape")
require("phangorn")
library("dplyr",verbose = FALSE)

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
    if(adj<=length(t$tip.label)) {
      adj <- which(oT$tip.label%in%t$tip.label[adj])
      parent <- oT$edge[which(oT$edge[,2]==adj),1]
      return(parent) 
    }
    else{
      return(adj)
    }
  })
  
  #Take in the distance calculated by 
  dist <- sapply(trees, function(t){
    nTip <- which(t$tip.label%in%nIds)
    dist <- t$edge.length[which(t$edge[,2]==nTip)]
    return(dist)
  })
  
  #Labelling new nodes consistantly and attatch this as growth info
  ntip.label <- sapply(trees, function(t){t$tip.label[which(t$tip.label%in%nIds)]})
  ntip.number <- length(oT$tip.label)+seq(1,length(nIds),1)
  growth <-  data.frame(ntip.label=nIds, ntip.number=as.numeric(ntip.number), node.assoc=as.numeric(node.assoc), distance=dist, stringsAsFactors = F)
  oT$gSum <- growth
  
  return(oT)
}

#Import Tree Data and annotate with sequence ID and Time
#Sequences must be dated with the date separated from the id by '_'. 
##TO-DO: Currently only accepts year dates. Work to allow more specific dates. 
impTree <-function(iFile, terminal=F){
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
  tips <- 1:length(t$ID)
  
  #Summarize internal branch length information
  t$dist <- dist.nodes(t)
  t$dist <- t$dist[tips, tips]
  t$tDiff <- sapply(1:length(t$ID), function(i){
    t1 <- t$Time[i]
    sapply(1:length(t$ID), function(j){ t1 - t$Time[j]})
  })
  
  #In this case, we set terminal branch lengths to 0. Meaning that they are not a factor for growth or clustering
  if(!terminal){
    parents <- unique(t$edge[which(t$edge[,2]%in%tips),1])
  }
  
  #Obtain a list of node numbers (matches indexes in the tree)
  nodes <- (length(t$ID)+1):(length(t$ID)+t$Nnode) 
  tipD <- t$dist[1:length(t$ID),1:length(t$ID)] #Subset of just terminal tip to tip distances
  
  #Obtain the descendants of each node
  decs <- lapply(nodes, function(x){Descendants(t,x,"all")})
  meanDist <- sapply(decs, function(x) {
    tips <- x[which(x<=length(t$ID))] 
    mean <- mean(sapply(tips, function(tip){mean(tipD[tip,tips])}))
    return(mean)
  })
  
  #Set up a node-summary data frame
  t$nSum <- data.frame(nID=nodes, mDist=meanDist)
  
  ##TO-DO: Add bootstrap here. Requires bootstrap to be added in tree creation.
  #bootstraps <- as.numeric(t$node.label)
  #t$nSum$Boot <- bootstraps
  
  return(t)
}

#A simple function, removing edges that sit above a maximum ancestral edge length (@param:maxD).
dFilt <- function(iT, maxD) {
  
  #Remove nodes with a greater mean distance than maxD from consideration.
  iT$nSum <- subset(iT$nSum, mDist<=maxD)
  
  #New cases above the max distance from their associated node are considered not to associate with that node
  if(length(which(iT$gSum$distance>=maxD))>0){
    iT$gSum[which(iT$gSum$distance>=maxD),]$node.assoc <- 0 
  }

  return(iT)
}

#Cluster using a modified subtree-based method
STClu <- function(iT) {

  #Obtain a list of node numbers (matches indexes in the tree)
  nodes <- iT$nSum$nID

  #Obtain the descendants of potential clustering node
  decs <- lapply(nodes, function(x){Descendants(iT,x,"all")})
  
  #Check the clustering conditions for all possible subtrees and ensure no cherry clusters are added
  filt <- sapply(seq_along(decs), function(i){
    
    #A cluster cannot be the subset of a larger cluster
    if(length(decs)>1){
      subsetCon <- sum(sapply(decs[-i], function(L) {all(decs[[i]] %in% L)}))>0
    }
    else {
      subsetCon <- F
    }

    #A cherry alone cannot be a cluster
    cherryCon <- length(decs[[i]])==2
    
    return((subsetCon||cherryCon))
  })
  
  #Remove non-clusters from the list
  decs[filt] <- NULL

  #Remove internal nodes for summary. Only cases are of interest
  clu <- lapply(decs, function(x){
    
    #Add new cases if they're associated with the internal and if 
    gNodes <- subset(iT$gSum, node.assoc%in%x)$ntip.number
    intNodes <- which(x>length(iT$ID))
                      
    x <- x[-intNodes]
    x <- c(x,gNodes)
    
    return(x)
  })
  
  #Add singletons and summarise clusters
  sing <- as.list(1:(length(iT$tip.label)+nrow(iT$gSum)))
  cTips <- unlist(clu)
  filt <- sapply(sing, function(tip){tip%in%cTips})
  sing[filt] <- NULL
  
  iT$clu <- c(clu,sing)
  
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

#Obtains some likelihood data in order to weight cases based on their recency 
likData <- function(iT, maxD){
  
  #Obtain some summarized information from the sub-Graph
  vTab <- table(iT$Time)
  times <- as.numeric(names(vTab))
  eTab <- sapply(times, function(t){
    tTips <- which(iT$Time==t)
    tDist <- unlist(unname(lapply(tTips, function(x){iT$dist[x, -x]})))
    length(which(tDist<=maxD))
  })
  names(eTab) <- names(vTab)
  
  #Obtain the minimum retrospective distance for each member sequence 
  ##TO-DO: Optimize, this should only happen once. It may belong in the tree-creation call ##
  minE <- bind_rows(sapply(1:length(iT$ID), function(i){
    retro <- which(iT$tDiff[i,]<0)
    row <- iT$dist[i,retro]
    
    #For the earliest year
    if(length(row)==0){
      NULL
    }
    else {
      close <- which(row==min(row))[[1]]
      data.frame(Tip=i, Distance=row[close], tDiff=abs(iT$tDiff[i,retro[close]]))  
    }
  }))
  
  #Take the total edge frequency data from the graph and format this information into successes and attempts
  #An edge to the newest year falling below the max distance is considered a success
  ageD <- bind_rows(lapply(times, function(t) {
    temp <- subset(minE, Tip%in%which(iT$Time==t))
    dfs <- split(temp, temp$tDiff)
    
    Positive <- sapply(dfs, function(df){length(which(df$Distance<=maxD))})
    vTotal <- sapply(dfs, function(df){vTab[as.character(t)]})
    tDiff <- as.numeric(names(Positive))
    oeDens <- sapply(tDiff, function(tD){
      oTime <- t-tD
      return(eTab[as.character(oTime)]/vTab[as.character(oTime)])
    })
    
    
    res <- data.frame(Positive=as.numeric(Positive), vTotal=as.numeric(vTotal), oeDens=oeDens, tDiff)
    return(res)
  }))
  

  return(ageD)

}

#Obtain GAIC at several different cutoffs
GAICRun <- function(iT) {
  
  #Clustering Test
  dists <- iT$nSum$mDist
  cutoffs <- seq(min(dists),max(dists),(max(dists)-min(dists))/50)
  
  res <- sapply(cutoffs, function(x){
    print(x)
    subT <- dFilt(iT, x)
    subT <- STClu(subT)
    ageD <- likData(subT, x)
    
    #Obtain a model of case connection frequency to new cases as predicted by individual case age
    #Use this to weight cases by age
    mod <- glm(cbind(Positive, vTotal) ~ tDiff+oeDens, data=ageD, family='binomial')
    subT$Weight <- predict(mod, type='response',
                             data.frame(tDiff=max(subT$Time)+1-subT$Time, 
                                        oeDens=as.numeric(eTab[as.character(subT$Time)]/vTab[as.character(subT$Time)])))
    
    #Filter clusters such that new singletons are not considered.
    clu <- subT$clu
    filt <- sapply(clu, function(c){length(which(c<=length(subT$ID)))==0})
    clu[filt] <- NULL  
    
    #Create clusters for this subgraph and measure growth
    cGrowth <- sapply(clu, function(c) {
      length(which(c>length(subT$ID)))
    })
    cPred <- sapply(clu, function(c) {
      sum(subT$Weight[c[which(c<=length(subT$ID))]])
    })
    
    cNPred <- sapply(clu, function(c) {
      length(c[which(c<=length(subT$ID))])
    })
    
    #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
    df1 <- data.frame(Growth = cGrowth, Pred = cPred)
    df2 <- data.frame(Growth = cGrowth, Pred = cNPred * (sum(cGrowth)/sum(cNPred)))
    fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
    fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
    
    #Save, gaic, model and age data as part of the output
    gaic <- fit1$aic-fit2$aic
    print(gaic)
    return(gaic)
  })
  
  
  plot(res)
  lines(res)
  
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
oT <- STClu(oT)

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
  dists <- oT$nSum$mDist
  cutoffs <- seq(min(dists),max(dists),(max(dists)-min(dists))/50)
  
  res <- lapply(cutoffs, function(x){
    subT <- dFilt(oT, x)
    subT <- STClu(subT)
    return(subT)
  })
  
  cnum <- sapply(res, function(c) {length(c)})
  
  csize <- sapply(res, function(c) {
    mean(sapply(c, function(x) {length(x)})) 
  })
  
  gmean <- sapply(res, function(c) {
    mean(sapply(c, function(x){
      if(length(x)>1) {
        length(which(x>length(oT$ID)))
      }
      else {
        0
      }
    }))
  })
  
  gcover <- sapply(res, function(c) {
    sum(sapply(c, function(x){
      if(length(x)>1) {
        length(which(x>length(oT$ID)))
      }
      else {
        0
      }
    }))/nrow(oT$gSum)
  })
  
  par(mfrow=c(2,2))
  plot(cutoffs, cnum, xlab="Threshold", ylab="Number of Clusters")
  lines(cutoffs, cnum)
  plot(cutoffs, csize, xlab="Threshold", ylab="Mean Cluster Size")
  lines(cutoffs, csize)
  plot(cutoffs, gmean, xlab="Threshold", ylab="Mean Cluster Growth")
  lines(cutoffs, gmean)
  plot(cutoffs, gcover, xlab="Threshold", ylab="Total Growth Coverage")
  lines(cutoffs, gcover)
}

