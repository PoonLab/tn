require("ape") # Possibly Unneeded
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

#Build a temporary tree from old sequences. Currently uses RaxML for compatibility with pplacer. 
tempTree <- function(sFile="~/Data/Seattle/SeattleB_PRO.fasta", oSFile, oTFile) {
  
  #Obtain Times from sequence file in order to specify the old partition
  times <- getT(sFile)
  keepT <- head(as.numeric(names(table(times))), -1)
  
  tFilt(sFile, keepT, oSFile)
  
  #Create trees, outputted to the temporary tree directory
  RaxMLCall(oSFile, oTFile)
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

#Obtains some likelihood data in order to weight cases based on their recency 
likData <- function(iT){
  
  #Obtain the minimum retrospective distance for each member sequence 
  ##TO-DO: Optimize, this should only happen once. It may belong in the tree-creation call ##
  f <- bind_rows(sapply(1:length(iT$v$ID), function(i){
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
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case age
  #Use this to weight cases by age
  mod <- glm(cbind(Positive, vTotal) ~ tDiff+oeDens, data=ageD, family='binomial')
  weight <- predict(mod, type='response',
                    data.frame(tDiff=max(iT$Time)+1-iT$Time, 
                               oeDens=as.numeric(eTab[as.character(iT$Time)]/vTab[as.character(iT$Time)])))
  
  
  return(weight)
  
}

###^^^^^CUT?^^^^^###

#Import Tree Data and annotate with sequence ID and Time
#Sequences must be dated with the date separated from the id by '_'. 
##TO-DO: Currently only accepts year dates. Work to allow more specific dates. 
impTree <-function(tFile){
  #@param iFile: The name/path of the input file (expecting a newick file)
  #@preturn: An ape tree object with associated lists of sequence ID and Time
  
  #Creating an ape phylogeny object from the newick file, store in a greater list "t" as "p" for phylogeny
  t <- list()
  t$p <- read.tree(tFile)
  
  #Obtain lists of sequence ID and Time
  temp <- sapply(t$p$tip.label, function(x) strsplit(x, '_')[[1]])
  ids <- unname(temp[1,])
  times <- as.numeric(temp[2,])
  t$v <- data.frame(ID=ids, Time=times, stringsAsFactors = F)
  tips <- 1:nrow(t$v)
  
  #Summarize internal branch length information
  t$e <- list()
  t$e$dist <- dist.nodes(t$p)[tips, tips]
  t$e$tDiff <- sapply(tips, function(i){
    t1 <- t$v$Time[i]
    sapply(tips, function(j){ t1 - t$v$Time[j] })
  })
  
  #Summarize information by node, representing the potential to cluster by subtree
  #Obtain the mean branch length under each node
  ##TO-DO: Add bootstrap here. Requires bootstrap to be added in tree creation.
  nodes <- (max(tips)+2):(max(tips)*2-2) 
  des <- lapply(nodes, function(x){Descendants(t$p,x,"all")})
  meanDist <- sapply(des, function(x) {
    x <- x[which(x%in%tips)] 
    mean <- mean(sapply(x, function(tip){t$e$dist[tip,tips[-tip]]}))
    return(mean)
  })
  
  #The minimum retrospective edge to 
  t$f <- bind_rows(lapply(which(t$v$Time>min(t$v$Time)), function(i){
    inc <- t$e$dist[i,-i]
    iTD <- t$e$tDiff[i,-i]
    ret <- inc[which(iTD<0)]
    
    minE <- which(ret==min(ret)[[1]]) 
    minD <- ret[minE]
    minTD <- (iTD[which(iTD<0)])[minE]
      
    df <- data.frame(Index=minE, Distance=minD, tDiff=minTD)
    return(df)
  }))
  
  #Set up a node-summary data frame
  t$n <- list()
  t$n$des <- des
  t$n$agg <- data.frame(ID=nodes, mDist=meanDist)
  names(t$n$des) <- t$n$agg$ID
  
  return(t)
}


#After simulating the growth of trees by placing recent sequences as tips on a fixed ML tree
growthSim <- function(iT, gFile) {
  
  #Obtain a set of trees with tips added
  ts <- read.tree(gFile)
  
  #The nodes in iT associated with new cases in t
  iT$g<- bind_rows(lapply(ts, function(t) {
    
    temp <- sapply(t$tip.label, function(x) strsplit(x, '_')[[1]])
    times <- as.numeric(temp[2,])
    ids <- temp[1,]
    nTip <- which(times>max(iT$v$Time))
    nAdj <- Siblings(t, nTip)
    
    dist <- t$edge.length[which(t$edge[,2]==nTip)]
    
    #If the sibling is a tip, returns the sibling's parent in iT
    if(nAdj<=length(t$tip.label)) {
      adj <- which(iT$v$ID%in%ids[nAdj])
      parent <- iT$p$edge[which(iT$p$edge[,2]==adj),1]

      return(data.frame(ID=ids[nTip], node=parent, dist=dist)) 
    }
    else{
      return(data.frame(ID=ids[nTip], node=nAdj, dist=dist))
    }
  }))
  
  #Calculate the "Breaking distance (ie. the distance required to p)
  iT$g$breakD <- sapply(node, function(x){
    
  })
  
  return(iT)
}

#Cluster using a modified subtree-based method
STClu <- function(iT, maxD) {

  #Identify potential clusters by some criterion
  cluNames <- as.character(subset(iT$n$agg, mDist<maxD)$ID)
  
  #Check the clustering conditions for all possible subtrees and ensure no cherry clusters are added
  subCon <- unname(sapply(cluNames, function(name){
    sum(sapply(setdiff(cluNames, name), function(L) {
      all(iT$n$des[[name]] %in% iT$n$des[[L]])
    }))>0
  }))
  
  #Remove non-clusters from the list
  cluNames <- cluNames[!subCon]

  #Remove internal nodes for summary. Only cases are of interest
  clu <- lapply(cluNames, function(i){
    
    #Obtain Tips already in Cluster
    iDes <- iT$n$des[[i]]
    iTipDes <- iT$v$ID[iDes[which(iDes<length(iT$v$ID))]]
    
    #Add new cases if they're associated with the internal
    gTips <- subset(iT$g, iT$g$node%in%c(iDes,as.numeric(i)))$ID
    return(c(iTipDes, gTips))
  })
  
  iT$c <- list()
  
  #Add singletons and summerize clusters
  sing <- as.list(setdiff(c(iT$v$ID, iT$g$ID), unlist(clu)))
  iT$c$membership <- c(clu,sing)
  temp <- sapply(iT$c$membership, function(x){
    newCases <- length(which(x%in%iT$g$ID))
    oldCases <- length(x) - newCases
    
    return(c(oldCases,newCases))
  }) 
  
  #For growth, record the number of old and new cases in the clusters which contain old cases
  colnames(temp) <- 1:length(temp[1,])
  iT$c$growth <- temp[,which(temp[1,]>0)]
    
  return(iT)
}

#Obtain GAIC at several different cutoffs
GAICRun <- function(iT) {
  
  #Clustering Test
  dists <- iT$nSum$mDist
  cutoffs <- seq(0,0.15,0.15/50)
  
  res <- sapply(cutoffs, function(x){
    print(x)
    subT <- dFilt(iT, x)
    subT <- STClu(subT)
    subT$Weight <- likData(subT, x)
    
    
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
  
  return(res)
  
}

###############Testing

tFile <- "~/Data/Seattle/SeattleB_RAxML/RAxML_result.Tree" 
gFile <- "~/Data/Seattle/SeattleB_pplacer/st.tre"
maxD <- 0.16

oT <- impTree(tFile) 
oT <- growthSim(oT, gFile)
oT <- STClu(oT, maxD)

#Plot Test
plotT <- F
if(plotT){
  
  #Clustering Test
  dists <- oT$nSum$mDist
  cutoffs <-seq(min(dists),0.20,(0.20-min(dists))/50)
  
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

