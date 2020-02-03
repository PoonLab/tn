require("ape") # Possibly Unneeded
require("phangorn")
library("dplyr",verbose = FALSE)

#Import Tree Data and annotate with sequence ID and Time
#Sequences must be dated with the date separated from the id by '_'. 
##TO-DO: Currently only accepts year dates. Work to allow more specific dates. 
impTree <-function(tFile){
  #@param iFile: The name/path of the input file (expecting a newick file)
  #@return: A list of 3 Objects. The ape phylo object, a vertex list, and a list of edge information
  
  #Obtaining and midpioint rooting an ape phylogeny object from the newick file, store in a greater list "t" as "p" for phylogeny
  t <- list()
  t$p <- midpoint(read.tree(tFile))
  
  #Obtain lists of sequence ID and Time
  temp <- sapply(t$p$tip.label, function(x) strsplit(x, '_')[[1]])
  ids <- unname(temp[1,])
  times <- as.numeric(temp[2,])
  t$v <- data.frame(ID=ids, Time=times, stringsAsFactors = F)
  tips <- 1:nrow(t$v)
  
  #Summarize internal branch length information 
  #Time differences, phenetic distance matrix, pairwise distances, and a table of timne differences
  t$e <- list()
  t$e$dist <- dist.nodes(t$p)
  t$e$tDiff <- sapply(tips, function(i){
    t1 <- t$v$Time[i]
    sapply(tips, function(j){ t1 - t$v$Time[j] })
  })
  t$e$tdTab <- table(t$e$el$tDiff)
  
  #Create a pairwise edgelist similair to tn93 output. 
  ##TO-DO: Currently acts as an alternative data structure for the storedge of edge info. Possibly redundant
  t$e$el <- bind_rows(lapply(tips, function(i){
    inc <- tips[which(tips>i)]
    tip <- rep(i, length(inc))
    dist <- t$e$dist[i, inc]
    tDiff <- t$e$tDiff[i, inc]
    return(data.frame(ID1=tip, ID2=inc, Distance=dist, tDiff=abs(tDiff)))
  }))
  
  #The minimum retrospective edge to each tip (except those at the earliest time point)
  t$v <- cbind(t$v, bind_rows(lapply(tips, function(i){
    
    if(t$v$Time[i]==min(t$v$Time)){
      return(data.frame(ClosestV=NA, Distance=NA, tDiff=NA))
    }
    else {
      inc <- t$e$dist[i,-i]
      iTD <- t$e$tDiff[i,-i]
      ret <- inc[which(iTD<0)]
      
      minE <- which(ret==min(ret))[[1]] 
      minD <- ret[minE]
      minTD <- (iTD[which(iTD<0)])[minE]
      
      df <- data.frame(ClosestV=minE, Distance=minD, tDiff=minTD)
      return(df)
    }

  })))
  
  return(t)
}

#Summarize information by node, representing the potential to cluster by subtree
#Obtain the mean branch length under each node
##TO-DO: Add bootstrap here. Requires bootstrap to be added in tree creation.
nodeInfo <- function(iT) {
  
  tips <- 1:nrow(iT$v)
  nodes <- (max(tips)+2):(max(tips)*2-2) 
  des <- lapply(nodes, function(x){Descendants(iT$p,x,"all")})
  
  #Obtain the average within subtree patristic distances
  meanDist <- sapply(des, function(x) {
    x <- sort(x[which(x%in%tips)]) 
    m <- mean(unlist(sapply(head(x,-1), function(tip){iT$e$dist[tip,x[which(x>tip)]]})))
    return(m)
  })
  
  #Obtain the average recency of subtrees
  meanRecency <- sapply(des, function(x) {
    desTips <- x[which(x%in%tips)]
    max(iT$v$Time)+1 - mean(iT$v[desTips,]$Time)
  })
  
  #Set up a node-summary data frame
  iT$n <- list()
  iT$n$des <- des
  iT$n$agg <- data.frame(ID=nodes, mDist=meanDist, mRec=meanRecency)
  names(iT$n$des) <- iT$n$agg$ID
  
  return(iT)
}

#After simulating the growth of trees by placing recent sequences as tips on a fixed ML tree
growthSim <- function(iT, gFile) {
  
  #Obtain a set of trees with tips added
  ts <- read.tree(gFile)
  
  #The nodes in iT associated with new cases in t
  iT$g<- bind_rows(lapply(ts, function(t) {
    
    #Find the terminal branch length of the new tip and the node it is adjacent to
    temp <- sapply(t$tip.label, function(x) strsplit(x, '_')[[1]])
    times <- as.numeric(temp[2,])
    ids <- temp[1,]
    sharedIds <- which(ids%in%iT$v$ID)
    nTip <- which(times==max(times))
    
    tips <- 1:length(t$tip.label)
    dists <- dist.nodes(t)[tips, sharedIds]
    eRet <- dists[nTip,-nTip]
    minRet <- which(eRet==min(eRet))[[1]]
    nAdj <- ids[[as.numeric(names(eRet)[minRet])]]
    parent <- iT$p$edge[which(iT$v$ID%in%nAdj),1]
    
    data.frame(ID=nAdj, Parent=parent, dist=min(eRet))

  }))
  
  return(iT)
}

#Cluster using a modified subtree-based method
STClu <- function(iT, maxD) {

  #Identify potential clusters by some criterion as well as tips within those potential clusters
  cluNames <- as.character(subset(iT$n$agg, mDist<maxD)$ID)
  cluTips <- unlist(iT$n$des[cluNames]) 
  cluTips <- unique(as.numeric(cluTips[cluTips<length(iT$v$ID)]))
  
  #A loop to remove any clusters which are subclusters of a larger cluster
  temp <- vector()
  
  while(length(cluTips)>0) {
    
    tip <- cluTips[1]
    tClu <- cluNames[sapply(iT$n$des[cluNames], function(iClu) {tip%in%iClu})]
    tCluS <- sapply(tClu, function(cluName) {length(iT$n$des[[cluName]])})
    tClu <- tClu[which(tCluS==max(tCluS))[[1]]]
    temp <- c(temp,tClu)
    cluTips <- cluTips[-which(cluTips%in%iT$n$des[[tClu]])]
    
  }
  
  #Mark down clusters in the trees 'n' list item
  iT$n$agg$Clustered <- rep(FALSE, nrow(iT$n$agg))
  iT$n$agg$Clustered[which(iT$n$agg$ID%in%cluNames)] <- TRUE
  
  cluNames <- temp

  #Remove internal nodes for summary. Only cases are of interest
  clu <- lapply(cluNames, function(i){

    #Obtain Tips already in Cluster and the mean distance of this cluster
    iDes <- iT$n$des[[i]]
    iTipDes <- iT$v$ID[iDes[which(iDes<length(iT$v$ID))]]
    mDist <- subset(iT$n$agg, ID==as.numeric(i))$mDist
    n <- length(iTipDes)
    
    #New cases associated with an internal node within the potential cluster
    gTips <- subset(iT$g, iT$g$node%in%c(iDes,as.numeric(i)))
    
    if(nrow(gTips)>0) {
      #Obtaining only new cases satisfying distance requirement
      breakCon <- sapply(1:nrow(gTips), function(j) {
        gTip <- gTips[j,]
        gDist <- iT$e$dist[gTip$node, iDes[which(iDes<length(iT$v$ID))]] + gTip$dist
        bDist <- (sum(gDist) + mDist*factorial(n)) / factorial(n+1)
        
        return(bDist>maxD)
      })
      
      gTips <- gTips[!breakCon,]
    }
    
    if(nrow(gTips)==0) {
      gTips <- NULL
    }
  
    return(c(iTipDes, gTips$ID))
  })
  
  iT$c <- list()
  
  #Add singletons and summerize clusters
  sing <- as.list(setdiff(c(iT$v$ID, iT$g$ID), unlist(clu)))
  iT$c$membership <- c(clu,sing)
  
  temp <- sapply(1:length(iT$c$membership), function(i){
    
    x <- iT$c$membership[[i]]
    newCases <- length(which(x%in%iT$g$ID))
    oldCases <- length(x) - newCases
    meanRec <- max(iT$v$Time)+1-mean(subset(iT$v, ID%in%x)$Time)
    
    return(c(oldCases,newCases, meanRec))
  }) 
  
  #For growth, record the number of old and new cases in the clusters which contain old cases
  colnames(temp) <- 1:length(temp[1,])
  iT$c$growth <- temp[,which(temp[1,]>0)]
    
  return(iT)
}

#Obtain GAIC at several different cutoffs
GAICRun <- function(iT, cutoffs) {
  
  #Clustering Test
  runRes <- lapply(cutoffs, function(maxD){
    print(maxD)
    subT <- STClu(iT, maxD)
    
    #Obtain successes (retrospective growth) and attempts (possible retrospective growths)
    tTab <- as.numeric(table(subT$v$Time))
    tDiffs <- as.numeric(names(subT$e$tdTab))
    ageD <- bind_rows(lapply(tDiffs[tDiffs>0] , function(x){
      data.frame(tDiff=x, 
                 Positive=nrow(subset(subT$e$el, Distance<maxD & tDiff==x)), 
                 vTotal=nrow(subset(subT$e$el, tDiff==x)))
    }))

    #Weighting cases using a model based on case recency
    mod <- glm(cbind(Positive, vTotal) ~ tDiff, data=ageD, family='binomial')
    cluWeights <- predict(mod, type='response', data.frame(tDiff=as.numeric(subT$c$growth[3,])))
    
    #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
    df1 <- data.frame(Growth = subT$c$growth[2,], Pred = subT$c$growth[2,]*cluWeights) 
    df2 <- data.frame(Growth = subT$c$growth[2,], Pred = subT$c$growth[2,]*(subT$c$growth[1,]/sum(subT$c$growth[1,])))
    fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
    fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
    
    res <- list()
    
    #Save, gaic, model and age data as part of the output
    res$gaic <- fit1$aic-fit2$aic
    res$propFit <- fit1
    res$nullFit <- fit2
    res$mod <- mod
    res$ageD <- ageD

    return(res)
  })
  
  return(runRes)
  
}

#Subsample a tree with given information 
subSample <- function(iT, ssize=1200) {
  
  i <- sample(1:length(iT$v$ID), ssize)
  
  iT$v <- iT$v[i,]
  iT$p <- keep.tip(iT$p, i)
  iT$e$dist <- dist.nodes(iT$p)
  iT$e$tDiff <- iT$e$tDiff[i,i]
  iT$e$el <- subset(iT$e$el, (ID1%in%i)&(ID2%in%i))
  iT$e$tdTab <- table(iT$e$el$tDiff)
  
  return(iT)
}

###############Testing

#tFile <- "~/Data/Seattle/SeattleB_PRO_iq_pplacer/st.refpkg/SeattleB_PRO_Filt.fasta.tree"
#gFile <- "~/Data/Seattle/SeattleB_PRO_iq_pplacer/st.tre"


tFile <- "~/Data/Tennessee/TennesseeB_Trim_Diag_IqTree_nm/TennesseeB_Trim_Diag_Filt.fasta.treefile" 
gFile <- "~/Data/Tennessee/TennesseeB_Trim_Diag_pplacer/tn.tre"

oT <- impTree(tFile) 

resList <- lapply(1:50, function(x) {
  print(x/50)
  
  sT <- subSample(oT)
  sT <- nodeInfo(sT)
  sT <- growthSim(sT, gFile)
  cutoffs <- seq(0,0.12,0.005)
  res <- GAICRun(sT, cutoffs)
  
  return(res)
})

saveRDS(resList, "~/Data/Tennessee/subTnsubB_RD.rds")

oT <- nodeInfo(oT)
oT <- growthSim(oT, gFile)

gaics <- lapply(resList, function(res){
  gaics <- sapply(res, function(x){x$a$gaic})
})

plot(oT$p, show.tip.label = F)
add.scale.bar(x=0, y=-40)
mean(oT$p$edge.length[which(oT$p$edge[,2]%in%1:nrow(oT$v))])

cutoffs <- seq(0,0.20,0.005)
res <- GAICRun(oT, cutoffs)

gaics <- sapply(res, function(x){x$a$gaic})
plot(cutoffs, gaics, xlab = "Mean Cutoff", ylab="GAIC")
lines(cutoffs, gaics)
modAIC <- sapply(res, function(x){x$a$propFit$aic})
nullAIC <- modAIC-gaics

ageDs <- sapply(res, function(x){x$ageD})


saveRDS(res, "~/Data/Tennessee/STtnsubB_GD_nm.rds")
