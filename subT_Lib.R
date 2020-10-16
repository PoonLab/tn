require("ape")
require("phangorn")
require("dplyr")
require("parallel")

#Import Tree Data and output an annotated tree with additional information to assist clustering
impTree <-function(tFile, reVars='/|\\|', varInd=c(5,6,2), dateFormat="%Y-%m-%d", varMan=NA, nCore=detectCores()){
  #@param iFile: The name/path of the input file (expecting a newick file)
  #@param iFile: The name/path of the input file (expecting tn93 output csv)
  #@param reVars: The regular expression used to extract variables from column headers. This is passed to strsplit, creating a vertex of values from the column header
  #@param varInd: A vector of numbers describing the order of variables in the split string. This should describe the index of the unique ID, the Timepoint and the location.
  #               ex. If the header would be split such that the 4th index is the Unique ID, then 4 should be the first number in this list
  #               ID and timepoint are currently required. If the location information is not available, it should be set as "0".
  #@param varMan: Variables can be assigned manually with a csv containing columns of ID, Time point, and Location, in that order. Again, location is not mandatory. 
  #               If this option is used, reVars and varInd, need not be provided. --CURRENTLY UNNUSED--
  #@return: An ape phylo object annotated with the additional data summarized below
  #    $v: A data frame storing vertex information ($ID, $Time, and, if given $Location)
  #    $e: A list of 2 matrices. $dist representing phenetic distance and $tDiff representing time difference.
  #        Additionally, if location information is provided, shared location is a variable stored in $lMatch. 
  #    $n: A list of each node's descendants ($Des), as well as information used to obtain clusters from nodes 
  #        This additional data is stored in ($Info) as $xDist for the longest distance $Bootstrap for the support value.
  #        also, $mTime for the mean date at which sequences were collected and $mtDiff
  
  
  #Obtaining and midpioint rooting an ape phylogeny object from the tree file, store in a greater list "t"
  t <- midpoint(read.tree(tFile))
  
  #Obtain lists of sequence ID and Time
  #Reformat edge list as data table object with predictors extracted from sequence header
  temp <- sapply(t$tip.label, function(x) strsplit(x, reVars)[[1]])
  t$v <- data.table(ID=temp[varInd[[1]],],  
                    Time=as.Date(temp[varInd[[2]],], format=dateFormat), 
                    stringsAsFactors = F)
  t$v <- t$v[order(t$v$Time)]
  
  #Summarize internal branch length information using a list of 2 matrices
  t$e <- list()
  t$e$Dist <- cophenetic.phylo(t)
  rownames(t$e$Dist) <- t$v$ID
  colnames(t$e$Dist) <- t$v$ID
  
  t$e$tDiff <- sapply(t$v$Time, function(x) {t$v$Time-x})
  rownames(t$e$tDiff) <- t$v$ID
  colnames(t$e$tDiff) <- t$v$ID
  
  #In the event that location information is also available
  if(length(varInd)>2) {
    t$v$Location <- temp[varInd[[3]],]
    
    t$e$lMatch <- sapply(t$v$Location, function(x) {
      matches <- rep(F, length(t$v$Location))
      matches[grep(x, t$v$Location)] <- T
      return(matches)
    })
    rownames(t$e$tDiff) <- t$v$ID
    colnames(t$e$tDiff) <- t$v$ID
  }
  
  #Obtain the tip and node names (as numbers)
  #Obtain the list of descendants
  tips <- 1:nrow(t$v)
  nodes <- (max(tips)+1):(max(tips)*2-1) 
  
  #Obtain the full set of descendants at each node
  t$n <- list()
  t$n$Des <- lapply(nodes, function(x){Descendants(t,x,"all")})
  names(t$n$Des) <- as.character(nodes)
  
  #Obtain information which can be used to build clusters and organize the data
  t$n$Info <- as.data.table(bind_rows(lapply(t$n$Des, function(x){
    des <- x[which(x%in%tips)]
    mTime <- mean(t$v$Time[des])
    mtDiff <- mean(abs(t$e$tDiff[des,des]))
    xDist <- max(t$e$Dist[des,des])
    data.frame(mTime=mTime, totTips=length(des), xDist=xDist)
  })))
  
  #Additional information regarding the label of each node, as well as the bootstrap support
  ##- TO-DO: Bootstrap support values are noted differently depending on software used. -## 
  t$n$Info$ID <- as.character(nodes)
  t$n$Info$BootStrap <- as.numeric(t$node.label)
  
  #A comprimise for the root
  t$n$Info$BootStrap[is.na(t$n$Info$BootStrap)]
  
  return(t)
}

#After simulating the growth of trees by placing recent sequences as tips on a fixed ML tree
growthSim <- function(iT, gFile) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param gFile: The growth file from a pplacer run for all new cases
  #@return: The input tree annotated with growth information stored as $g.
  #    $nID: The ID if the new node. This will be a full tip label
  #    $xDist: The distance between this node and it's the most distant tip in the tree that is formed
  #    $oNode: The index corresponding to the neighbour of the newly added tip.
  #            If this is "NA", the new tip was added to a terminal branch.
  
  #Obtain a set of trees with new tips added
  #This is one tree for each new case
  ts <- read.tree(gFile)
  df <- bind_rows(lapply(ts, function(t){
    
    #Extract the subtree which contains the 
    nID <- setdiff(t$tip.label, iT$tip.label)
    nTip <- which(t$tip.label%in%nID)
    p <- t$edge[which(t$edge[,2]%in%nTip), 1]
    t <- extract.clade(t, p)
    
    #Obtain the node corresponding to old nodes within this subtree
    #This index is named based on that nodes number in the old tree
    oIDs <- setdiff(t$tip.label, nID)
    oTips <- which(iT$tip.label%in%oIDs)
    oNode <- mrca.phylo(iT, oTips)
    
    #To catch a singleton connecting to a new tip
    if(is.null(oNode)){
      oNode <- NA
    }
  
    data.frame(nID=nID, xDist=max(dist.nodes(t)), oNode=oNode)
  }))  
  
  #Add this information as a data table under g.
  iT$g <- as.data.table(df)
  
  return(iT)
}

#Cluster using a modified subtree-based method
#Clusters are defined as subtrees with a mean tip-to tip distance under some maximum
STClu <- function(iT, maxD, minB=0, meanD=F) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param maxD: The maximum distance criterion defining clusters
  #@return: The tree annotated with cluster size, growth and membership
  #         Growth is summarized as a matrix of old cases, new cases and predictor values
  
  #Identify potential clusters by some criterion as well as tips within those potential clusters
  if(meanD){
    cluNames <- as.character(subset(iT$n$agg, (mDist<=maxD)&(BootStrap>=minB))$ID)
  }else {
    cluNames <- as.character(subset(iT$n$agg, (xDist<=maxD)&(BootStrap>=minB))$ID)
  }

  
  cluTips <- unlist(iT$n$des[cluNames]) 
  cluTips <- unique(as.numeric(cluTips[cluTips<length(iT$v$ID)]))
  
  #A loop to remove any clusters which are subclusters of a larger cluster
  temp <- vector()
  
  #Cycle through all clustered and assure that no cluster is simply part of a larger cluster
  #Only the clusters which are not parts of a larger cluster are reasssigned to temp
  while(length(cluTips)>0) {
    tip <- cluTips[1]
    tClu <- cluNames[sapply(iT$n$des[cluNames], function(iClu) {tip%in%iClu})]
    tCluS <- as.numeric(sapply(tClu, function(cluName) {length(iT$n$des[[cluName]])}))
    tClu <- tClu[which(tCluS==max(tCluS))[[1]]]
    temp <- c(temp,tClu)
    cluTips <- cluTips[-which(cluTips%in%iT$n$des[[tClu]])]
  }
  
  #Track clusters in the trees 'n' list item
  iT$n$agg$Clustered <- rep(FALSE, nrow(iT$n$agg))
  iT$n$agg$Clustered[which(iT$n$agg$ID%in%cluNames)] <- TRUE
  cluNames <- temp
  
  iT$c$cluNames <- cluNames

  #Obtain clusters in the form of membership lists, including new cases
  clu <- lapply(cluNames, function(i){
    
    #Obtain Tips, ids and edges already in Cluster 
    iDes <- iT$n$des[[i]]
    iDesTs <- iDes[which(iDes<length(iT$v$ID))]
    iDesEs <- which(iT$p$edge[,2]%in%iDes)
    
    #New cases associated with an internal node within the potential cluster
    gTips <- subset(iT$g, (oE%in%iDesEs)) #&(Bootstrap>=minB))
    
    #To catch the event that no new tips are added
    if(nrow(gTips)>0) {
      #To Find any new tips that need to be added to membership
      #These tips will not be added to a cluster if they would increase the mean distance above dMax
      ##TO-DO: Test function, this currently catches no breaks in real data (they could be rare)
      gTips$breakCon <- sapply(1:nrow(gTips), function(j) {
        gTip <- gTips[j,]
        mDist <- subset(iT$n$agg, ID==as.numeric(i))$mDist
        adj <- iT$p$edge[gTip$oE,2]
        
        #Recalculate mean edge distance
        oDists <- rep(mDist, choose(length(iDesTs),2))
        nDists <- iT$e$dist[adj, setdiff(iDesTs, adj)]-gTip$pendLength+gTip$distLength
        nmDist <- mean(c(oDists,nDists, gTip$pendLength+gTip$distLength))
 
        return(nmDist>=maxD)
      })
    
      
      if(meanD) {
        gTips <- subset(gTips, !breakCon) 
      }else {
        gTips <- subset(gTips, (distLength+pendLength)<=maxD)
      }
    }
    
    #To catch the event that no new tips are added
    if(nrow(gTips)==0) {
      gTips <- NULL
    }
    
    return(c(iT$v$ID[iDesTs], gTips$nID))
  })
  
  #Obtain old and new singletons.
  oSing <- as.list(setdiff(iT$v$ID, unlist(clu)))
  nSing <- as.list(setdiff(iT$g$nID, unlist(clu)))
  
  #If an old singleton would form a cluster with an unclustered new singleton we capture that too
  clu <- c(clu, lapply(oSing, function(i){
    
    #New cases associated with an internal node within the potential cluster
    iE <- which(iT$p$edge[,2]==which(iT$v$ID%in%i))
    gTips <- subset(iT$g, (iT$g$oE==iE)) #&(Bootstrap>=minB))
    
    if(nrow(gTips)>0){
      gTips <- subset(gTips, (distLength+pendLength)<=maxD)
    }
    return(c(i, gTips$nID))
  }))
  
  #Add old and new singletons to updated clusters
  oSing <- as.list(setdiff(iT$v$ID, unlist(clu)))
  nSing <- as.list(setdiff(iT$g$nID, unlist(clu)))
  
  iT$c$membership <- c(clu,oSing,nSing)
  
  #Return a Data frame of important cluster growth information
  temp <- bind_rows(lapply(1:length(iT$c$membership), function(i){
    
    x <- iT$c$membership[[i]]
    oldCases <- length(which(x%in%iT$v$ID))
    
    #For growth, record the number of old and new cases in the clusters which contain old cases
    if(oldCases==0) {
      return(data.frame(Old=NULL,New=NULL,mRec=NULL,cSize=NULL, mDist=NULL))
    }
    else {
      newCases <- length(x) - oldCases
      mRec <- max(iT$v$Time)+1-mean(subset(iT$v, ID%in%x)$Time)
      
      #Calculate the total length in branches that is taken up by a given subtree
      if(i<=length(cluNames)){
        totLength <- subset(iT$n$agg, ID%in%cluNames[i])$mDist
        mDist <- subset(iT$n$agg, ID%in%cluNames[i])$mDist
      }
      
      #For singletons
      else {
        totLength <- iT$p$edge.length[which(iT$p$edge[,2]%in%which(iT$v$ID%in%x))]
        mDist <- 0
      }
      
      return(data.frame(Old=oldCases,New=newCases,mRec=mRec,totLength=totLength, mDist=mDist))
    }
  }))
  
  iT$c$growth <- temp
  return(iT)
}

#Obtain GAIC at several different cutoffs
GAICRun <- function(iT, maxDs, minBs=0, rand=F) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param cutoffs: A set of maximum distance criterion defining clusters
  #@return: A list of analysis based on a fitting actual new data to predicted cluster growth
  #         This most importantly includes GAIC
  
  #Obtain the subtree and a list of analysis 
  runRes <- lapply(maxDs, function(maxD){
    x <- lapply(minBs, function(minB){
      
      #Create subTree
      subT <- STClu(iT, maxD, minB)
      
      #Check positives now that we have params
      subT$f$Positive <- (subT$f$xDist<=maxD) #&(subT$f$BootStrap>=minB)
      
      #Layout of Time and time lag information
      tTab <- table(iT$v$Time)
      tDiffs <- sort(unique(abs(as.vector(iT$e$tDiff))))
      tdTab <- rep(nrow(iT$v)-1, length(tDiffs))
      names(tdTab) <- tDiffs
      largeTD <- tDiffs[which(tDiffs>(ceiling(max(tDiffs/2))))]
      leftOut <- sapply(1:length(largeTD), function(i){
        td <- largeTD[[i]]
        ts <- as.numeric(names(tTab))[between(as.numeric(names(tTab)), (max(iT$v$Time)-td+1), (min(iT$v$Time)+td-1))]
        sum(tTab[as.character(unique(ts))]) 
      })
      
      #Obtains the "attempts". Or how many tips find a given tDiff possible
      #For example, it may be impossible for central time points to see the largest time difference in the set
      tdTab[as.character(largeTD)] <- tdTab[as.character(largeTD)]-leftOut
      
      #Sort Data into Age Data
      names(tdTab) <- tDiffs
      posTab <- rep(0,length(tDiffs))
      names(posTab) <- tDiffs
      tempTab <- table(round(subset(subT$f, Positive)$tDiff))
      posTab[names(tempTab)] <- as.numeric(tempTab)
      ageD <- data.frame(tDiff=tDiffs, Positives=as.numeric(posTab), 
                         Total=as.numeric(tdTab))
      
      #Weighting whole clusters based on their size and mean recency
      mod <- glm(cbind(Positives,Total)~tDiff, data=ageD, family='binomial')
      
      #Individual node weighting
      if(rand){ 
        #Random Weight Tests
        subT$v$Weight <- sample(1:20, nrow(subT$v), replace=T)
      } else {
        #Tips are weighted based on recency
        subT$v$Weight <- predict(mod, type='response', data.frame(tDiff=max(subT$v$Time)-subT$v$Time+1))
      }
        
      #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
      df1 <- data.frame(Growth = subT$c$growth$New, Pred = sapply(subT$c$membership[which(subT$c$growth$Old>0)], function(x){
        members <- subset(subT$v, ID%in%x)
        sum(members$Weight)
      })) 
      df2 <- data.frame(Growth = subT$c$growth$New, Pred = sapply(subT$c$membership[which(subT$c$growth$Old>0)], function(x){
        members <- subset(subT$v, ID%in%x)
        nrow(members)
      }))
        
      #Cluster growht prediction model
      fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
      fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
      
      
      #Whole Cluster Weighting 
      ##Possibly Not approp.
      if(F){
        fit1 <- glm(New ~ mRec+Old, data = subT$c$growth, family = "poisson")
        fit2 <- glm(New ~ Old, data = subT$c$growth, family = "poisson")
      }

      #Save, gaic, model and age data as part of the output
      res <- list()
      res$gaic <- fit1$aic-fit2$aic
      res$propFit <- fit1$aic
      #res$nullFit <- fit2
      #res$mod <- mod
      res$par <- c(maxD, minB)
      res$ageD <- ageD
      res$growth <- subT$c$growth
      
      #print(res$par)
      #print(res$gaic)
      #print(length(subT$c$cluNames))
      print(res$par)
      print(res$gaic)
      
      return(res)
    })
  })
  
  return(runRes)
}

if(F){
  tFile <- "~/Data/Seattle/strefpackages/st1.refpkg/sttree1.nwk"
  gFile <- "~/Data/Seattle/strefpackages/growthFiles/growthFile1.tree"
  oT <- impTree(tFile, reVars='_', varInd = c(1,2), dateFormat = "%Y") 
  oT <- growthSim(oT, gFile)
  
  iT <- oT
  minB <- 1.0
  maxD <- 0.02
  meanD <- F
  rand <- F
  
  cutoffs <- seq(0,0.2,0.002)
  #cutB <- seq(1,0,-0.02)
  res <- GAICRun(oT, cutoffs)
  
  #saveRDS(res, "resST")
  
  gaics <- sapply(res, function(x){x[[1]]$gaic})
  plot(cutoffs, gaics, xlab = "Mean Cutoff", ylab="GAIC")
  lines(cutoffs, gaics, col="red")
  
  
  blahs <- sapply(res, function(x){nrow(x[[1]]$growth)})
  plot(cutoffs, blahs, xlab = "Mean Cutoff", ylab="GAIC")
  lines(cutoffs, blahs)

  
  
  blahs <- sapply(res, function(x){(x[[1]]$nullFit$aic)})
  plot(cutoffs, blahs, xlab = "Mean Cutoff", ylab="GAIC")
  lines(cutoffs, blahs)
  blahs <- sapply(res, function(x){(x[[1]]$propFit$aic)})
  lines(cutoffs, blahs)
  
  
  #res <- GAICRun(oT, cutoffs,rand=F)
  #modAIC <- sapply(res, function(x){x[[1]]$propFit$aic})
  #nullAIC <- modAIC-gaics
}