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
  tips <- 1:length(ids)
  
  tEs <- which(t$p$edge[,2]%in%tips)
  t$v <- data.frame(ID=ids, Time=times, stringsAsFactors = F,
                    termDist=t$p$edge.length[tEs])
  #Summarize internal branch length information 
  #Time differences, phenetic distance matrix, pairwise distances, and a table of time differences
  t$e <- list()
  t$e$dist <- dist.nodes(t$p)
  t$e$tDiff <- sapply(tips, function(i){
    t1 <- t$v$Time[i]
    sapply(tips, function(j){ t1 - t$v$Time[j] })
  })
  
  el <- t$e$dist[1:nrow(t$v),1:nrow(t$v)]
  t$el <- as.numeric(unlist(lapply(1:nrow(t$v), function(x){el[x,x:nrow(t$v)]})))
  
  #Unnused
  #t$e$tdTab <- table(abs(t$e$tDiff))

  return(t)
}

#Summarize information by node, this is used to cluster by subtree and add this to the tree file
#Obtain the mean branch length under each node
##TO-DO: Reassess any of this that could be moved to initial tree creation. Or sub-sampling ##
##TO-DO: Simplify. Needs to be re-run for robust tests. ##
nodeInfo <- function(iT) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@return: The tree annotated with 2 Objects. The list of descendant tips under each given node  
  #         and a dataframe which summarizes information relevant to clustering
  
  #Obtain the tip and node names (as numbers)
  #Obtain the list of descendants
  tips <- 1:nrow(iT$v)
  nodes <- (max(tips)+1):(max(tips)*2-1) 
  
  iT$n <- list()
  iT$n$des <- lapply(nodes, function(x){Descendants(iT$p,x,"all")})
  names(iT$n$des) <- as.character(nodes)
  
  #Obtain information for each node
  #This includes information around the subtree defined by each node
  iT$n$agg <- bind_rows(lapply(iT$n$des, function(x){
    tips <- x[which(x<=length(iT$v$ID))]
    
    mTime <- mean(iT$v$Time[tips])
    mRec <- max(iT$v$Time)+1 - mTime
    mDist <- mean(unlist(lapply(tips, function(tip){iT$e$dist[tip,tips[tips>tip]]})))
    xDist <- max(unlist(lapply(tips, function(tip){iT$e$dist[tip,tips[tips>tip]]})))
    mtDiff <- mean(unlist(lapply(tips, function(tip){abs(iT$e$tDiff[tip,tips[tips>tip]])})))
    totEdge <- sum(iT$p$edge.length[which(iT$p$edge[,2]%in%x)])
    
    data.frame(mDist=mDist, mtDiff=mtDiff, totEdge=totEdge, mRec=mRec, mTime=mTime, totTips=length(tips), xDist=xDist)
  }))
  iT$n$agg$ID <- as.character(nodes)
  iT$n$agg$BootStrap <- as.numeric(iT$p$node.label)/100
  
  ##TO-DO: Confirm, this is statistaically valid.
  iT$n$agg$BootStrap[is.na(iT$n$agg$BootStrap)] <- 1 # A comprimise for the root.
  
  #Based on Cherry nodes, obtain some frequency info
  if(T) {
    #Obtain Cherry nodes for the sake of frequency function analysis
    tEs <- which(iT$p$edge[,2]%in%tips)
    
    #Analyze things about each cherry node
    iT$f <- bind_rows(lapply(tEs, function(tE){
      
      cN <- (iT$p$edge[,1])[tE]
      cT <- (iT$p$edge[,2])[tE]
      cA <- (iT$p$edge[,2])[setdiff(which(iT$p$edge[,1]==cN), tE)]
      
      if(cN!=(nrow(iT$v)+1)){
        #The length of the edge of interest (without the node that cT makes)
        eL <- (iT$p$edge.length[which(iT$p$edge[,2]==cA)] +
                 iT$p$edge.length[which(iT$p$edge[,2]==cN)])
      }
      
      #If adjacent 'cA' is a node
      if(cA>length(tips)){
        nNmRec <- iT$n$agg$mTime[which(as.numeric(iT$n$agg$ID)==cA)]
        tDiff <- abs(nNmRec-iT$v$Time[cT])
      } else {
        tDiff <- abs(iT$e$tDiff[cT,cA])
      }
      
      
      infoR <- subset(iT$n$agg, ID%in%as.character(cN))
      
      if(cN==(nrow(iT$v)+1)) {
        return(data.frame(Node=NULL, Tip=NULL, tDiff=NULL, BootStrap=NULL, mDist=NULL, xDist=NULL,  branchL=NULL))
      } else{
        return(data.frame(Node=cN, Tip=cT, tDiff=tDiff, BootStrap=infoR$BootStrap, mDist=infoR$mDist, xDist=infoR$xDist, branchL=eL))          
      }
            
    }))
  }
  
  #Alternate - Placement based analysis
  if(F) {
    #Frequencey information model, placement weighting
    iT$f <- bind_rows(lapply(tips, function(tip){
      print(tip)
      
      pth <- as.character(nodepath(iT$p, tip, nodes[[1]])[-1])
      invpth <- as.character(setdiff(nodes, pth))
      
      mem <- subset(iT$n$agg, ID%in%pth) %>% select(mDist, mTime, BootStrap, totEdge, totTips, ID)
      invmem <- subset(iT$n$agg, ID%in%invpth) %>% select(mDist, mTime, BootStrap, totEdge, totTips, ID)
      
      mem$tDiff <- abs((mem$mTime*mem$totTips-iT$v$Time[tip])/(mem$totTips-1)-iT$v$Time[tip])
      invmem$tDiff <- abs(invmem$mTime-iT$v$Time[tip])
      
      mem$mem <- 1
      invmem$mem <- 0
      
      return(rbind(mem,invmem))
    }))
  }

  return(iT)
}

#After simulating the growth of trees by placing recent sequences as tips on a fixed ML tree
growthSim <- function(iT, gFile) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param gFile: The growth file from a pplacer run for all new cases
  #@return: The tree annotated with growth information (ie. their closest neighbour and the distance)
  
  #Obtain a set of trees with new tips added
  #This is one tree for each new case
  ts <- read.tree(gFile)
  
  #Summarize a table of information for new tips placed on the tree
  df <- bind_rows(lapply(ts, function(t){
    
    #Obtain potential placements
    nIDs <- setdiff(t$tip.label, iT$p$tip.label)
    bss <- sapply(nIDs, function(s){
      sl <-  strsplit(s, '_')[[1]]
      s <- sl[[length(sl)]]
      bs <- as.numeric(substr(s, 3, nchar(s)))
    })
    nIDi <- which.max(bss)
    
    #Find most certain placement by bootstrap value
    bs <- bss[[nIDi]]
    nID <- nIDs[[nIDi]]
    nT <- which(t$tip.label%in%nID)
    sT <- unlist(Siblings(t, nT))
    
    #Find the edges of interest on the old Tree 
    oE <- sapply(sT, function(x){
      
      #If our new tip is adjacent to another tip, use that tips tip label
      if(x<=length(t$tip.label)){
        sID <- t$tip.label[x]
        oT <- which(iT$p$tip.label%in%sID)
        return(which(iT$p$edge[,2]%in%oT))  
      }
      
      #If our new tip is adjacent to an internal node, 
      #We must consider the node's position in a node list list without new nodes
      else{
        nN <- which(t$node.label%in%"")+length(t$tip.label)
        noN <- setdiff((length(t$tip.label)+1):(length(t$node.label)+length(t$tip.label)), nN)
        oN <- which(noN%in%x)+length(t$tip.label)
        return(which(iT$p$edge[,2]%in%oN))
      }
    })
    
    #Obtain edge lengths of the new tip and it's neighbour leading to the created parent node
    distLength <- t$edge.length[which(t$edge[,2]%in%sT)]
    pendLength <- t$edge.length[which(t$edge[,2]%in%nT)]
    nID <- (strsplit(nIDs[[nIDi]], '_')[[1]])[[1]]
    
    return(data.frame(nID=nID, oE=oE, Bootstrap=bs, distLength=distLength, pendLength=pendLength))
  }))
  
  iT$g <- df
  
  return(iT)
}

#Cluster using a modified subtree-based method
#Clusters are defined as subtrees with a mean tip-to tip distance under some maximum
##TO-DO: Simplify. Currently the Root of Speed issues. ##
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
      res$propFit <- fit1
      res$nullFit <- fit2
      res$mod <- mod
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

#Subsample a tree with given information 
subSample <- function(iT, ssize=1200) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param ssize: How many cases to take 
  #@return: The subsample of sample size "ssize"
  
  i <- sample(1:length(iT$v$ID), ssize)
  
  #Insure that the tips sampled are represented in v,p and e.
  iT$v <- iT$v[i,]
  iT$e$f <- iT$e$f[i,]
  iT$p <- keep.tip(iT$p, i)
  #iT$e$dist <- dist.nodes(iT$p)
  #iT$e$tDiff <- iT$e$tDiff[i,i]
  
  return(iT)
}

if(T){
  tFile <- "~/Data/Seattle/IqTree_Bootstrap/st.refpkg/SeattleB_PRO_Filt.fasta.treefile"
  gFile <- "~/Data/Seattle/IqTree_Bootstrap/st.tre"
  oT <- impTree(tFile) 
  oT <- nodeInfo(oT)
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