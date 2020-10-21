require("ape")
require("phangorn")
require("dplyr")
require("parallel")

#Import Tree Data and output an annotated tree with additional information to assist clustering
impTree <-function(tFile, reVars='/|\\|', varInd=c(5,6,2), dateFormat="%Y-%m-%d", varMan=NA, nCore=detectCores()){
  #####
  #@param iFile: The name/path of the input file (expecting a newick file)
  #@param iFile: The name/path of the input file (expecting tn93 output csv)
  #@param reVars: The regular expression used to extract variables from column headers. This is passed to strsplit, creating a vertex of values from the column header
  #@param varInd: A vector of numbers describing the order of variables in the split string. This should describe the index of the unique ID, the Timepoint and the location.
  #               ex. If the header would be split such that the 4th index is the Unique ID, then 4 should be the first number in this list
  #               ID and timepoint are currently required. If the location information is not available, it should be set as "0".
  #@param varMan: Variables can be assigned manually with a csv containing columns of ID, Time point, and Location, in that order. Again, location is not mandatory. 
  #               If this option is used, reVars and varInd, need not be provided. --CURRENTLY UNNUSED--
  #@param nCore: The number of cores used for multi-threading. --CURRENTLY UNNUSED--
  #@return: An ape phylo object annotated with the additional data summarized below
  #    $v: A data frame storing vertex information ($ID, $Time, and, if given $Location)
  #    $e: A list of 2 matrices. $dist representing phenetic distance and $tDiff representing time difference.
  #        Additionally, if location information is provided, shared location is a variable stored in $lMatch. 
  #    $n: A list of each node's descendants ($Des), as well as information used to obtain clusters from nodes 
  #        This additional data is stored in ($Info) as $xDist for the longest distance $Bootstrap for the support value.
  #        also, $mTime for the mean date at which sequences were collected and $mtDiff
  #    $f: A list of data used to train a predictive model. This focuses on instances of terminal branches joining the tree
  ....
  #####
  
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
  
  #A compromise for the root
  t$n$Info$BootStrap[is.na(t$n$Info$BootStrap)] <- 1

  #Obtain the direct parent node of each tip 
  #This information is stored in $f
  termE <- t$edge[which(t$edge[,2]<length(t$tip.label)),]
  t$f <- data.table(Tip = termE[,2], tNode = termE[,1])
  
  #Obtain the neighbour node to each tip (this could be an internal node or another tip)
  #This also stores the time difference and maximum distance given this neighbour node
  temp <- sapply(1:nrow(t$f), function(i){
    tNode <- t$f[i, (tNode)]
    cTips <- t$edge[which(t$edge[,1]%in%tNode),2]
    neighbour <- setdiff(cTips, t$f[i,(Tip)])
    
    #Obtain the largest distance and the time difference
    #There are 2 cases for time difference calculation
    xDist <- t$n$Info[ID %in% as.character(tNode), (xDist)]
    tDiff <- ifelse(neighbour>length(t$tip.label),
                    t$v$Time[i] - t$n$Info[ID %in% as.character(tNode), (mTime)],
                    t$e$tDiff[i, neighbour])
    
    return(c(neighbour, xDist, tDiff))
  })
  
  t$f[,"Neighbour" := temp[1,]]
  t$f[,"xDist" := temp[2,]]
  t$f[,"tDiff" := temp[3,]]
  
  return(t)
}

#After simulating the growth of trees by placing recent sequences as tips on a fixed ML tree
growthSim <- function(iT, gFile) {
  #####
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param gFile: The growth file from a pplacer run for all new cases
  #@return: The input tree annotated with growth information stored as $g.
  #    $nID: The ID if the new node. This will be a full tip label
  #    $xDist: The distance between this node and it's the most distant tip in the tree that is formed
  #    $oConn: The index corresponding to the neighbour of the newly added tip. This can be an internal node or a tip.
  .....
  #####
  
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
      oNode <- which(iT$tip.label%in%oIDs)
    }
  
    data.frame(nID=nID, xDist=max(dist.nodes(t)), oConn=oNode)
  }))  
  
  #Add this information as a data table under g.
  iT$g <- as.data.table(df)
  
  return(iT)
}

#Cluster using a modified subtree-based method
#Clusters are defined as subtrees with all tip-to tip distances under some maximum
#A bootstrap criterion can also be applied to these subtrees (ie. confidence in the parent node)
STClu <- function(iT, maxD, minB=0) {
  #####
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param minBL The minimum bootstrap criterion for clusters
  #@param maxD: The maximum distance criterion defining clusters
  #@return: The tree annotated with $c, which contains cluster membership information ($Membership).
  #         This also contains additional information stored in the data frame $c$Info...
  #    $ID: The ID of either the parent node representing the cluster, or the tip (if orig was a singleton)
  #    $Old: The number of cases (not new) in the original cluster
  #    $New: The number of new cases added to the cluster. This represents cluster growth
  #    $mTime: The mean time of the old cases in the cluster 
  ....
  #####
  
  #Identify potential clusters by a maximum internal distance criterion and bootstrap criterion
  #These are tracked in the info data table 
  cluNodes <- iT$n$Info[(BootStrap>=minB) & (xDist<=maxD)]
  cluIDs <- cluNodes[order(totTips, decreasing = T), (ID)]
  iT$n$Info$Clustered <- F
  iT$n$Info[ID%in%cluIDs, "Clustered" := T]
  
  #This loop populates temp with only clusters that are not part of a larger cluster
  temp <- vector()
  while(length(cluIDs)>0) {
    cluID <- cluIDs[1]
    cluDes <- as.character(iT$n$Des[[cluID]])
    cluIDs <- setdiff(cluIDs, c(cluDes,cluID))
    temp <- c(temp, cluID)
  }
  
  #Obtain clusters, with internal nodes now removed
  clus <- iT$n$Des[sort(temp)]
  clus <- lapply(clus, function(c) {c[which(c)<length(iT$tip.label)]})
  
  #Obtain old singletons. These are added as clusters of size 1
  oSing <- as.list(setdiff(1:length(iT$tip.label), unlist(clus)))
  names(oSing) <- unlist(oSing)
  clus <- c(oSing, clus)
  
  #The cluster that each new tip joins is saved as a column in $g
  iT$g$Clu <- sapply(iT$g[,(oConn)], function(x){
    hostClu <- names(clus[sapply(clus, function(clu){x%in%clu})])
    ifelse(length(hostClu)==0, NA, hostClu)
  })
  
  #Correcting for new tips which would create a distance over the maximum clustering distance allowed
  iT$g[(xDist)>maxD, Clu := NA]
  
  #Sort cluster membership. Showing the lists of tip labels for each cluster
  #This is sorted with growth info to be appended onto the final tree
  iT$c <- list()
  
  #Summarizes growth information as a data table
  iT$c$Info <- data.table(ID = names(clus), Old = sapply(clus, function(x) {length(x)}),
                      mTime = c(iT$v[unlist(oSing), (Time)], iT$n$Info[ID%in%temp, (mTime)]))
  

  #This function ensures members are listed by their tip labels.
  #New members are also added here
  iT$c$Membership <- lapply(names(clus), function(cluName){
    clu <- clus[[cluName]]
    oMem <- iT$tip.label[clu]
    nMem <- iT$g[Clu%in%cluName, (nID)]
    return(c(oMem, nMem))
  })

  #Cluster information is finalized and added to $c
  cInfo$New <- sapply(iT$c$Membership, function(x){length(x)}) - cInfo$Old
  
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
  
  #Set inputs for test
  reVars <- '_'
  varInd <- c(1,2)
  dateFormat <- '%Y'
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