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
  #@param nCore: The number of cores used for multi-threading. --CURRENTLY UNNUSED--
  #@return: An ape phylo object annotated with the additional data summarized below
  #    $v: A data frame storing vertex information ($ID, $Time, and, if given $Location)
  #    $e: A list of 2 matrices. $dist representing phenetic distance and $tDiff representing time difference.
  #        Additionally, if location information is provided, shared location is a variable stored in $lMatch. 
  #    $n: A list of each node's descendants ($Des), as well as information used to obtain clusters from nodes 
  #        This additional data is stored in ($Info) as $xDist for the longest distance $Bootstrap for the support value.
  #        also, $mTime for the mean date at which sequences were collected and $mtDiff
  #    $f: A list of data used to train a predictive model. This focuses on instances of terminal branches joining the tree

  
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
  termi <- sapply(tips, function(tip){which(t$edge[,2]==tip)})
  termDist <- t$edge.length[termi]
  termE <- t$edge[termi,]
  pNode <- termE[,1]
  
  t$v[,"TermDist" := termDist]
  t$v[,"ParentNode" := pNode]
  
  t$f <- data.table(Tip = termE[,2], tNode = pNode, termDist = termDist)
  
  #Obtain the neighbour node to each tip (this could be an internal node or another tip)
  #This also stores the time difference and maximum distance given this neighbour node
  temp <- sapply(1:nrow(t$f), function(i){
    tNode <- t$f[i, (tNode)]
    cTips <- t$edge[which(t$edge[,1]%in%tNode),2]
    neighbour <- setdiff(cTips, t$f[i,(Tip)])
    
    #Obtain the pendant length (the length of the edge leading to the neighbour node) and the time difference
    #There are 2 cases for time difference calculation
    penDist <- t$edge.length[which(t$edge[,2]==neighbour)]
    tDiff <- abs(ifelse(neighbour>length(t$tip.label),
                    t$v$Time[i] - t$n$Info[ID %in% as.character(tNode), (mTime)],
                    t$e$tDiff[i, neighbour]))
    
    return(c(penDist, neighbour, tDiff))
  })
  
  t$f[,"penDist" := temp[1,]]
  t$f[,"Neighbour" := temp[2,]]
  t$f[,"tDiff" := temp[3,]]
  
  return(t)
}

#After simulating the growth of trees by placing recent sequences as tips on a fixed ML tree
growthSim <- function(iT, gFile) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param gFile: The growth file from a pplacer run for all new cases
  #@return: The input tree annotated with growth information stored as $g.
  #    $nID: The ID if the new node. This will be a full tip label
  #    $xDist: The distance between this node and it's the most distant tip in the tree that is formed
  #    $oConn: The index corresponding to the neighbour of the newly added tip. This can be an internal node or a tip.
  #    $penDist: The length of the pendant branch (ie. from new node to neighbour)
  
  #Obtain a set of trees with new tips added
  #This is one tree for each new case
  ts <- read.tree(gFile)
  df <- bind_rows(lapply(ts, function(t){
    
    #Obtain the newest tip and it's edge
    nID <- setdiff(t$tip.label, iT$tip.label)
    nTip <- which(t$tip.label%in%nID)
    nE <- which(t$edge[,2]%in%nTip)
    
    #Obtain the parent node created by the new tip
    #As well as the branch length
    p <- t$edge[nE, 1]
    penE <- setdiff(which(t$edge[,1]%in%p), nE)
    termDist <- t$edge.length[nE]
    penDist <- t$edge.length[penE]
    
    #Extract the subtree which contains the newest tip
    t <- extract.clade(t, p)
    
    #Obtain the node corresponding to the neighbour based on that node's number in the old tree
    oIDs <- setdiff(t$tip.label, nID)
    oTips <- which(iT$tip.label%in%oIDs)
    oNode <- mrca.phylo(iT, oTips)
    
    #To catch a singleton connecting to a new tip
    if(is.null(oNode)){
      oNode <- which(iT$tip.label%in%oIDs)
    }
  
    data.frame(nID=nID, TermDist=termDist, PenDist=penDist, oConn=oNode)
  }))  
  
  #Add this information as a data table under g.
  iT$g <- as.data.table(df)
  
  return(iT)
}

#Clusters are defined as a series of tips diverging from a series of quickly branching nodes
STClu <- function(iT, maxD) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param maxD: The maximum distance criterion defining clusters
  #@return: The tree annotated with $c, which contains cluster membership information ($Membership).
  #         This also contains additional information stored in the data frame $c$Info...
  #    $ID: The ID of either the parent node representing the cluster, or the tip (if orig was a singleton)
  #    $Old: The number of cases (not new) in the original cluster
  #    $New: The number of new cases added to the cluster. This represents cluster growth
  
  #Initiate the position of each tip and node
  #An if statement checks if the parent of an edge is the root
  pathEdge <- sapply(1:(length(iT$tip.label)+iT$Nnode), function(x) {
    ifelse(x==(length(iT$tip.label)+1), 
           NA, which(iT$edge[,2]==x))
  })
  pathLens <- iT$edge.length[pathEdge]
  pathSteps <- which(pathLens<=maxD)
  
  #Step up any indivual edges that are short enough 
  #This loop terminates when either all tips have reached the root
  while(length(pathSteps)>0){
    
    #PathEdge is updated with the parent edge of the parent node for any short original edges
    pathEdge[pathSteps] <- sapply(pathEdge[pathSteps], function(x){
      p <- iT$edge[x,1]
      ifelse(p==(length(iT$tip.label)+1), 
             NA, which(iT$edge[,2]==p))
    })
    
    #These two are upated based on pathEdge
    pathLens <- iT$edge.length[pathEdge]
    pathSteps <- which(pathLens<=maxD)
  }
  
  #Each tip has the cluster it joins (ie. the highest node reached in an uninterrupted path of short edges) 
  #The same is true for internal nodes, which travel by the same logic
  pathEdge[is.na(pathEdge)] <- which(iT$edge[,1]==length(iT$tip.label)+1)
  clusIDs <- iT$edge[pathEdge,2]
  iT$v$Cluster <- clusIDs[1:length(iT$tip.label)]
  iT$n$Info$Cluster <- clusIDs[(length(iT$tip.label)+1):(length(iT$tip.label)+iT$Nnode)]

  #The cluster that each new tip joins is saved as a column in $g
  iT$g[,"Cluster" := c(iT$v$Cluster, iT$n$Info$Cluster)[(oConn)]]
  iT$g[(TermDist>maxD)|(PenDist>maxD), Cluster := 0 ]

  #Sort cluster membership. Showing the lists of tip labels for each cluster
  #This is sorted with growth info to be appended onto the final tree
  iT$c <- list()
  old <- table(iT$v$Cluster)
  temp <- table(iT$g$Cluster)
  new <- old-old
  new[names(new)%in%names(temp)] <- temp[names(temp)%in%names(new)]

  #Summarizes growth information as a data table
  #This excludes any clusters with 
  iT$c$Info <- data.table(ID = as.numeric(names(old)), Old = as.numeric(old), New = as.numeric(new))

  #This function ensures members are listed by their tip labels.
  #New members are also added here
  iT$c$Membership <- lapply(iT$c$Info[,(ID)], function(x){
    oMem <- iT$v[(Cluster)%in%x, (ID)]
    nMem <- iT$g[(Cluster)%in%x, (nID)]
    return(c(oMem, nMem))
  })
  
  return(iT)
}

#Obtain GAIC at several different cutoffs
GAICRun <- function(iT, maxDs, minBs=0, rand=F) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param cutoffs: A set of maximum distance criterion defining clusters
  #@return: A list of analysis based on a fitting actual new data to predicted cluster growth
  #         This most importantly includes GAIC
  
  #There are multiple GAIC calculations at different maximum distances. 
  distRun <- lapply(maxDs, function(maxD){
    
    #A calculation at a given maximum distance is made up of several runs at different minimum bootstraps
    bootRun <- lapply(minBs, function(minB){
      
      #Obtain clusters stored in $c 
      subT <- STClu(iT, maxD, minB)
      
      #Check positives now that we have params
      subT$f$Positive <- (subT$f$xDist<=maxD) 
      
      times <- 
      
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
  oT <- STClu(oT, maxD = 0.0)
  
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