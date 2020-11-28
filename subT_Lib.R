require("ape")
require("dplyr")
require("data.table")
require("phytools") #Literally just for midpoint rooting. Use APE for this in future if possible
#require("parallel") #Likely unnecessary

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
  #    $n: A list of each node's descendants ($Des), as well as information used to obtain clusters from nodes 
  #        This additional data is stored in ($Info) as $cDist for the longest branch length between the two child branches.
  #        $cDist will later be correlated with variables such as $tDiff (time difference between nodes).  
  
  #Obtaining and midpioint rooting an ape phylogeny object from the tree file, store in a greater list "t"
  t <- midpoint.root(read.tree(tFile))
  tips <- 1:length(t$tip.label)
  nodes <- (max(tips)+1):(max(tips)*2-1) 
  
  #Obtain lists of sequence ID and Time
  #Reformat edge list as data table object with predictors extracted from sequence header
  temp <- sapply(t$tip.label, function(x) {(strsplit(x, reVars)[[1]])[varInd]})
  t$v <- data.table(ID=temp[1,],  
                    Time=as.Date(temp[2,], format=dateFormat), 
                    stringsAsFactors = F)
  
  #Obtain the full set of descendants at each node
  t$n <- list()
  t$n$Des <- lapply(nodes, function(x){getDescendants(t,x)})
  names(t$n$Des) <- as.character(nodes)
  
  #Information necessary for clustering each node and building a growth model is stored in an info data table
  #See descriptions of the $n output above for more detail on each item
  t$n$Info <- data.table()
  t$n$Info[,"ID" := nodes] 
  t$n$Info[,"TipCount" := sapply(t$n$Des, function(x){length(x[x%in%tips])})]
  t$n$Info[,"mTime" := sapply(t$n$Des, function(x){mean(t$v$Time[(x[x%in%tips])])})]
  t$n$Info[,"Bootstrap" := as.numeric(t$node.label)]
  
  #Set root node bootstrap support to 1 and adjust these values to be fractions out of 1
  #Different tree-building methods display bootstrap support differently
  t$n$Info[is.na(Bootstrap)|((Bootstrap)%in%c("", "Root")), "Bootstrap" := 10^ceiling(log10(max(t$n$Info$Bootstrap[!is.na(t$n$Info$Bootstrap)])))]
  t$n$Info[,"Bootstrap" := (t$n$Info$Bootstrap)/max(t$n$Info$Bootstrap)]
  
  #This block obtains the time difference between each child for each parent node
  #As well as the largest branch length obtained by a specific child branch
  temp <- sapply(t$n$Info[,(ID)], function(x){
    cE <- which(t$edge[,1]==as.numeric(x))
    c <- t$edge[cE,2]
    times <- c(t$v$Time, t$n$Info$mTime)[c]
    c(abs(times[1]-times[2]), max(t$edge.length[cE]))
  })
  t$n$Info[,"tDiff" := temp[1,]]
  t$n$Info[,"cDist" := temp[2,]]
  
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
    
    #A nested loop to obtain placement information for each tip
    bind_rows(lapply(1:length(p), function(i){

      #Extract the subtree which contains the newest tip      
      c <- extract.clade(t, p[i])
      
      #Obtain the node corresponding to the neighbour based on that node's number in the old tree
      oIDs <- setdiff(c$tip.label, nID[i])
      oTips <- which(iT$tip.label%in%oIDs)
      oNode <- getMRCA(iT, oTips)
      
      #To catch a singleton connecting to a new tip
      if(is.null(oNode)){
        oNode <- which(iT$tip.label%in%oIDs)
      }
      
      #Obtain Bootstrap Certainty and clarify node IDs
      #All new node placements should share the same ID
      temp <- strsplit(nID[i], "_#")[[1]]
      ID <- temp[[1]]
      if(length(temp)>1) {
        b <- strsplit(temp[[2]], "=")[[1]]
        b <- as.numeric(b[length(b)])  
      } else {
        b <- 1
      }
      
      data.frame(nID=ID, TermDist=termDist[i], PenDist=penDist[i], oConn=oNode, Bootstrap = b)
    }))
  }))  
  
  #Add this information as a data table under g.
  iT$g <- as.data.table(df)
  
  return(iT)
}

#Clusters are defined as a series of tips diverging from a series of quickly branching nodes
STClu <- function(iT, maxD, minB=0) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param maxD: The maximum distance criterion defining clusters
  #@return: The tree annotated with $c, which contains cluster membership information ($Membership).
  #         This also contains additional information stored in the data frame $c$Info...
  #    $ID: The ID of either the parent node representing the cluster, or the tip (if orig was a singleton)
  #    $Old: The number of cases (not new) in the original cluster
  #    $New: The number of new cases added to the cluster. This represents cluster growth
  
  #Initiate the position of each tip and node
  #An if statement checks if the parent of an edge is the root
  pathStack <- lapply(1:(length(iT$tip.label)+iT$Nnode), function(x) {
    ifelse(x==(length(iT$tip.label)+1), 
           NA, which(iT$edge[,2]==x))
  })
  pathEdge <- sapply(pathStack, function(x){x[[1]]})
  pathLens <- iT$edge.length[pathEdge]
  pathSteps <- which(pathLens<=maxD)

  #Step up any individual edges that are short enough 
  #This loop terminates when either all tips have reached the root
  while(length(pathSteps)>0){
    
    #PathEdge is updated with the parent edge of the parent node for any short original edges
    pathStack[pathSteps] <- lapply(pathStack[pathSteps], function(x){
      pN <- iT$edge[x[[1]],1]
      pE <- ifelse(pN==(length(iT$tip.label)+1), 
             NA, which(iT$edge[,2]==pN))
      return(c(pE, x))
    })
    
    #These two are updated based on pathEdge
    pathEdge <- sapply(pathStack, function(x){x[[1]]})
    pathLens <- iT$edge.length[pathEdge]
    pathSteps <- which(pathLens<=maxD)
  }
  
  #Each tip has the cluster it joins (ie. the highest node reached in an uninterrupted path of short edges) 
  #The same is true for internal nodes, which travel by the same logic
  pathEdge <- sapply(pathStack, function(x){x[[1]]})
  clusIDs <- iT$edge[pathEdge,2]
  clusIDs <- replace(clusIDs, which(is.na(clusIDs)), length(iT$tip.label)+1)
  
  #Check a bootstrap requirement by walking back down the tree
  #Parent Edges with low bootstrap certainty are stepped back to the next point in path which meets bootstrap requirement
  stepDownI <- which(clusIDs%in%iT$n$Info[(Bootstrap)<minB, (ID)])
  rescuedI <- which(clusIDs%in%iT$g[((Bootstrap)>=minB)&((PenDist)<=maxD)&((TermDist)<=maxD), (oConn)])
  stepDownI <- setdiff(stepDownI, rescuedI)
  
  if(length(stepDownI)>1) {
    pathEdge[stepDownI] <- sapply(pathStack[stepDownI], function(x){
      nodes <- iT$edge[x,2]
      cluI <- which(nodes%in%iT$n$Info[Bootstrap>=minB, (ID)])
      ifelse(length(cluI)>0, x[cluI], x[length(x)])
    })  
  }

  clusIDs <- iT$edge[pathEdge,2]
  clusIDs <- replace(clusIDs, which(is.na(clusIDs)), length(iT$tip.label)+1)
  iT$v$Cluster <- clusIDs[1:length(iT$tip.label)]
  iT$n$Info$Cluster <- clusIDs[(length(iT$tip.label)+1):(length(iT$tip.label)+iT$Nnode)]

  #The cluster that each new tip joins is saved as a column in $g
  iT$g[,"Cluster" := c(iT$v$Cluster, iT$n$Info$Cluster)[(oConn)]]
  iT$g[(TermDist>maxD)|(PenDist>maxD), Cluster := 0 ]

  #Sort cluster membership. Showing the lists of tip labels for each cluster
  #This is sorted with growth info to be appended onto the final tree
  iT$c <- list()
  old <- table(iT$v$Cluster)
  
  #Sum bootstrap values within a given cluster
  #Only new cases which join a particular cluster with enough total certainty are considered growth.
  nidG <- sapply(unique(iT$g[,(nID)]), function(id) {
    idG <- iT$g[(nID)%in%id,]
    totB <- sapply(unique(idG[,(Cluster)]), function(x){
      sum(idG[(Cluster)==x, (Bootstrap)])
    })
    ifelse(max(totB)<minB, 0,
           idG[which.max(totB)[[1]], (Cluster)]) 
  })
  temp <- table(nidG)
  
  #Initialize and fill table representing new case attachment (ie. growth)
  new <- old-old
  new[names(new)%in%names(temp)] <- temp[names(temp)%in%names(new)]

  #Summarizes growth information as a data table
  #This excludes any clusters with 
  iT$c$Info <- data.table(ID = as.numeric(names(old)), Old = as.numeric(old), New = as.numeric(new))

  #This function ensures members are listed by their tip labels.
  #New members are also added here
  iT$c$Membership <- lapply(iT$c$Info[,(ID)], function(x){
    iMem <- iT$n$Info[(Cluster)%in%x, ((ID))]
    oMem <- iT$v[(Cluster)%in%x, (ID)]
    nMem <- names(nidG[nidG%in%x])
    return(c(iMem, oMem, nMem))
  })
  
  return(iT)
}

#Obtain GAIC at several different cutoffs
GAICRun <- function(iT, maxDs=NA, minB=0, monitor=T, runID=0) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param maxDs: The maximum distance criteria defining clusters
  #@param minB: The minimum bootstrap criterion for clustering
  #@param runID: An identifier to stash this particular run and compare it to others
  #@param monitor: A switch to determine if there should be a readout of this run
  #@return: A vector of GAIC measurements obtained with different clustering criteria
  
  #Initialize a set of cutoffs to observe (based on the branch-length distribution)
  ##-UNTESTED-##
  if(is.na(maxDs)) {
    steps <- head(hist(iT$edge.length, plot=FALSE)$breaks, -5)
    maxDs <- seq(0 , max(steps), max(steps)/50) 
  }
  
  #This function runs through severel comparisons of a model weighted by predictors, to a model without those variables
  df <- lapply(maxDs, function(d) {
    
    #Obtain clusters
    t <- STClu(iT, d, minB)
    
    #Obtain recency (a sum value of all member tips collection date recency) for each cluster.
    t$c$Info[, "Recency" := sapply(t$c$Membership, function(x) {
      mean(as.numeric(t$v[ID%in%x, (Time)]) - min(as.numeric(t$v[,(Time)])))
    })]
    
    #A randomly weighted model, can also act as a useful comparison
    t$c$Info[, "Rando" := sample(1:00, nrow(t$c$Info), replace=T)]
    
    #Compares clusters with weights obtained through variables to clusters with even weights
    fit1 <- glm(formula = New~Old+Recency, data = t$c$Info, family = "poisson")
    fit2 <- glm(formula = New~Old, data = t$c$Info, family = "poisson")
    fitR <- glm(formula = New~Old+Rando, data = t$c$Info, family = "poisson")

    if(monitor) {
      print(fit1$aic - fit2$aic)
    }
    
    #GAIC is the difference between the AIC of two models
    #Put another way, this is the AIC loss associated with predictive variables
    #Other descriptive data characteristics
    data.frame(modAIC=fit1$aic, nullAIC=fit2$aic, GAIC=(fit1$aic-fit2$aic), randAIC=fitR$aic ,
               GrowthTot=nrow(t$g[Cluster>0,]), Singletons=nrow(t$c$Info[(Old)==1,]), MeanSize=mean(t$c$Info[,(Old)]),
               GrowthMax=max(t$c$Info[,(New)]), GrowthMaxID=t$c$Info[which.max((New)), (ID)],
               SizeMax=max(t$c$Info[,(Old)]), SizeMaxID=t$c$Info[which.max((Old)), (ID)])
  })
  
  dt <- as.data.table(bind_rows(df))
  dt[,"RunID" := runID]
  
  return(dt)
}

#Obtain results for a much larger set of sub-sampled trees.
#See pplacer_utils, for how this data is set up on a larger scale
multiGAICRun <- function(resDir, maxDs, minB=0, 
                         reVars='_', varInd = c(1,2), dateFormat = "%Y") {
  #@param resDir: A directory containing a set of trees previously made as well 
  #               as a matching set of growth file placements
  #@param maxDs: The maximum distance criteria defining clusters
  #
  #@return: A data table of multiple runs worth of cluster information.
  #         Each run will be labelled with a particular run ID
  
  #Obtain a pair of lists - trees and growth files
  tfs <- list.files(paste0(resDir, "trees"))
  tfs <- tfs[order(tfs)]
  gfs <- list.files(paste0(resDir, "growthFiles"))
  gfs <- gfs[order(gfs)]
  
  #Multi GAIC Run based on pre-made directory
  #Again, see pplacer_utils for multi-directory creation
  dt <- bind_rows(lapply(1:length(tfs), function(i){
    
    #Pair tree and growth file from list
    tf <- paste0(resDir, "trees/", tfs[[i]])
    gf <- paste0(resDir, "growthFiles/", gfs[[i]])
  
    print(i)
    
    #Run GAIC run on smaller tree
    sampT <- impTree(tf, reVars, varInd, dateFormat)
    sampT <- growthSim(sampT, gf)
    GAICRun(sampT, maxDs, runID=i, monitor=F)
  }))
  
  return(dt)
}

#Test scripts
if(F){
  
  #Set inputs for test
  reVars <- '_'
  varInd <- c(1,2)
  dateFormat <- "%Y"
  
  tFile <- "~/Data/Seattle/IqTree_Bootstrap/SeattleB_PRO_Filt.fasta.treefile"
  gFile <- "~/Data/Seattle/IqTree_Bootstrap/st.tre"
  
  #tFile <- "~/Data/Tennessee/tn_ColTreeData/tn.refpkg/tn_Coltree.nwk"
  #gFile <- "~/Data/Tennessee/tn_ColTreeData/tn_ColGrowth.tre"
  
  #tFile <- "~/Data/Tennessee/tn_DiagTreeData/tn.refpkg/tn.tre"
  #gFile <- "~/Data/Tennessee/tn_DiagTreeData/tnTGrowth.tre"
  
  #tFile <- "~/Data/Seattle/stKingTrees/stKing_PRO_H_Filt.fasta.treefile"
  #gFile <- "~/Data/Seattle/stKingTrees/stKing.tre"

  #tFile <- "~/Data/NAlberta/IqTree_Bootstrap/na.refpkg/NorthAlbertaB_PRO_Filt.fasta.treefile"
  #gFile <- "~/Data/NAlberta/IqTree_Bootstrap/na.tre"
  
  oT <- impTree(tFile, reVars='_', varInd = c(1,2), dateFormat = "%Y")
  oT <- growthSim(oT, gFile)
     
  maxDs <- seq(0, 0.04, 0.001)
  res <- GAICRun(oT, maxDs, minB=0.90)
  
  plot(maxDs, res$GAIC, xlab="Thresholds", ylab="AIC Loss")
  lines(maxDs, res$GAIC)
  lines(maxDs, res$randAIC-res$nullAIC, col="red")
}