require("ape")
require("dplyr")
require("data.table")
require("phytools")
require("parallel") 

#Import Tree Data and output an annotated tree with additional information to assist clustering
impTree <-function(tFile, reVars="_", varInd=c(1,2), dateFormat="%Y",  addVarN=NA, rootID=NA){
  #@param iFile: The name/path of the input file (expecting a newick file)
  #@paran rootID: The rootID which can be manually used to root the tree. If NA - the tree is midpoint rooted
  #@param reVars: The regular expression used to extract variables from column headers. This is passed to strsplit, creating a vertex of values from the column header
  #@param varInd: A vector of numbers describing the order of variables in the split string. This should describe the index of the unique ID, the Timepoint and the location.
  #               ex. If the header would be split such that the 4th index is the Unique ID, then 4 should be the first number in this list
  #               ID and timepoint are currently required. If the location information is not available, it should be set as "0".
  #@oaram addvarN: The names of additional variables beyond the second.
  #@return: An ape phylo object annotated with the additional data summarized below
  #    $Des: A list of each descendant for each node
  #    $n: Several pieces of info, including $cDist for the longest branch length between the two child branches.
  #        $cDist will later be correlated with variables such as $tDiff (time difference between nodes).  
  
  #Import the truncated tree from the tree file
  #By default, root the tree by using mmidpoint root, alternatively, the rootID can be provided.
  t <- read.tree(tFile)
  if(is.na(rootID)) {
    t <- midpoint.root(t)
  } else{
    t <- root(t, outgroup = rootID)
  }
  nodes <- 1:(2*length(t$tip.label)-1)
  
  #Obtain lists of sequence ID and Time
  #Reformat edge list as data table object with predictors extracted from sequence header
  splitHeaders <- sapply(t$tip.label, function(x) {(strsplit(x, reVars)[[1]])[varInd]})
  seqInfo <- data.table(ID=splitHeaders[1,],  
                        Time= as.Date(splitHeaders[2,], format=dateFormat), 
                        stringsAsFactors = F)
  
  #Obtain the full set of descendants at each node
  t$Des <- lapply(nodes, function(x){getDescendants(t,x)})
  
  #Information necessary for clustering each node and building a growth model is stored in an info data table
  #See descriptions of the $n output above for more detail on each item
  t$n <- data.table()
  t$n[,"ID" := c(seqInfo$ID, nodes[(nrow(seqInfo)+1):length(nodes)])] 
  t$n[,"TipCount" := sapply(t$Des, function(x){length(x[x<=nrow(seqInfo)])})]
  t$n[,"mTime" := sapply(t$Des, function(x){mean(seqInfo$Time[(x[x<=nrow(seqInfo)])])})]
  t$n[,"Bootstrap" := c(rep(100, nrow(seqInfo)), as.numeric(t$node.label))]
    
  #Set root node bootstrap support to 1 and adjust these values to be fractions out of 1
  #Different tree-building methods display bootstrap support differently
  t$n[is.na(Bootstrap), "Bootstrap" := 10^ceiling(log10(max(t$n$Bootstrap[!is.na(t$n$Bootstrap)])))]
  t$n[,"Bootstrap" := (t$n$Bootstrap)/max(t$n$Bootstrap)]
  
  #This block obtains the time difference between each child for each parent node
  #As well as the largest branch length obtained by a specific child branch
  temp <- sapply(t$n[-(1:nrow(seqInfo)),(ID)], function(x){
    cE <- which(t$edge[,1]==as.numeric(x))
    c <- t$edge[cE,2]
    times <- t$n$mTime[c]
    c(abs(times[1]-times[2]), max(t$edge.length[cE]))
  })
  t$n[,"tDiff" := c(rep(NA, nrow(seqInfo)), temp[1,])]
  t$n[,"cDist" := c(rep(NA, nrow(seqInfo)), temp[2,])]
  
  #In the event of additional variables
  if(!is.na(addVarN)){
    
    #Add additional variables
    for(i in 1:length(addVarN)) {
      
      #Add variables named based on character strings in addVarN
      addVars <- splitHeaders[2+i,]
      
      #To catch the event that the inputted variable is likely a date
      #This is determined if over half of the input can be converted to a date
      if(sum(sapply(sample(addVars,100,replace =F), function(x) {
        o <- tryCatch({as.Date(x, format=dateFormat)}, 
                      error=function(cond){return(NA)})
        return(is.na(o))
      }))/100<0.5) {
        addVars <- as.Date(addVars, dateFormat)
      } else {
        
        #If not a date, type.convert automatically converts the variable typing
        addVars <- type.convert(addVars)
      }
      
      #Obtain typing, this is better than typeof() in that Date and Factor values are captured
      addVarT <- (strsplit(capture.output(str(addVars)), " |\\[")[[1]])[2]
      
      #Assign variable
      seqInfo[, (addVarN[i]):=addVars]
      
      #Check Means under each node if applicable
      if(addVarT%in%c("num", "logi", "int", "Date")){
        nm <- paste0("m", addVarN[i])
        t$n[, (nm) := sapply(t$Des, function(x){
          mean(seqInfo[x[x<=nrow(seqInfo)], (addVarN[i])])
        })]
      }
      
      #Check predominant Level if applicable, as well as number of levels under each node
      if(addVarT%in%"Factor"){
        nm1 <- paste0("Dominant", addVarN[i])
        nm2 <- paste0("NumberOf", addVarN[i])
        
        temp <- sapply(t$Des, function(x){
          tb <- table(seqInfo[x[x<=nrow(seqInfo)], get(addVarN[i])])
          c(names(tb)[which.max(tb)], length(tb[tb>0]))
        })
        t$n[,(nm1) := as.factor(temp[1,])]
        t$n[,(nm2) := as.numeric(temp[2,])]
      }
    }
  }
  
  #Obtain the path info from every tip to the root node for future clustering
  #This includes the step length (from nodes to their parents) and bootstrap values of nodes
  depLens <- node.depth.edgelength(t)
  t$pathInfo <- lapply(1:nrow(t$n), function(x){
    nodes <- nodepath(t, x, nrow(seqInfo)+1)
    pathBoot <- t$n$Bootstrap[nodes]
    stepLens <- c((depLens[nodes[-length(nodes)]] - depLens[nodes[-1]]),NA)
    m <- matrix(c(nodes,pathBoot,stepLens), ncol=length(nodes), nrow=3, byrow = T)
    rownames(m) <- c("Nodes", "Boots", "StepLength")
    return(m)
  })
  
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
  
  #Find the Nodes which represent clusters by finding the first step length greater than maxD in the path from node to root
  #Each node then has the cluster it joins (ie. the highest node reached in an uninterrupted path of short edges) 
  temp <- sapply(iT$pathInfo, function(p) {
    h <- which(p["StepLength",]>maxD)[1]
    c(p[,h], h)
  })
  rownames(temp) <- c("Nodes", "Boots", "StepLength", "Height")
  temp["Nodes", is.na(temp["Nodes",])] <- length(iT$tip.label)+1
  
  #Check a bootstrap requirement of these nodes
  #Clusters are "Rescued" if a new node is added such that the internal node above them meets the bootstrap requirement
  stepDownI <- which(temp["Boots",]<minB)
  rescuedI <- which(temp["Nodes",]%in%iT$g[((Bootstrap)>=minB)&((PenDist)<=maxD)&((TermDist)<=maxD), (oConn)])
  stepDownI <- setdiff(stepDownI, rescuedI)
  
  #If the parent node of a cluster does not meet bootstrap requirements, step back down path until it does
  if(length(stepDownI)>1) {
    temp[,stepDownI]  <- sapply(stepDownI, function(i){
      x <- iT$pathInfo[[i]]
      h <- which(x["Boots", 1:temp["Height", i]]>=minB)[1]
      return(c(x[, h], h))
    })
  }

  iT$n$Cluster <- temp["Nodes",]

  #The cluster that each new tip joins is saved as a column in $g
  iT$g[,"Cluster" := c(iT$n$Cluster)[(oConn)]]
  iT$g[(TermDist>maxD)|(PenDist>maxD), Cluster := 0 ]

  #Sort cluster membership. Showing the lists of tip labels for each cluster
  #This is sorted with growth info to be appended onto the final tree
  iT$c <- list()
  old <- table(iT$n$Cluster[1:length(iT$tip.label)])
  
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
    oMem <- iT$n[(Cluster)%in%x, ((ID))]
    nMem <- names(nidG[nidG%in%x])
    return(c(oMem, nMem))
  })
  names(iT$c$Membership) <- iT$c$Info$ID
  
  return(iT)
}

#Obtain GAIC at several different cutoffs
GAICRun <- function(iT, maxDs=NA, minB=0, runID=0, nCores=1, modFormula=(New~Old+Recency)) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param modFormula: The predictive model formula. This may be changed with additional variables
  #                   Recency is always calculated for all clusters.
  #@param maxDs: The maximum distance criteria defining clusters
  #@param minB: The minimum bootstrap criterion for clustering
  #@param runID: An identifier to stash this particular run and compare it to others
  #@param nCores: Number of cores for parallel processing
  #@return: A data table of each runs cluster information.
  #         Both null and proposed model AIC values, as well as the AIC loss ($nullAIC, $modAIC and $GAIC)
  #         The max size, average size and number of singletons ($SizeMax, $MeanSize and $Singletons)
  #         The total growth, largest growth and number of growing clusters ($GrowthTot, $GrowthMax, and $nGrowing)
  #         The ID of the largest cluster and the cluster with the highest growth ($SizeMaxID and $GrowthMaxID)
  #         The effect ratio of mean recency in growing clusters over mean recency in non-growing clusters ($xMag)
  
  
  #Initialize a set of cutoffs to observe (based on the branch-length distribution)
  if(is.na(maxDs)) {
    steps <- head(hist(iT$edge.length, plot=FALSE)$breaks, -5)
    maxDs <- seq(0 , max(steps), max(steps)/50) 
  }
  
  #This function runs through severel comparisons of a model weighted by predictors, to a model without those variables
  df <- mclapply(maxDs, function(d) {
    
    #Obtain clusters
    t <- STClu(iT, d, minB)
    
    #Obtain recency (a sum value of all member tips collection date recency) for each cluster.
    t$c$Info[, "Recency" := sapply(t$c$Membership, function(x) {
      tips <- x[which(x%in%t$n$ID[1:length(t$tip.label)])]
      mean(as.numeric(t$n[ID%in%tips, (mTime)]) - min(as.numeric(t$n[,(mTime)])))
    })]
    
    #A randomly weighted model, can also act as a useful comparison
    t$c$Info[, "Rando" := sample(1:100, nrow(t$c$Info), replace=T)]
    
    #Compares clusters with weights obtained through variables to clusters with even weights
    fit1 <- glm(formula = modFormula, data = t$c$Info, family = "poisson")
    fit2 <- glm(formula = New~Old, data = t$c$Info, family = "poisson")
    fitR <- glm(formula = New~Old+Rando, data = t$c$Info, family = "poisson")
    
    #GAIC is the difference between the AIC of two models
    #Put another way, this is the AIC loss associated with predictive variables
    #Other descriptive data characteristics
    data.frame(modAIC=fit1$aic, nullAIC=fit2$aic, GAIC=(fit1$aic-fit2$aic), randAIC=fitR$aic ,
               GrowthTot=nrow(t$g[Cluster>0,]), Singletons=nrow(t$c$Info[(Old)==1,]), MeanSize=mean(t$c$Info[,(Old)]),
               GrowthMax=max(t$c$Info[,(New)]), GrowthMaxID=t$c$Info[which.max((New)), (ID)],
               SizeMax=max(t$c$Info[,(Old)]), SizeMaxID=t$c$Info[which.max((Old)), (ID)])
  }, mc.cores = nCores)
  
  dt <- as.data.table(bind_rows(df))
  dt[,"RunID" := runID]
  dt[,"MaxDistance" := maxDs]
  dt[,"minBootstrap" := minB]
  
  return(dt)
}

#Obtain results for a much larger set of sub-sampled trees.
#See pplacer_utils, for how this data is set up on a larger scale
##-UNTESTED WITH PPLACER UTILS AND ADDVARN-##
multiGAICRun <- function(resDir, maxDs, minB=0, modFormula=(New~Old+Recency),
                         reVars='_', varInd = c(1,2), dateFormat = "%Y", addVarN=NA) {
  #@param resDir: A directory containing a set of trees previously made as well 
  #               as a matching set of growth file placements
  #@param modFormula: The predictive model formula. This may be changed with additional variables
  #                   Recency is always calculated for all clusters. This is passed to GAICRun()
  #@param maxDs: The maximum distance criteria defining clusters
  #@param minB: The minimum bootstrap criterion for defining clusters
  #@param reVars: The regular expression used to extract variables from column headers. This will be passed to impTree().
  #               The value is then used by strsplit, creating a vertex of values from the column header. 
  #@param varInd: A vector of numbers describing the order of variables in the split string. This should describe the index of the unique ID, the Timepoint and the location.
  #               ex. If the header would be split such that the 4th index is the Unique ID, then 4 should be the first number in this list. This will be passed to impTree()
  #               ID and timepoint are currently required. If the location information is not available, it should be set as "0". 
  #@oaram addvarN: The names of additional variables beyond the second. This will be passed to impTree()
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
    GAICRun(sampT, maxDs, runID=i, modFormula=modFormula)
  }))
  
  return(dt)
}

#Test scripts
if(F){
  
  #Set inputs for test
  reVars <- '_'
  varInd <- c(1,2,3)
  dateFormat <- "%Y-%m-%d"
  addVarN <- "SubType"
  
  tFile <- "~/Data/NAlberta/naFullTree/old.treefile"
  gFile <- "~/Data/NAlberta/naFullTree/old_growth.tree"
  
  oT <- impTree(tFile, reVars='_', varInd = c(1,2,3), dateFormat = "%Y-%m-%d", addVarN = "SubType")
  oT <- growthSim(oT, gFile)
  
  maxDs <- seq(0, 0.04, 0.001)
  res <- GAICRun(oT, maxDs, minB=0.90, nCores=8)
  
  plot(maxDs, res$GAIC, xlab="Thresholds", ylab="AIC Loss", ylim =c(-100, 3))
  lines(maxDs, res$GAIC)
  lines(maxDs, res$randAIC-res$nullAIC, col="red")
  
  #tFile <- "~/Data/Seattle/IqTree_Bootstrap/SeattleB_PRO_Filt.fasta.treefile"
  #gFile <- "~/Data/Seattle/IqTree_Bootstrap/st.tre"
  
  #tFile <- "~/Data/Tennessee/tn_ColTreeData/tn.refpkg/tn_Coltree.nwk"
  #gFile <- "~/Data/Tennessee/tn_ColTreeData/tn_ColGrowth.tre"
  
  #tFile <- "~/Data/Tennessee/tn_DiagTreeData/tn.refpkg/tn.tre"
  #gFile <- "~/Data/Tennessee/tn_DiagTreeData/tnTGrowth.tre"
  
  #tFile <- "~/Data/Seattle/stKingTrees/stKing_PRO_H_Filt.fasta.treefile"
  #gFile <- "~/Data/Seattle/stKingTrees/stKing.tre"

  #tFile <- "~/Data/NAlberta/IqTree_Bootstrap/na.refpkg/NorthAlbertaB_PRO_Filt.fasta.treefile"
  #gFile <- "~/Data/NAlberta/IqTree_Bootstrap/na.tre"
  
  oT <- impTree(tFile, reVars='_', varInd = c(2,3), dateFormat = "%Y-%m-%d")
  oT <- growthSim(oT, gFile)
     
  maxDs <- seq(0, 0.04, 0.001)
  res <- GAICRun(oT, maxDs, minB=0.90, nCores=8, modFormula=New~Old+Recency+DominantSubType)
  
  plot(maxDs, res$GAIC, xlab="Thresholds", ylab="AIC Loss")
  lines(maxDs, res$GAIC)
  #lines(maxDs, res$randAIC-res$nullAIC, col="red")
}