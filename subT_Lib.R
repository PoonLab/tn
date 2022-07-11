require("ape")
require("dplyr")
require("data.table")
require("phytools")
require("parallel") 
require("pROC")


#' Import Tree Data and output an annotated tree with additional information to 
#' assist clustering
#'
#' @param iFile: The name/path of the input file (expecting a newick file)
#' @param rootID: The rootID which can be manually used to root the tree. If NA 
#' - the tree is midpoint rooted
#' @param reVars: The regular expression used to extract variables from column 
#' headers. This is passed to strsplit, creating a vertex of values from the 
#' column header
#' @param varInd: A vector of numbers describing the order of variables in the 
#' split string. This should describe the index of the unique ID, the Timepoint 
#' and the location.
#' ex. If the header would be split such that the 4th index is the Unique ID, 
#' then 4 should be the first number in this list ID and timepoint are currently 
#' required. If the location information is not available, it should be set as 
#' "0".
#' @param addvarN: The names of additional variables beyond the second.
#' @return: An ape phylo object annotated with the additional data summarized 
#' below
#'    $Des: A list of each descendant for each node
#'    $n: Several pieces of info, including $cDist for the longest branch length
#'    between the two child branches.
#'    $cDist: will later be correlated with variables such as $tDiff (time 
#'    difference between nodes).  
import.tree <-function(iFile, varInd=c(1,2), dateFormat="%Y", reVars='_', 
                       addVarN=character(0), rootID=NA, priorityQ=0.80) {
  # Import the truncated tree from the tree file
  # By default, root the tree by using midpoint root, alternatively, the rootID 
  # can be provided.
  t <- read.tree(iFile)
  if(is.na(rootID)) {
    t <- midpoint.root(t)
  } else{
    t <- root(t, outgroup = rootID)
  }
  t <- multi2di(t)
  nodes <- 1:(2*length(t$tip.label)-1)
  
  # Obtain lists of sequence ID and Time
  # Reformat edge list as data table object with predictors extracted from sequence header
  splitHeaders <- sapply(t$tip.label, function(x) {(strsplit(x, reVars)[[1]])[varInd]})
  rownames(splitHeaders) <- c("ID", "Time", addVarN)
  t$seqInfo <- data.table(ID=splitHeaders["ID",],  
                        Time= as.Date(splitHeaders["Time",], format=dateFormat), 
                        stringsAsFactors = F)
  
  if(length(addVarN>0)){
    for(i in 1:length(addVarN)) {
      #Add variables named based on character strings in addVarN
      addVars <- splitHeaders[2+i,]
      
      #To catch the event that the inputted variable is likely a date
      #This is determined if over half of the input can be converted to a date
      #If not a date, type.convert automatically converts the variable typing
      if(sum(sapply(sample(addVars,100,replace =F), function(x) {
        o <- tryCatch({as.Date(x, format=dateFormat)}, 
                      error=function(cond){return(NA)})
        return(is.na(o))
      }))/100<0.5) {
        addVars <- as.Date(addVars, dateFormat)
      } else {
        addVars <- type.convert(addVars)
      }
      
      #Obtain typing, this is better than typeof() in that Date and Factor values are captured
      addVarT <- (strsplit(capture.output(str(addVars)), " |\\[")[[1]])[2]
      
      #Assign variable
      t$seqInfo[, (addVarN[i]):=addVars]
    }
  }
  
  #Information necessary for clustering each node and building a growth model is stored in an info data table
  #See descriptions of the $n output above for more detail on each item
  t$n <- data.table()
  t$n[,"ID" := c(t$seqInfo$ID, nodes[(nrow(t$seqInfo)+1):length(nodes)])] 
  t$n[, "Bootstrap" := sapply(c(rep("-1", nrow(t$seqInfo)), t$node.label), function(x){
    ifelse(grepl("[0-9|/.]", x), as.numeric(x), -1)
  })] 
  t$n[,"Bootstrap" := t$n[,"Bootstrap"]/max(t$n[,"Bootstrap"])]
  t$n[Bootstrap < 0,"Bootstrap" := 1]

  t$n[,"Bootstrap" := (t$n$Bootstrap)/max(t$n$Bootstrap)]
  
  #Obtain the path info from every tip to the root node for future clustering
  #This includes the step length (from nodes to their parents) and bootstrap values of nodes
  depLens <- node.depth.edgelength(t)
  pNodes <- nodepath(t)
  names(pNodes) <- sapply(pNodes, function(p){p[length(p)]})

  #Extend nodepath() to internal nodes
  l <- 0
  while(length(pNodes) < nrow(t$n)+1) {
    newNs <- sapply(pNodes[(l+1):length(pNodes)], function(x){x[-length(x)]})
    newNs <- newNs[match(unique(newNs),newNs)]  
    names(newNs) <- sapply(newNs, function(p){p[length(p)]})
    l <- length(pNodes)
    pNodes <- c(pNodes, newNs[which(!(names(newNs)%in%names(pNodes)))])
  }
  pNodes <- pNodes[which(!sapply(pNodes, function(p) {length(p)})==0)]
  pNodes <- pNodes[order(as.numeric(names(pNodes)))]
  
  pBoots <- lapply(pNodes, function(x){t$n$Bootstrap[x]})
  pLens <- lapply(pNodes, function(x){c((depLens[x[-1]]-depLens[x[-length(x)]]),NA) })
  
  t$pathInfo <- lapply(1:nrow(t$n), function(i) {
    m <- matrix(ncol=length(pNodes[[i]]), nrow=3)
    rownames(m) <- c("Nodes", "Boots", "StepLength")
    m[1,] <- rev(pNodes[[i]])
    m[2,] <- rev(pBoots[[i]])
    m[3,] <- rev(pLens[[i]])
    return(m)
  })
  
  return(t)
}

#More Node stuff
#This is relevent for alternative clustering methods, or results, but takes up time
extendInfo <- function(iT) {
  
  oN <- iT$n
  
  ds <- dist.nodes(iT)
  des <- Descendants(iT, type = "all")
  desDs <- sapply(des, function(x) {ds[x[x<=nrow(iT$seqInfo)],x[x<=nrow(iT$seqInfo)]]})
  
  oN$MaxD <- sapply(desDs, function(x){max(x)})
  oN$MeanD <- sapply(desDs, function(x){mean(x[x>0])})
  oN$Size <- sapply(des, function(x){length(x)})
  oN$des <- des
  
  for(nm in colnames(iT$seqInfo)[-1]){
    oN[,(nm) := lapply(oN$des, function(x){iT$seqInfo[x[x<=nrow(iT$seqInfo)],get(nm)]})]
  }
  oN[, "Membership" := lapply(oN$des, function(x){
    iT$seqInfo[x[x<=nrow(iT$seqInfo)],(ID)]
  })]
  
  return(oN)
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
  i <-1 
  df <- bind_rows(lapply(ts, function(t){
    #print(i)
    i<-i+1
    
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
  g <- as.data.table(df)
  
  return(g)
}

#Clusters are defined as a series of tips diverging from a series of quickly branching nodes
STClu <- function(iT, maxD, minB=0, setID=0) {
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
    c(p[, h], h)
  })
  
  rownames(temp)[4] <- "Height"
  temp["Nodes", is.na(temp["Nodes",])] <- length(iT$tip.label)+1
  
  #Check a bootstrap requirement of these nodes
  stepDownI <- which(temp["Boots",]<minB)
  
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
  iT$g[,"Cluster" := (iT$n$Cluster)[(oConn)]]
  iT$g[(TermDist>maxD)|(PenDist>maxD), Cluster := 0]

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
  clu <- data.table(ID = as.numeric(names(old)), Old = as.numeric(old), New = as.numeric(new))
  clu[,"Membership" := lapply(clu$ID, function(id){c(iT$n[(Cluster)%in%id, (ID)], iT$g[(Cluster)%in%id, (nID)])})]
     
  clu[,"MaxD" := maxD]
  clu[,"MinB" := minB]
  clu[,"SetID" := setID]
  return(clu)
}

#Use cluster-pickers clustering method.
CPClu <- function(iT, maxD, minB=0, setID=0) {
  
  iT$n[,"Clustered":=F]
  iT$n[(MaxD<=maxD)&(Bootstrap>=minB), "Clustered":=T]
  cNodes <- as.numeric(unlist(iT$n[(Clustered), .(des)]))
  iT$n[ID%in%cNodes, "Clustered" := F]

  #Sort out singleton clustering
  iT$n[1:length(iT$tip.label), "Clustered" := F]
  cTips <- unlist(iT$n[(Clustered), .(des)])
  cTips <- as.numeric(cTips[cTips<=nrow(iT$seqInfo)])
  iT$n[1:length(iT$tip.label), "Clustered" := T]
  iT$n[cTips, "Clustered" := F]
  
  #Find parent clusters
  iT$n[,"Cluster":=0]
  for(i in which(iT$n$Clustered)) {
    if(i<=nrow(iT$seqInfo)){
      iT$n[i, "Cluster" := i]
    } else{
      iT$n[iT$n$des[[i]], "Cluster" := iT$n$ID[[i]]] 
    }
  }
  
  #The cluster that each new tip joins is saved as a column in $g
  iT$g[,"Cluster" := (iT$n$Cluster)[(oConn)]]
  maxDs <- iT$n[iT$g$oConn, (maxD)]+iT$g$TermDist+iT$g$PenDist
  iT$g[which(maxDs<maxD), "Cluster" := 0]
  
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
  old <- table(iT$n$Cluster[1:length(iT$tip.label)])
  new <- old-old
  new[names(new)%in%names(temp)] <- temp[names(temp)%in%names(new)]
  
  #Summarizes growth information as a data table
  #This excludes any clusters with 
  clu <- data.table(ID = as.numeric(names(old)), Old = as.numeric(old), New = as.numeric(new))
  clu[,"Membership" := lapply(clu$ID, function(id){c(iT$n[(Cluster)%in%id, (ID)], iT$g[(Cluster)%in%id, (nID)])})]
  
  clu[,"MaxD" := maxD]
  clu[,"MinB" := minB]
  clu[,"SetID" := setID]
  
  return(clu)
}

multiSTClu <- function(iT, maxDs, minBs=0, nCores=1, cluFun=STClu) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param modFormula: The predictive model formula. This may be changed with additional variables
  #                   Recency is always calculated for all clusters.
  #@param maxDs: The maximum distance criteria defining clusters
  #@param minB: The minimum bootstrap criterion for clustering
  #@param runID: An identifier to lable this particular run and compare it to others
  
  #Building all Clusters
  if(length(minBs)==1){
    clus <- bind_rows(mclapply(1:length(maxDs), function(i) {dt <- cluFun(iT, maxD=maxDs[i], minB=minBs, setID=i)}, 
                               mc.cores = nCores))
  } else{
    clus <- bind_rows(mclapply(1:length(minBs), function(i) {dt <- cluFun(iT, maxD=maxDs, minB=minBs[i], setID=i)}, 
                               mc.cores = nCores))
  }
  
  clus[,"Growing" := F]
  clus[(New>0),"Growing" := T]
  
  #Attaching Model Data
  modD <- bind_rows(mclapply(clus$Membership, function(x) {
    mRow <- iT$seqInfo[(ID)%in%x, -1]
    sRow <- lapply(colnames(mRow) ,function(nm) {list(mRow[,get(nm)])})
    dt <- as.data.table(sRow)
    colnames(dt) <- colnames(mRow)
    return(dt)
  }, mc.cores = nCores))
  clus <- cbind(clus, modD)
  
  return(clus)
}

#Obtain GAIC at several different cutoffs
GAICRun <- function(clus, runID=0, nCores=1, modFormula=New~Old+Time,
                    propVar="Time", propTrans=list(function(x){mean(x)})) {
  #@param nCores: Number of cores for parallel processing
  #@return: A data table of each runs cluster information.
  #         Both null and proposed model AIC values, as well as the AIC loss ($nullAIC, $modAIC and $GAIC)
  #         The max size, average size and number of singletons ($SizeMax, $MeanSize and $Singletons)
  #         The total growth, largest growth and number of growing clusters ($GrowthTot, $GrowthMax, and $nGrowing)
  #         The ID of the largest cluster and the cluster with the highest growth ($SizeMaxID and $GrowthMaxID)
  #         The effect ratio of mean recency in growing clusters over mean recency in non-growing clusters ($xMag)
  
  setIDs <- unique(clus$SetID)
  dt <- bind_rows(mclapply(setIDs, function(id) {
    
    dt <- clus[(SetID)==id, c("ID", "Old", "New", "Growing", "MaxD", "MinB")]
    dt[, (propVar) := lapply(1:length(propVar), function(i){
      sapply(clus[(SetID)==id, get(propVar[i])], propTrans[[i]])
    })]
    
    #Compares clusters with weights obtained through variables to clusters with even weights
    fit1 <- glm(formula = modFormula, data=dt, family = "poisson")
    fit2 <- glm(formula = New~Old, data=dt, family = "poisson")
    
    if(length(unique(dt$Growing))<2){
      AUCs <- rep(NA, 1+length(propVar))
    } else {
      tfFormula <- as.formula(do.call("substitute",list(modFormula, list(New=quote(Growing)))))
      fitROC <- roc(formula=tfFormula, data=dt)
      AUCs <- sapply(fitROC, function(x){x$auc})
    }
    res <- data.table(modAIC=fit1$aic, nullAIC=fit2$aic, GAIC=(fit1$aic-fit2$aic),
                      GrowthTot=sum(dt[,(New)]), Singletons=nrow(dt[(Old)==1,]), MeanSize=mean(dt[,(Old)]),
                      GrowthMax=max(dt[,(New)]), GrowthMaxID=dt[which.max((New)), (ID)], Growing=sum(dt$Growing),
                      SizeMax=max(dt[,(Old)]), SizeMaxID=dt[which.max((Old)), (ID)],
                      MaxD=dt[1,(MaxD)], MinB=dt[1,(MinB)])
    nmAUC <- sapply(names(dt[,!c("ID","New", "Growing", "MinB", "MaxD")]), function(s){paste0(s, "AUC")})
    res[, (nmAUC) :=  as.list(AUCs)]
    
    return(res)
  }, mc.cores=nCores))
  
  dt[,"RunID" := runID]
  
  return(dt)
}

#Obtain results for a much larger set of sub-sampled trees.
#See pplacer_utils, for how this data is set up on a larger scale
##-UNTESTED WITH PPLACER UTILS AND ADDVARN-##
multiGAICRun <- function(sampsDir, maxDs, minB=0, nCores=1, reVars='_', varInd = c(1,2), dateFormat = "%Y", addVarN=NA,
                         modFormula=New~Old+Time, propVar="Time", propTrans=list(function(x){mean(x)})) {
  #@param sampsDir: A directory containing a set of trees previously made as well 
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
  sampsDir <- gsub("$|/$", "/", sampsDir)
  short <- (strsplit(list.files(sampsDir, pattern = ".fasta$"), "_")[[1]])[1]
  
  tfs <- sapply(1:length(list.files(sampsDir, pattern = ".fasta$")), function(i) {
    paste0(sampsDir, short, "_refpackages/",
           short, "_refpkg", i, "/",
           short, "_tree", i, ".nwk")
  })
  tfs <- tfs[order(tfs)]
  
  gfs <- sapply(1:length(list.files(sampsDir, pattern = ".fasta$")), function(i){
    paste0(sampsDir, short, "_GrowthFiles/",
           short, "_growth", i, ".tree")
  })
  gfs <- gfs[order(gfs)]
  
  print("Preparing Trees...")
  trees <- mclapply(tfs, function(tf){
    print(tf)
    import.tree(iFile=tf, reVars=reVars, varInd=varInd, dateFormat=dateFormat)
  }, mc.cores=nCores)
  
  print("Simulating Growth...")
  trees <- mclapply(1:length(trees), function(i){growthSim(trees[[i]], gfs[[i]])}, mc.cores=nCores)
  
  #Multi GAIC Run based on pre-made directory
  #Again, see pplacer_utils for multi-directory creation
  print("Running GAIC Analysis...")
  dt <- bind_rows(lapply(1:length(trees), function(i){
    print(i)
    t <- trees[[i]]
    GAICRun(clus)
  }))
  
  return(dt)
}
