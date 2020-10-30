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
  #    $n: A list of each node's descendants ($Des), as well as information used to obtain clusters from nodes 
  #        This additional data is stored in ($Info) as $cDist for the longest branch length between the two child branches.
  #        $cDist will later be correlated with variables such as $tDiff (time difference between nodes).  
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
  
  #Obtain the tip and node names (as numbers)
  #Obtain the list of descendants
  tips <- 1:nrow(t$v)
  nodes <- (max(tips)+1):(max(tips)*2-1) 
  
  #Obtain the full set of descendants at each node
  t$n <- list()
  t$n$Des <- lapply(nodes, function(x){Descendants(t,x,"all")})
  names(t$n$Des) <- as.character(nodes)
  
  #Information necessary for clustering each node and building a growth model is stored in an info data table
  #See descriptions of the $n output above for more detail on each item
  t$n$Info <- data.table()
  t$n$Info[,"ID" := nodes] 
  t$n$Info[,"TipCount" := sapply(t$n$Des, function(x){length(x[x%in%tips])})]
  t$n$Info[,"mTime" := sapply(t$n$Des, function(x){mean(t$v$Time[(x[x%in%tips])])})]
  
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
    iMem <- iT$n$Info[(Cluster)%in%x, ((ID))]
    oMem <- iT$v[(Cluster)%in%x, (ID)]
    nMem <- iT$g[(Cluster)%in%x, (nID)]
    return(c(iMem, oMem, nMem))
  })
  
  return(iT)
}

#Obtain GAIC at several different cutoffs
GAICRun <- function(iT, maxDs) {
  #@param iT: The input tree file, annotated with vertex and edge information
  #@param maxD: The maximum distance criteria defining clusters
  #@return: A vector of GAIC measurements obtained with different clustering criteria
  
  #This function runs through severel comparisons of a model weighted by predictors, to a model without those variab;s
  gaics <- sapply(maxDs, function(d) {
    
    #Obtain clusters
    t <- STClu(iT, d)
    
    #Obtain data used to weight individual tips within a given cluster
    ##TO-DO: Allow for generalization with other predictors -##
    dt <- data.table(tDiff=t$n$Info[,(tDiff)])
    dt[,"Positive" := F]
    dt[(1:nrow(t$n$Info))[which(t$n$Info$cDist<=d)],Positive := T]
    
    #Train model to create a set of weights for each sequence
    mod <- glm(formula = Positive~tDiff, data = dt, family = "binomial")
    weights <- predict(mod, type = 'response', 
                       data.table(tDiff=as.numeric(max(t$v$Time)-c(t$v$Time, t$n$Info$mTime))))
    names(weights) <- c(t$v$ID, t$n$Info$ID)
    
    #Apply those weights to the members of each cluster for cluster predictive models 
    t$c$Info[, "Weight" := sapply(t$c$Membership, function(x) {sum(weights[x])})]
    
    #Compares clusters with weights obtained through variables to clusters with even weights
    fit1 <- glm(formula = New~Weight, data = t$c$Info, family = "poisson")
    fit2 <- glm(formula = New~Old, data = t$c$Info, family = "poisson")
    
    #GAIC is the difference between the AIC of two models
    #Put another way, this is the AIC loss associated with predictive variables
    gaic <- fit1$aic-fit2$aic
  })
  
  return(gaics)
}

#Test scripts (seattle subset)
if(F){
  
  #Set inputs for test
  reVars <- '_'
  varInd <- c(1,2)
  dateFormat <- '%Y'
  tFile <- "~/Data/Seattle/strefpackages/st1.refpkg/sttree1.nwk"
  gFile <- "~/Data/Seattle/strefpackages/growthFiles/growthFile1.tree"
  
  oT <- impTree(tFile, reVars='_', varInd = c(1,2), dateFormat = "%Y") 
  oT <- growthSim(oT, gFile)

  maxDs <- seq(0.001, 0.05, 0.001)
  gaics <- GAICRun(oT, maxDs)
  
  plot(gaics, xlab="Thresholds", ylab="AIC Loss")
  
  
}