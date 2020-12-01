# Functions necessary for efficient data storage and processing
library(dplyr)
library(data.table)
library(parallel)
#library(speedglm) --POSSIBLY UNNECESSARY

#Creates a set of data-tables representing a graph of sequences, with the edges between those sequences representing the TN93 Distance.
#The time and location associated with the sequence can be taken either directly from the sequence header, or provided separately in a .csv file
#This set of data tables also includes the set of minimum retrospective edges from sequences at the newest time point.
impTN93 <- function(iFile, reVars='/|\\|', varInd=c(5,6,2), varMan=NA, dateFormat="%Y-%m-%d", partQ=0.95){
  #@param iFile: The name/path of the input file (expecting tn93 output csv)
  #@param reVars: The regular expression used to extract variables from column headers. This is passed to strsplit, creating a vertex of values from the column header
  #@param varInd: A vector of numbers describing the order of variables in the split string. This should describe the index of the unique ID, the Timepoint and the location.
  #               ex. If the header would be split such that the 4th index is the Unique ID, then 4 should be the first number in this list
  #               ID and timepoint are currently required. If the location information is not available, it should be set as "0".
  #@param varMan: Variables can be assigned manually with a csv containing columns of ID, Time point, and Location, in that order. Again, location is not mandatory. 
  #               If this option is used, reVars and varInd, need not be provided. --CURRENTLY UNNUSED--
  #@param partQ: The proportion of the set that is to define "known" cases for the purposes of cluster growth. The remaining quantile is marked as new cases.
  #@return: A list of 3 Data frames. An edge list (weighted by TN93 genetic distance), a vertex list, 
  #         and a list of minimum edges, for the future establishment of a timepoint-based model
  
  #Obtain tn93 edge list from file
  idt <- fread(iFile)
  
  #Reformat edge list as data table object with predictors extracted from sequence header
  temp1 <- sapply(idt$ID1, function(x) (strsplit(x,reVars)[[1]]))
  temp2 <- sapply(idt$ID2, function(x) (strsplit(x,reVars)[[1]]))

  el <- data.table(ID1=as.character(temp1[varInd[[1]],]), t1=as.Date(temp1[varInd[[2]],], format=dateFormat),
                   ID2=as.character(temp2[varInd[[1]],]), t2=as.Date(temp2[varInd[[2]],], format=dateFormat),
                   Distance = as.numeric(idt$Distance), stringsAsFactors= F)
  
  #Obtain the maximum time and time difference between the head and tail of each edge
  el[,"tMax" := pmax(t1,t2)] 
  el[,"tDiff" := (t1-t2)]
  
  if(length(varInd)>2) {
    el[,"l1" := temp1[varInd[[3]],]]
    el[,"l2" := temp2[varInd[[3]],]]
    el[,"lMatch" := .SD$l1 %in%.SD$l2, .(ID1, ID2)]
  }
  
  #Obtain list of unique sequences (also data table)
  vl <- data.table(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), Location = c(el$l1, el$l2), stringsAsFactors=F)
  vl <- vl[match(unique(vl$ID), vl$ID),]
  
  #Vertices and edges lists together start a graph object
  g <- list(v=vl[order(vl$Time),], e=el[order(el$tMax),])
  g$e[, "I" := .I]
  g$v[, "I" := .I]
  
  #Set a "New Point", such that cases after this point are considered new for the purposes of growth
  newPoint <- quantile(as.numeric(g$v$Time), partQ)
  
  #Split into testing and training partitions
  #Label new vertices as inherently New, new edges are those which contain a new vertex
  g$v[, "New" := F]
  g$e[, "New" := F]
  g$v[as.numeric(Time)>=newPoint, "New" := T]  
  g$e[as.numeric(tMax)>=newPoint, "New" := T]
  
  #This is the complete List of edges, however, edges are filtered by default if they lie between two new vertices
  #In the future, more edges may be filtered based on a distance cutoff requirement
  g$e[,"Filtered" := F]
  
  #Obtain minimum retrospective edges for new cases
  retE <- g$e[(New)&((t1<newPoint)|(t2<newPoint)),]
  iMRE <- sapply(g$v[(New), (ID)], function(id){
    iE <- retE[(ID1%in%id)|(ID2%in%id), ]
    iE[which.min(iE$Distance)[[1]], I] 
  })
  
  #Only the minimum retrospective edges for all new cases remain unfiltered
  g$e[(New), "Filtered" := T]
  g$e[iMRE, "Filtered" := F]
  
  return(g)
}

#Create clusters based on component clustering by some measure of genetic distance
#If a new case clusters with a known case, this is considered growth
compClu <- function(iG, maxD) {
  #@param iG: The inputted graph. Expecting all vertices, but some edges filtered by distance.
  #@param maxD: The maximum tn93 distance. Edges higher than this distance will be filtered out of the graph
  #@return: The inputted graph, annotated with a cluster size summary and case membership in the vertices section
  
  iG <- copy(iG)
  
  #Filter cases distant edges
  iG$e[(Distance>maxD), "Filtered" := T]
  
  #Simplify the list of unsorted vertices (just id's) and edges (just head and tail id's)
  #These are separated by whether or not they are new
  vid <- iG$v[!(New), (ID)]
  adj <- as.matrix(iG$e[!(Filtered)&!(New),c("ID1","ID2")])
  adjN <- as.matrix(iG$e[!(Filtered)&(New),c("ID1","ID2")])
  
  #Initialize the first cluster name and a column for cluster membership.
  iG$v$Cluster <- vector(mode="numeric", length=nrow(iG$v))

  #The search vertex becomes the first member of the first cluster and is removed from the searchable set of cluster names
  ci <- 1
  srchV <- vid[ci]
  memV <- srchV
  vid <- setdiff(vid, memV)
  growth <- integer(0)
  
  #Assigning Cluster Membership
  repeat {
    
    #Remove edges internal to search query and list outgoing edges
    adj <- adj[!((adj[,"ID1"]%in%srchV) & (adj[,"ID2"]%in%srchV)),,drop=F]
    exE <- adj[((adj[,"ID1"]%in%srchV) | (adj[,"ID2"]%in%srchV)),,drop=F]
    
    #Find all neighboring vertices to the search vertex (or search vertices) through external edges
    #These are then added to the list of member vertices and removed from the list of searchable vertices
    nbV <- setdiff(c(exE[,"ID1"],exE[,"ID2"]), srchV)
    memV <- c(memV, nbV) 
    vid <- setdiff(vid, nbV)
    
    #If there are no more neigbours to the search vertices, the cluster is completed and we reset the search parameters
    if (length(nbV)==0) {
      
      #Update the growth of this cluster
      memVN <- c(adjN[adjN[,"ID1"]%in%memV,"ID2"], adjN[adjN[,"ID2"]%in%memV,"ID1"])
      growth <- c(growth, length(memVN))
      
      #Update the vertices with their cluster membership and update the list of 
      iG$v[(ID)%in%memV, "Cluster" := ci]

      #The end condition, catching the event that there are no vertices to assign to clusters
      if (length(vid)==0) {break}
      
      #Reset search parameters
      ci <- ci+1
      srchV <- vid[1]
      memV <- srchV
      vid <- setdiff(vid, memV)
      
      next
    }
    
    #Remove all edges within the current cluster from the adjacency list
    adj <- adj[!((adj[,"ID1"]%in%srchV) | (adj[,"ID2"]%in%srchV)),,drop=F]
    srchV <- nbV
  }
  
  #Add the overall size of clusters (before calculating their connectivity to new cases).
  #Also add the connectivity to new clusters under the variable "growth"
  iG$c <- list()
  iG$c$Membership <- lapply(1:ci, function(x){iG$v[(Cluster)==x, (ID)]})
  iG$c$Info <- data.table(Old = as.numeric(table(iG$v[!(New), (Cluster)])), New = growth)
  
  return(iG)
}

#Run across a set of several subGraphs created at various filters, analyzing GAIC at each with clusterAnalyze
GAICRun <- function(iG, maxDs=NA, runID=0, nCores=1) {
  #@param iG: Expecting the entire Graph, but in some cases may take a subset  
  #@param maxDs: A list of cutoff thresholds
  #@param runID: An identifier to stash this particular run and compare it to others
  #@param nCores: The number of cores for parallel functionality
  #@return: A data frame of each runs cluster information.
  
  #Initialize a set of cutoffs to observe (based on the genetic distance distribution)
  if  (is.na(maxDs)) {
    steps <- head(hist(subset(iG$e, Distance<0.05)$Distance, plot=FALSE)$breaks, -5)
    maxDs <- seq(0 , max(steps), max(steps)/50) 
  }
  
  #This function runs through severel comparisons of a model weighted by predictors, to a model without those variables
  df <- mclapply(maxDs, function(d) {
    
    #Obtain clusters
    subG <- compClu(iG, d)
    
    #Obtain recency (a sum value of all member tips collection date recency) for each cluster.
    subG$c$Info[, "Recency" := sapply(subG$c$Membership, function(x) {
      mean(as.numeric(subG$v[ID%in%x, (Time)]) - min(as.numeric(subG$v[,(Time)])))
    })]
    
    #Compares clusters with weights obtained through variables to clusters with even weights
    fit1 <- glm(formula = New~Old+Recency, data = subG$c$Info, family = "poisson")
    fit2 <- glm(formula = New~Old, data = subG$c$Info, family = "poisson")
    
    #GAIC is the difference between the AIC of two models
    #Put another way, this is the AIC loss associated with predictive variables
    #Other descriptive data characteristics
    data.frame(modAIC=fit1$aic, nullAIC=fit2$aic, GAIC=(fit1$aic-fit2$aic),
               GrowthTot=sum(subG$c$Info[,(New)]), Singletons=nrow(subG$c$Info[(Old)==1,]), MeanSize=mean(subG$c$Info[,(Old)]),
               GrowthMax=max(subG$c$Info[,(New)])[[1]], GrowthMaxID=which.max(subG$c$Info[,(New)]),
               SizeMax=max(subG$c$Info[,(Old)]), SizeMaxID=which.max(subG$c$Info[,(Old)]))
  }, mc.cores=nCores)
  
  dt <- as.data.table(bind_rows(df))
  dt[,"RunID" := runID]
  
  return(dt)
}

#Obtain several run results for a much larger set of sub-sampled graphs
multiGAICRun <- function(iG, n, maxDs=NA, prop=0.80) {
  #@param iG: Expecting the entire Graph  
  #@param maxDs: A list of cutoff thresholds
  #@param n: The number of subsample runs to determine
  #@param prop: The proportion of the data set sampled for each sub-sample.
  #@return: A data table of multiple runs worth of cluster information.
  #         Each run will be labelled with a particular run ID

  #Initialize a set of cutoffs to observe (based on the genetic distance distribution)
  if  (is.na(maxDs)) {
    steps <- head(hist(subset(iG$e, Distance<0.05)$Distance, plot=FALSE)$breaks, -5)
    maxDs <- seq(0 , max(steps), max(steps)/50) 
  }
  
  #Run n different runs, each labelled with their own runID
  dt <- bind_rows(lapply(1:n, function(i){
    
    #Copy a sample instance of the inserted graph
    sampG <- copy(iG)
    
    #Sample new and old IDs 
    sampIDs <- c(sample(iG$v[!(New),(ID)], round(prop*nrow(iG$v[!(New)]))),
                 iG$v[(New), (ID)])
    
    sampG$v <- iG$v[(ID)%in%sampIDs,]
    sampG$e <- iG$e[((ID1)%in%sampIDs)&((ID2)%in%sampIDs),]
    
    print(i)
    
    GAICRun(sampG, maxDs, runID=i)
  }))
  
  return(dt)
}