# Functions necessary for efficient data storage and processing
library(dplyr)
library(data.table)
library(parallel)
#library(speedglm) --POSSIBLY UNNECESSARY

#Creates a set of data-tables representing a graph of sequences, with the edges between those sequences representing the TN93 Distance.
#The time and location associated with the sequence can be taken either directly from the sequence header, or provided separately in a .csv file
#This set of data tables also includes the set of minimum retrospective edges from sequences at the newest time point.
impTN93 <- function(iFile, reVars='/|\\|', varInd=c(5,6,2), varMan=NA, dateFormat="%Y-%m-%d", partQ=0.95, nCores=detectCores()){
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
  g$e[, "i" := .I]
  g$v[, "i" := .I]
  
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
  g$e[(as.numeric(t1)>=newPoint)&(as.numeric(t2)>=newPoint), "Filtered" := T]

  ### - CURRENTLY REMOVED TO EXPERIMENT WITH THE LACK OF MINIMUM RETROSPECTIVE EDGE REQUIREMENT - ###
  if (F) {
    #Minimum retrospective edges (saved as "f" component of graph object)
    temp <- g$e
    
    #Append reversed edges for the use of a search by first edge
    temp[,c("ID1", "t1")] <- g$e[,c("ID2", "t2")]
    temp[,c("ID2", "t2")] <- g$e[,c("ID1", "t1")]
    if(length(varInd)>2) {
      temp[,c("l1")] <- g$e[,c("l2")]
      temp[,c("l2")] <- g$e[,c("l1")]
    }
    temp$tDiff <- -(g$e$tDiff)
    dt <- rbindlist(list(g$e,temp))
    
    pos <- dt[, list(ID2=.SD$ID2[which.min(.SD[tDiff<0]$Distance)],
                     Distance=min(.SD[tDiff<0]$Distance), 
                     tDiff=.SD$tDiff[which.min(.SD[tDiff<0]$Distance)],
                     lMatch=ifelse((length(varInd)>2), .SD$lMatch[which.min(.SD[tDiff<0]$Distance)], NA)), .(ID1)]
    
    #Exclude any edges that were excluded from initial analysis (ie. TN93 calculation threshold)
    #If locations were not present, lMatch is excluded as well
    g$f <- g$f[!is.na(ID2)]
    
    if(is.na(g$f$lMatch[[1]])) { 
      g$f <- g$f[,c("ID1","ID2", "Distance", "tDiff")]
    }
  }
  
  return(g)
}

#Create clusters based on component clustering by some measure of genetic distance
#If a new case clusters with a known case, this is considered growth
compClu <- function(iG) {
  #@param iG: The inputted graph. Expecting all vertices, but some edges filtered by distance.
  #@return: The inputted graph, annotated with a cluster size summary and case membership in the vertices section
  
  #Simplify the list of unsorted vertices (just id's) and edges (just head and tail id's)
  #These are separated by whether or not they are new
  vid <- iG$v[!(New)]$ID
  adj <- iG$e[!(Filtered)&!(New),c("ID1","ID2")]
  adjN <- iG$e[!(Filtered)&(New),c("ID1","ID2")]
  
  #Initialize the first cluster name and a column for cluster membership.
  iG$v$Cluster <- vector(mode="numeric", length=nrow(iG$v))

  #The search vertex becomes the first member of the first cluster and is removed from the searchable set of cluster names
  i <- 1
  srchV <- vid[i]
  memV <- srchV
  vid <- setdiff(vid, memV)
  growth <- integer(0)
  
  #Assigning Cluster Membership
  repeat {
    
    #Remove edges internal to search query and list outgoing edges
    #Currently slow
    adj <- subset(adj, !(ID1%in%srchV & ID2%in%srchV))
    exE <- subset(adj,(ID1%in%srchV | ID2%in%srchV))
    
    #Find all neighboring vertices to the search vertex (or search vertices) through external edges
    #These are then added to the list of member vertices and removed from the list of searchable vertices
    nbV <- setdiff(c(exE$ID1,exE$ID2), srchV)
    memV <- c(memV, nbV) 
    vid <- setdiff(vid, nbV)
    
    #If there are no more neigbours to the search vertices, the cluster is completed and we reset the search parameters
    if (length(nbV)==0) {
      print(i)
      #Update the growth of this cluster
      memVN <- c(adjN[ID1%in%memV]$ID2,  adjN[ID2 %in% memV]$ID1)
      growth <- c(growth, length(memVN))
      
      #Update the vertices with their cluster membership and update the list of 
      iG$v[(ID%in%memV), "Cluster" := i]

      #The end condition, catching the event that there are no vertices to assign to clusters
      if (length(vid)==0) {break}
      
      #Reset search parameters
      i <- i+1
      srchV <- vid[1]
      memV <- srchV
      vid <- setdiff(vid, memV)
      
      next
    }
    
    #Remove all edges within the current cluster from the adjacency list
    adj <- subset(adj, !(ID1%in%srchV|ID2%in%srchV))
    srchV <- nbV
  }
  
  #Add the overall size of clusters (before calculating their connectivity to new cases).
  #Also add the connectivity to new clusters under the variable "growth"
  iG$g <- data.table(Old = as.numeric(table(iG$v[!(New), "Cluster"])), New = growth)
  
  return(iG)
}
    
#A simple function, removing edges that sit above a maximum reporting distance.
dFilt <- function(iG, maxD) {
  subG <- iG
  subG$e[(Distance>maxD), "Filtered" := T]
  return(iG)
}

#A simple function, removing vertices that do not exist in the range of keepT
tFilt <- function(iG, keepT) {
  iG$v <- iG$v[Time%in%keepT]
  iG$e <- iG$e[tMax%in%keepT]
  return(iG)
}

#Analyze a given Clustered Graph to establish the difference between the performance of two different models
#Performance is defined as the ability for cluster growth to fit a predictive model.
compAnalyze <- function(subG) {
  #@param subG: A subGraph cut based on a threshold distance, expecting a member of the multiGraph set
  #@return: A graph annotated with growth, cluster info and level of predictive performance (measured through GAIC)
  
  #Downsample the number of negative outcomes to 
  posi <- subG$e[!(New) & !(Filtered)]$i
  negi <- subG$e[!(New) & (Filtered)]$i
  negi <- sample(negi, ifelse(length(negi)<=length(posi), length(negi), length(posi)))
  df <-subG$e[c(posi, negi),]
  
  mod <- glm(!Filtered ~ abs(as.numeric(tDiff)), df, family = 'binomial')
  subG$v$Weight <- predict(mod, type='response', data.frame(tDiff=max(subG$v[(New)]$Time)-subG$v$Time))
  subG$g$OW <- subG$v[!(New), sum(.SD$Weight) , by="Cluster"]$V1
  
  #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
  fit1 <- glm(New ~ OW, data = subG$g, family = "poisson")
  fit2 <- glm(New ~ Old, data = subG$g, family = "poisson")

  a <- list()
  a$gaic <- fit1$aic - fit2$aic
  a$mod <- mod
  a$fit1 <- fit1
  a$fit2 <- fit2

  return(a)
}

#Run across a set of several subGraphs created at various filters, analyzing GAIC at each with clusterAnalyze
##-NOT YET REDONE - ##
gaicRun <- function(iG, cutoffs=NA) {
  #@param iG: Expecting the entire Graph, but in some cases may take a subset  
  #@return: A data frame of each runs cluster information (clusterAnalyze output)
  
  #Initialize a set of cutoffs to observe (based on the genetic distance distribution)
  if  (is.na(cutoffs)) {
    steps <- head(hist(subset(iG$e, Distance<0.05)$Distance, plot=FALSE)$breaks, -5)
    cutoffs <- seq(0 , max(steps), max(steps)/50) 
  }

  #A set of several graphs created at different cutoffs
  gs <- lapply(cutoffs, function(d) {dFilt(iG, d)})
  
  #Generate cluster data for each subGraph in gs
  res <- lapply(gs, function(subG) {compAnalyze(subG)})
  names(res) <- cutoffs
  
  return(res)
}

if(F){
  g <- impTN93(iFile="Data/tn93st.txt", reVars='_', varInd=c(1,2), dateFormat="%Y",partQ=1-0.06819591)
  subG <- dFilt(g, 0.015)
}