# Functions necessary for efficient data storage and processing
library(dplyr,verbose = FALSE)
library(data.table)

#Creates a set of data-tables representing a graph of sequences, with the edges between those sequences representing the TN93 Distance.
#The time and location associated with the sequence can be taken either directly from the sequence header, or provided separately in a .csv file
#This set of data tables also includes the set of minimum retrospective edges from sequences at the newest time point.
impTN93 <- function(iFile, reVars='/|\\|', varInd=c(5,6,2), varMan=NA, dateFormat="%Y-%m-%d", partQ=0.95){
  #@param iFile: The name/path of the input file (expecting tn93 output csv)
  #@param reVars: The regular expression used to extract variables from column headers. This is passed to strsplit, creating a vertex of values from the column header
  #@param varInd: A vector of numbers describing the order of variables in the split string. This should describe the index of the unique ID, the Timepoint and the location.
  #-------------  ex. If the header would be split such that the 4th index is the Unique ID, then 4 should be the first number in this list
  #-------------  ID and timepoint are currently required. If the location information is not available, it should be set as "0".
  #@param varMan: Variables can be assigned manually with a csv containing columns of ID, Time point, and Location, in that order. Again, location is not mandatory. 
  #-------------  If this option is used, reVars and varInd, need not be provided. --CURRENTLY UNNUSED--
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
  
  if((length(varInd)>2) | !is.na(varInd[[3]]) | !is.na(varInd[[3]])==0) {
    el[,"l1" := temp1[varInd[[3]],]]
    el[,"l2" := temp2[varInd[[3]],]]
    el[,"lMatch" := .SD$l1 %in%.SD$l2, .(ID1, ID2)]
  }
  
  #Obtain list of unique sequences (also data table)
  vl <- data.table(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), Location = c(el$l1, el$l2), stringsAsFactors=F)
  vl <- vl[match(unique(vl$ID), vl$ID),]
  
  #Vertices and edges lists together start a graph object
  g <- list(v=vl[order(vl$Time),], e=el[order(el$tMax),])
  
  #Set a "New Point", such that cases after this point are considered new for the purposes of growth
  newPoint <- quantile(as.numeric(g$v$Time), partQ)
  
  #Split into testing and training partitions
  #Label new vertices as inherently New
  g$v[, "New" := F]
  g$v[as.numeric(Time)>=newPoint]$New <- T 
  
  #Remove edges from edge list if they are between two new vertices
  g$e <- g$e[!((ID1%in%g$v$ID[(New)]) & (ID2%in%g$v$ID[(New)]))]
  
  #Minimum retrospective edges (saved as "f" component of graph object)
  temp <- g$e
  
  ##--CURRENTLY ASSUMES LOCATION DATA --- ##
  temp[,c("ID1", "t1", "l1")] <- g$e[,c("ID2", "t2", "l2")]
  temp[,c("ID2", "t2", "l2")] <- g$e[,c("ID1", "t1", "l1")]
  temp$tDiff <- -(g$e$tDiff)
  dt <- rbindlist(list(g$e,temp))
  
  g$f <- dt[, list(ID2=.SD$ID2[which.min(.SD[tDiff<0]$Distance)],
                   Distance=min(.SD[tDiff<0]$Distance), 
                   tDiff=.SD$tDiff[which.min(.SD[tDiff<0]$Distance)],
                   lMatch=ifelse(((length(varInd)>2) | !is.na(varInd[[3]]) | !is.na(varInd[[3]])==0),
                                 .SD$lMatch[which.min(.SD[tDiff<0]$Distance)], NA)), .(ID1)]
  
  #Exclude any edges that were excluded from initial analysis (ie. TN93 calculation threshold)
  #If locations were not present, lMatch is excluded as well
  g$f <- g$f[!is.na(ID2)]
  
  if(is.na(g$f$lMatch[[1]])) { 
    g$f <- g$f[, 1:4] 
  }
  
  return(g)
}

#Create clusters based on component clustering by some measure of genetic distance
###-- CURRENTLY SLOW. Data tables, are fast for a small number of searches into large data sets, not a large number of searches --###
compClu <- function(iG) {
  #@param iG: The inputted graph. Expecting all vertices, but some edges filtered by distance.
  #@return: The inputted graph, annotated with a cluster size summary and case membership in the vertices section
  
  #Simplify the list of unsorted vertices (just id's) and edges (just head and tail id's)
  vid <- iG$v$ID
  adj <- iG$e[,c("ID1","ID2")]
  
  #Initialize the first cluster name and a column for cluster membership. c0 will be reserved for all singletons if sing=T
  iG$v$Cluster <- vector(mode="numeric", length=nrow(iG$v))
  
  #The search vertex becomes the first member of the first cluster and is removed from the searchable set of cluster names
  i <- 1
  srchV <- vid[i]
  memV <- srchV
  vid <- setdiff(vid, memV)

  #Assigning Cluster Membership
  repeat {
    
    #Remove edges internal to search query and list outgoing edges
    adj <- subset(adj, !(ID1%in%srchV & ID2%in%srchV))
    exE <- subset(adj,(ID1%in%srchV | ID2%in%srchV))
    
    #Find all neighboring vertices to the search vertex (or search vertices) through external edges
    #These are then added to the list of member vertices and removed from the list of searchable vertices
    nbV <- setdiff(c(exE$ID1,exE$ID2), srchV)
    memV <- c(memV, nbV) 
    vid <- setdiff(vid, nbV)
    
    #If there are no more neigbours to the search vertices, the cluster is completed and we reset the search parameters
    if (length(nbV)==0) {
      
      iG$v$Cluster[iG$v$ID%in%memV] <- i  
      
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
    adj <- subset(adj, !(ID1%in%srchV | ID2%in%srchV))
    srchV <- nbV
    
  }
  
  #Add some summary information regarding clusters
  iG$c <- table(iG$v$Cluster)
  
  return(iG)
}
    
#A simple function, removing edges that sit above a maximum reporting distance.
dFilt <- function(iG, maxD) {
  iG$e <- iG$e[Distance<=maxD]
  return(iG)
}

#A simple function, removing vertices that do not exist in the range of keepT
tFilt <- function(iG, keepT) {
  iG$v <- iG$v[Time%in%keepT]
  iG$e <- iG$e[tMax%in%keepT]
  return(iG)
}

#Simulate the growth of clusters, showing the difference in cluster size between the newest and the penultimate time point
#The frame of reference for clusters is the penultimate year, simulating one making forcasting decisions based on one time point and validating them with the next
simGrow <- function(iG) {
  #@param: The inputted graph. Expecting all vertices, but some edges filtered by distance.
  #@return: The same cluster annotated with the actual growth and cluster information
  
  #Obtain clusters at the new time point
  nG <- iG
  nSing <- nG$v[!((ID%chin%nG$e$ID1)|(ID%chin%nG$e$ID2))&(New)]
  nG$v <- fsetdiff(nG$v, nSing)
  nG <- compClu(nG)
  
  #obtain clusters at an old time point
  keepT <- min(iG$v$Time):max(iG$v[!(New)]$Time)
  oG <- compClu(tFilt(iG, keepT))
  
  #Define growth as the difference in cluster size between new and old graphs 
  #After clsFilter(), nG will have the same number of clusters as oG, and similar membership
  iG$g <- nG$c-oG$c
  iG$c <- oG$c
  iG$v <- rbind(nG$v, nSing)
  
  return(iG)
}

#Obtains some likelihood data in order to weight cases based on their recency  
likData <- function(iG) {
  #@param iG: The inputted graph. Expecting the entire Graph with new year included.
  #@return: A data frame of "Positives" (related cases) between one time point and another, annotated with the number of possible total positives.
  #         Each positive also carries it's Genetic Distance measurement and time point difference (between time points)
  
  #Take in total graph without the newest time point (otherwise, we include the validation set in or data set)
  keepT <- tail(head(as.numeric(names(table(iG$v$Time))),-1),-1)
  subG <- tFilt(iG, keepT)
  
  #Obtain the closest retrospective edge of every vertex
  f <- bind_rows(lapply(2:nrow(subG$v), function(i){
    
    v <- subG$v[i,]
    
    incE <- subset(iG$e, (ID1%in%v$ID)|(ID2%in%v$ID))
    retE <- subset(incE, (tMax==v$Time)&(tDiff>0))
    minE <- retE[which(retE$Distance==min(retE$Distance))[[1]],]
    df <- data.frame(Distance=minE$Distance, tMax=minE$tMax, tDiff=minE$tDiff, vTotal=nrow(subset(iG$v, Time==v$Time)))
    
    return(df)
    
  }))
  
  return(f)
}

#Analyze a given Clustered Graph to establish the difference between the performance of two different models
#Performance is defined as the ability for cluster growth to fit a predictive model.
compAnalyze <- function(subG) {
  #@param subG: A subGraph cut based on a threshold distance, expecting a member of the multiGraph set
  #@return: A graph annotated with growth, cluster info and level of predictive performance (measured through GAIC)
  
  #Obtain successes (retrospective growth) and attempts (possible retrospective growths)
  tTab <- as.numeric(table(subG$v$Time))
  tDiffs <-(max(subG$f$tMax))-as.numeric(names(table(subG$f$tMax)))
  maxD <- max(subG$e$Distance)
  ageD <- bind_rows(lapply(tDiffs[tDiffs>0] , function(x){
    data.frame(tDiff=x, 
               Positive=nrow(subset(subG$f, tDiff==x & Distance<=maxD)), 
               vTotal=sum(tail(tTab,-x)) )
  }))
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case age
  #Use this to weight cases by age
  mod <- glm(cbind(Positive, vTotal) ~ tDiff, data=ageD, family='binomial')
  subG$v$Weight <- predict(mod, type='response', data.frame(tDiff=max(subG$v$Time)-subG$v$Time))   
  
  #Create clusters for this subgraph and measure growth
  subG <- simGrow(subG)
  cPred <- subset(subG$v, Time<max(Time))[,c("Weight", "Cluster")]

  #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
  df1 <- data.frame(Growth = as.numeric(subG$g), Pred = sapply(names(subG$c), function(x) { sum(subset(cPred, Cluster==as.numeric(x))$Weight) }))
  df2 <- data.frame(Growth = as.numeric(subG$g), Pred = as.numeric(subG$c) * (sum(as.numeric(subG$g))/sum(as.numeric(subG$c))))
  fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
  fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
  
  #For the analysis
  a <- list()
  
  #Save, gaic, model and age data as part of the output
  subG$a$gaic <- fit1$aic-fit2$aic
  subG$a$mod <- mod
  subG$a$propFit <- fit1
  subG$a$nullFit <- fit2
  subG$a$f <- ageD

  return(subG)
}

#Run across a set of several subGraphs created at various filters, analyzing GAIC at each with clusterAnalyze
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
  
}