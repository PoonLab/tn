library(dplyr,verbose = FALSE)

#Creates a set of data-frames representing a graph of sequences, with the edges between those sequences representing the TN93 Distance.
#Sequences must be dated with the date separated from the id by '_'. 
##TO-DO: Currently only accepts year dates. Work to allow more specific dates. 
impTN93 <- function(iFile, mtD, minNS=63){
  #@param iFile: The name/path of the input file (expecting tn93 output csv)
  #@param mtD: the filename for a dataframe of associated metadata (optional).
  #@param minNS: The minimum number of acceptible new Sequences. By default we keep this high.
  #@return: A list of 3 Data frames. An edge list (weighted by TN93 genetic distance), a vertex list, 
  #         and a list of minimum edges, for the future establishment of a timepoint-based model
  
  #Reading input file as a data frame. This will essentially act as an edgelist
  idf <- read.csv(iFile, stringsAsFactors = F)
  temp1 <- sapply(idf$ID1, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp2 <- sapply(idf$ID1, function(x) (strsplit(x,'_')[[1]])[[2]])
  temp3 <- sapply(idf$ID2, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp4 <- sapply(idf$ID2, function(x) (strsplit(x,'_')[[1]])[[2]])
  
  #Create a data frame from the imported edge list. 
  el <- data.frame(ID1=as.character(temp1), t1=as.numeric(temp2), 
                   ID2=as.character(temp3), t2=as.numeric(temp4), 
                   Distance = as.numeric(idf$Distance), stringsAsFactors= F)
  
  #Obtain the maximum time and time difference between the head and tail of each edge
  el$tMax <- pmax(el$t1,el$t2)
  el$tDiff <- abs(el$t1-el$t2)
  
  ##TO-DO: Add mtD meta data info here
  #Create a list of vertices based on the edge list. Original sequences with no edges will not be considered.
  vl <- unique(data.frame(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), stringsAsFactors=F))
  
  #Order edges and vertices by time point and sort them into a single, larger list
  #If the newest timepoint contains a small number of sequences, we remove the newest year from consideration
  g <- list(v=vl[order(vl$Time),], e=el[order(el$tMax),], f=el[order(el$tMax),])
  while(nrow(subset(g$v,Time==max(Time)))<=minNS) {
    g <- tFilt(g, max(g$v$Time)-1)
  }
  
  #Permanently remove edges from the new year such that only the closest edge between new vertices and old vertices remains
  #'g$f' represents the closest edge for all vertices (excluding edges within that vertex's own year)
  g <- clsFilt(g)
  g$f <- likData(g)
  
  return(g)
}

#Create clusters based on component clustering by some measure of genetic distance
compClu <- function(iG) {
  #@param iG: The inputted graph. Expecting all vertices, but some edges filtered by distance.
  #@return: The inputted graph, annotated with a cluster size summary and case membership in the vertices section
  
  #Simplify the list of unsorted vertices (just id's) and edges (just head and tail id's)
  vid <- iG$v[,"ID"]
  adj <- iG$e[,c("ID1","ID2")]
  
  #Initialize the first cluster name and a column for cluster membership. c0 will be reserved for all singletons if sing=T
  iG$v$Cluster <- vector(mode="numeric", length=nrow(iG$v))
  
  #The search vertex becomes the first member of the first cluster and is removed from the searchable set of cluster names
  i <- 1
  srchV <- vid[1]
  memV <- srchV
  vid <- setdiff(vid, memV)
  
  #Assigning Cluster Membership
  repeat {
    
    #Remove edges internal to search query and list outgoing edges
    adj <- subset(adj, !(ID1%in%srchV & ID2%in%srchV))
    exE <- subset(adj, ID1%in%srchV | ID2%in%srchV)

    #Find all neighbouring vertices to the search vertex (or search vertices) through external edges
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

#A simple function, removing edges that sit above a maximum reporting distance (@param:maxD).
dFilt <- function(iG, maxD) {
  iG$e <- subset(iG$e, Distance<=maxD)
  return(iG)
}

#A simple function, removing vertices that sit above a maximum time point (@param: maxT)
tFilt <- function(iG, maxT) {
  iG$v <- subset(iG$v, Time<=maxT)
  iG$e <- subset(iG$e, tMax<=maxT)
  return(iG)
}

#Filter the edges coming from new cases such that the new cases have no edges to eachother and only one edge leading from them to old cases
#This is a simplification to help resolve the merging involved in cluster growth and to prevent growth by whole clusters.
clsFilt <- function(iG){
  #@param iG: The inputted graph. Expecting the entire Graph with new year included.
  #@return: The inputted graph with the new year connected by only its shortest edges per vertex
  
  #Obtain new vertices and remove any internal edges within the new vertices 
  #We are not interested in completely new clustering
  nV <- subset(iG$v, Time==max(Time))$ID
  iG$e <- subset(iG$e, !(t1==max(tMax) & t2==max(tMax))) 
  
  #Obtain edges leading to newer vertices
  nE <- subset(iG$e, tMax==max(tMax))
  
  #Obtain the minimum of these edges for each new vertex
  minEi <- sapply(nV, function(v) {
    incEi <- c(which(nE$ID1%in%v), which(nE$ID2%in%v))
    
    if (length(incEi)>0) {
      incE <- nE[incEi,]
      return(incEi[which(incE$Distance==min(incE$Distance))[[1]]])
    }else{0}
  })
  
  #Obtain minimum distance edges and replace the set of new edges
  minE <- nE[minEi[minEi>0],]
  iG$e <- subset(iG$e, tMax!=max(tMax))
  iG$e <- rbind(iG$e, minE)
  
  return(iG)
} 

#Simulate the growth of clusters, showing the difference in cluster size between the newest and the penultimate time point
#The frame of reference for clusters is the penultimate year, simulating one making forcasting decisions based on one time point and validating them with the next
simGrow <- function(iG) {
  #@param: The inputted graph. Expecting all vertices, but some edges filtered by distance.
  #@return: The same cluster annotated with the actual growth and cluster information
  
  #Obtain clusters at the new time point, after removing singletons
  nG <- iG
  nSing <- subset(nG$v, (!ID%in%c(nG$e$ID1,nG$e$ID2) & Time==max(Time)))
  nG$v <- subset(nG$v, !(!ID%in%c(nG$e$ID1,nG$e$ID2) & Time==max(Time)))
  nG <- compClu(nG)
  
  #obtain clusters at an old time point
  oG <- compClu(tFilt(iG, as.numeric(tail(names(table(iG$v$Time)),2))[[1]]))
  
  #Define growth as the difference in cluster size between new and old graphs 
  #After clsFilter(), nG will have the same number of clusters as oG, and similar membership
  iG$g <- nG$c-oG$c
  iG$c <- oG$c
  
  #Re-Add the singletons, citing new singletons as members of the cluster 0
  if (nrow(nSing)>0){
    nSing$Cluster <- 0
  }
  iG$v <- rbind(nG$v, nSing)
  
  return(iG)
}

#Obtains some likelihood data in order to weight cases based on their recency  
likData <- function(iG) {
  #@param iG: The inputted graph. Expecting the entire Graph with new year included.
  #@return: A data frame of "Positives" (related cases) between one time point and another, annotated with the number of possible total positives.
  #         Each positive also carries it's Genetic Distance measurement and time point difference (between time points)
  
  #Take in total graph without the newest time point (otherwise, we include the validation set in or data set)
  iG <- tFilt(iG, as.numeric(tail(names(table(iG$v$Time)),2))[[1]])
  tTab <- tail(table(iG$v$Time),-1) 
  
  #Obtain subsets of "Positives" (related cases) between one time point and another, annotated with the number of possible total positives.
  #Positives are only considered if, from the perspective of the newest time point, there are no old cases with a closer TN93 distance measurement
  ageD <- lapply(as.numeric(names(tTab)), function(t){ 
    
    #Close filter cases from the current time point t. Treating it as though it is the new time point
    tG <- tFilt(iG, t)
    tdTab <- tail(table(tG$e$tDiff),-1) 
    tG <- clsFilt(tG)
    
    #To catch the event of a complete edgelist with no edges to year t
    if(nrow(tG$e)>0){
      tG$e$Total <- tTab[as.character(t)]
    } else {tG$e$Total <- integer(0)}
    
    #Given the filtered subgraph, record any edges to/from the time point of interest.
    #Notes only the Distance and time difference (between the older point and  the new) 
    #All of tdG is annotated The total number of cases in time t (ie. the total possible edges for that tMax value)
    tdD <- bind_rows(lapply(as.numeric(names(tdTab)), function(td) {
      es <- subset(tG$e, tDiff==td & tMax==t)
      es[,c("Distance", "tMax", "tDiff", "Total")]
    })) 
    
    return(tdD)
  }) 
  
  return(ageD)
}

#Analyze a given Clustered Graph to establish the difference between the performance of two different models
#Performance is defined as the ability for cluster growth to fit a predictive model.
compAnalyze <- function(subG) {
  #@param subG: A subGraph cut based on a threshold distance, expecting a member of the multiGraph set
  #@return: A graph annotated with growth, cluster info and level of predictive performance (measured through GAIC)
  
  #Obtain the frequency of edges between two years based on the time difference between those years
  #Annotate edge information with total possible edges to a given newer year
  dMax <- max(subG$e$Distance)
  
  #Take the total edge frequency data from the graph and format this information into successes and attempts
  #An edge to the newest year falling below the max distance is considered a success
  ageD <- bind_rows(lapply(subG$f, function(t) {
    tDiffs <- as.numeric(levels(as.factor(t$tDiff)))
    pos <- sapply(tDiffs, function(td) {
      nrow(subset(t, tDiff==td & Distance<=dMax))
    })
    data.frame(Positive=pos, Total=t$Total[[1]], tDiff=tDiffs)
  }))
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case age
  #Use this to weight cases by age
  mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageD, family='binomial')
  subG$v$Weight <- predict(mod, data.frame(tDiff=max(subG$v$Time)-subG$v$Time), type='response')
  
  #Create clusters for this subgraph and measure growth
  subG <- simGrow(subG)
  cPred <- subset(subG$v, Time<max(Time))[,c("Weight", "Cluster")]

  #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
  df1 <- data.frame(Growth = as.numeric(subG$g), Pred = sapply(names(subG$c), function(x) { sum(subset(cPred, Cluster==as.numeric(x))$Weight) }))
  df2 <- data.frame(Growth = as.numeric(subG$g), Pred = as.numeric(subG$c) * (sum(as.numeric(subG$g))/sum(as.numeric(subG$c))))
  fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
  fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
  
  #Save, gaic, model and age data as part of the output
  subG$gaic <- fit1$aic-fit2$aic
  subG$mod <- mod
  subG$f <- ageD
  
  return(subG)
}

#Run across a set of several subGraphs created at various filters, analyzing GAIC at each with clusterAnalyze
gaicRun <- function(iG) {
  #@param iG: Expecting the entire Graph, but in some cases may take a subset  
  #@return: A data frame of each runs cluster information (clusterAnalyze output)
  
  #Initialize a set of cutoffs to observe (based on the genetic distance distribution)
  steps <- head(hist(iG$e$Distance, plot=FALSE)$breaks,-5)
  cutoffs <- seq(0 , max(steps), max(steps)/50)
  
  #A set of several graphs created at different cutoffs
  gs <- lapply(cutoffs, function(d) {dFilt(iG, d)})
  
  #Generate cluster data for each subGraph in gs
  res <- lapply(gs, function(subG) {compAnalyze(subG)})
  names(res) <- cutoffs
  
  return(res)
}

