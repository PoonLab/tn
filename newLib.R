iFile <- "stD.txt"
#iFile <- "Data/Seattle/tn93StsubB.txt"
iMaxT <- 0
mtD <- NA
g <- impTN93(iFile,iMaxT, mtD)

#Creates a pair of of data frames based on TN93 output files.
#One data frame represents an edgelist or case interactions, while the other represents the cases themselves
impTN93 <- function(iFile, iMaxT, mtD){
  #@param iFile: The name/path of the input file (expecting tn93 output csv)
  #@param iFilt: Will manually drop x of the most recent time points from the total data set based on this input
  #@param mtD: the filename for a dataframe of associated metadata (optional).
  #@return: A pair of data frames representing edge and vertex information for a distance-based graph
  
  #From the input file, a tn93 output file. This
  idf <- read.csv(iFile, stringsAsFactors = F)
  temp1 <- sapply(idf$ID1, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp2 <- sapply(idf$ID1, function(x) (strsplit(x,'_')[[1]])[[2]])
  temp3 <- sapply(idf$ID2, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp4 <- sapply(idf$ID2, function(x) (strsplit(x,'_')[[1]])[[2]])
  
  #Create data frame of edges (ie. Vertex interactions)
  el <- data.frame(ID1=as.character(temp1), t1=as.numeric(temp2), ID2=as.character(temp3), t2=as.numeric(temp4), 
                   Distance = as.numeric(idf$Distance), stringsAsFactors= F)
  el$tMax <- pmax(el$t1,el$t2)
  el$tDiff <- abs(el$t1-el$t2)
  
  ##TO-DO: Add mtD meta data info here
  #Create data frame of Vertex (ie. Case by case data)
  vl <- unique(data.frame(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), stringsAsFactors=F))
  
  #Order both list elements by time point
  g <- list(v=vl[order(vl$Time),], e=el[order(el$tMax),], f=el[order(el$tMax),])
  
  #Filter out newest years for the sake of sample size
  sMaxT <- as.numeric(max(names(which(table(g$v$Time)>63))))
  g <- tFilt(g, sMaxT)
  
  while(iMaxT>0) {
    g <- tFilt(g, max(g$v$Time)-1)
    sMaxT <- as.numeric(max(names(which(table(g$v$Time)>63))))
    g <- tFilt(g, sMaxT)
  }
  
  #Close Filter the overall graph at this point to save future time complexity
  #g <- clsFilt(g)
  
  #Save a copy of the complete list of minimum edges
  g$f <- bpeFreq(g)
  
  return(g)
}

#Create clusters based on component clustering by some measure of genetic distance
compClu <- function(iG) {
  #@param iG: A subgraph cut based on threshold distance and .
  #@param sing: An option determining whether or not singletons are to be ignored.
  #             If True, singletons are not considered clusters of size one and their growth is not counted
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
  iG$f <- subset(iG$f, Distance<=maxD)
  return(iG)
}

#A simple function, removing vertices that sit above a maximum time point (@param: maxT)
tFilt <- function(iG, maxT) {
  iG$v <- subset(iG$v, Time<=maxT)
  iG$e <- subset(iG$e, tMax<=maxT)
  iG$f <- subset(iG$f, tMax<=maxT)
  return(iG)
}

#Filter the edges coming from new cases such that the new cases have no edges to eachother and only one edge leading from them to old cases
#This is a simplification to hep resolve the merging involved in cluster growth and to prevent growth by whole clusters.
clsFilt <- function(iG){
  #@param iG: A subgraph cut based on threshold distance at a specific year.
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

#Simulate the growth of clusters, showing the difference in cluster size between two years (with the same set of clusters)
simGrow <- function(iG) {
  #@param iG: The complete subgraph cut based on threshold distance
  #@return: The same cluster annotated with the actual growth
  
  #Obtain clusters at the new time point, after removing singletons
  nG <- iG
  nG$v <- subset(nG$v, !(!ID%in%c(nG$e$ID1,nG$e$ID2) & Time==max(Time)))
  nG <- compClu(nG)
  
  #obtain clusters at an old time point
  oG <- compClu(tFilt(iG, max(iG$v$Time)-1))
  
  #Define growth as the difference in cluster size between new and old graphs (after min filtering)
  iG$g <- nG$c-oG$c
  
  return(iG$g)
}

#Obtain some frequency data regarding number of linkages from a given year to the newest year in the input graph.
bpeFreq <- function(iG) {
  #@param iG: The total graph
  #@return: A data frame counting number of new year connections, number of possible new year connections and time difference (from cases to new year)
  
  #Take in total graph without the newest year
  iG <- tFilt(iG, max(iG$v$Time)-1)
  iG$e <- subset(iG$e, tDiff>0)
  
  #Close filter all vertices in the total graph
  minEi <- sapply(iG$v$ID, function(v) {
    
    #Obtain indices representing edges incident on the vertex
    incEi <- c(which(iG$e$ID1%in%v), which(iG$e$ID2%in%v))
    incE <- iG$e[incEi,]
    
    #Obtain the minimum of the edges found by incVi or assign 0 for new singletons
    if (nrow(incE)>0){incEi[which(incE$Distance==min(incE$Distance))[[1]]]
    } else {0}
  })
  
  #Obtain the minimum edges external edges of the total edgelist
  f <- iG$e[unname(minEi[minEi>0]),]
  
  return(f)
}

############Point of Review

#Analyze a given Graph to establish the difference between the performance of two different models
#Performance is defined as the ability for cluster growth to fit a predictive model.
clusterAnalyze <- function(subG) {
  #@param subG: A subGraph cut based on a threshold distance, expecting a member of the multiGraph set
  #@return: A graph annotated with growth, cluster info and level of predictive performance (measured through GAIC)
  
  #subG <- dFilt(g, 0.015)
  
  #Obtain the frequency of edges, for a series of subgraphs cut to each possible year (excluding the newest year)
  ageDi <- table(subG$f$tDiff)
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case age
  mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')
  
  #Assign a predicted growth value to each member of the graph based on the model
  subG$v$Freq <- predict(mod, data.frame(tDiff=max(subG$v$Time)-subG$v$Time), type='response')
  
  #Create clusters for this subgraph, make predictions and measure growth
  clu <- simGrow(subG, sing)
  cSize <- clu[[2]]$cSum
  cPred <- sapply(tail(clu[[1]],-1), function(c) {sum(subG$v$Freq[which(subG$v$ID%in%c)])})
  cGrowth <- clu[[2]]$gSum
  
  if (length(cGrowth)>0){
    #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
    df1 <- data.frame(Growth = cGrowth, Pred = cPred)
    df2 <- data.frame(Growth = cGrowth, Pred = cSize * (sum(cGrowth)/sum(cSize)))
    
    #Determine fit when
    fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
    fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
    
    #Save, gaic, model and age data as part of the output
    clu$gaic <- fit1$aic-fit2$aic
    
  } else {
    clu$gaic <- 0
  }

  clu$mod <- mod
  clu$ageD <- ageDi
  
  return(clu)
}

#Run across a set of several graphs from multiGraph, analyzing GAIC at each with clusterAnalyze
gaicRun <- function(inG, cutoffs) {
  #@param gs: A set of graphs, each created with slightly different parameters
  #@param cutoffs: A list of cutoffs which we will build graphs based off of
  #@param threads: To define how many threads to use for parallel functionality
  #@return: A data frame of each runs cluster information (clusterAnalyze output)
  
  cutoffs <- seq(0, 0.04, 0.0008)
  
  #This script will give warnings because of poor null model fit and lack of convergence
  options(warn=-1)
  
  gs <- lapply(cutoffs, function(d) {dFilt(inG, d)})
  
  #Generate cluster data for each graph in gs
  res <- lapply(gs, function(subG) {
    clusterAnalyze(subG, sing=T)
  })
  gaics <- sapply(res, function(i){i$gaic})
  plot(gaics)
  
  return(res)
}






