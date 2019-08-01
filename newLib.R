iFile <- "stD.txt"
#iFile <- "Data/Seattle/tn93StsubB.txt"
iMaxT <- 0
mtD <- NA

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
  
  
  ##TO-DO: Add mtD meta data info here
  #Create data frame of Vertex (ie. Case by case data)
  vl <- unique(data.frame(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), stringsAsFactors=F))
  
  #Order both list elements by time point
  g <- list(v=vl[order(vl$Time),], e=el[order(el$tMax),])
  
  #Filter out newest years for the sake of sample size
  sMaxT <- as.numeric(max(names(which(table(g$v$Time)>63))))
  g <- tFilt(g, sMaxT)
  
  while(iMaxT>0) {
    g <- tFilt(g, max(g$v$Time)-1)
    sMaxT <- as.numeric(max(names(which(table(g$v$Time)>63))))
    g <- tFilt(g, sMaxT)
  }
  
  return(g)
}

#Create clusters based on component clustering by some measure of genetic distance
compClu <- function(iG) {
  #@param iG: A subgraph cut based on threshold distance.
  #@param sing: An option determining whether or not singletons are to be ignored.
  #             If True, singletons are not considered clusters of size one and their growth is not counted
  #@return: The inputted graph, annotated with a cluster size summary and case membership in the vertices section

  #iG <- dFilt(g, 0.02)
  
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
  
  #Labelling Clusters through membership
  repeat {
    
    #Remove edges internal to search query and list outgoing edges
    inE <- adj[intersect(which(adj$ID1%in%srchV),which(adj$ID2%in%srchV)),]
    adj <- setdiff(adj,inE)
    exE <- adj[union(which(adj$ID1%in%srchV),which(adj$ID2%in%srchV)),]

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
    adj <- setdiff(adj,exE)
    srchV <- nbV
  }

  #Add some summary information regarding clusters
  iG$c <- table(iG$v$Cluster)
  
  outG <- inG
  
  return(outG)
}

#Obtain some frequency data regarding number of linkages from a given year to the newest year in the input graph.
bpeFreq <- function(inG) {
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A data frame of Number of positives (edges from one year to the newest year) 
  #         with total possible edges and time difference (in years) between the two years
  
  #Obtain the range of years to plug into an edge-counting function
  maxT <- max(inG$v$Time)
  minT <- min(inG$v$Time)
  ys <- seq(minT, (maxT-1), 1)
  nV <- inG$v$ID[inG$v$Time == maxT] # nodes in more recent year of subgraph
  fromNew <- which(inG$e$ID1%in%nV)
  toNew <- which(inG$e$ID2%in%nV)
  nE <- inG$e[union(toNew, fromNew),]
  
  #Runs through each year, counting the number of edges from that year to the newest year. 
  #Also counts the number of possible edges between those two years
  frequency <- sapply(ys, function(x) {
    pV <- inG$v$ID[inG$v$Time == x]
    fromP <- which(nE$ID1%in%pV)
    toP <- which(nE$ID2%in%pV)
    pos <- length(c(fromP,toP))
    tot <- length(nV) 
    return(c(pos,tot))
  })
  
  #Assign a time difference to each year (the number of years between this year and the newest year)
  tDiff <- sapply(ys, function(x) maxT-x)
  
  #Create a data frame of case attachment frequency and case age
  df <- data.frame(tDiff = tDiff, Positive = frequency[1,], Total = frequency[2,])
  
  return(df)
}


simGrow <- function(inG) {
  
  #inG <- dFilt(tFilt(g, 2012),0.02)
  
  maxT <- max(inG$v$Time)
  inG <- clsFilt(inG)
  
  oldG <- tFilt(inG, (maxT-1))
  newG <- inG
  
  oldClu <- compClu(oldG, sing) 
  oldC <- oldClu[[1]]
  oldG <- oldClu[[2]]
  
  if (sing){
    singL <- oldC$c0
    newG$e <- newG$e[!newG$e$ID1%in%singL,]
    newG$e <- newG$e[!newG$e$ID2%in%singL,]
  }
  newG <- compClu(newG, sing)[[2]]
  
  if (nrow(inG$e[inG$e$tMax,])==0)   {
    growth <- newG$cSum
  } else {
    growth <- newG$cSum - oldG$cSum
  }
  
  oldClu[[2]]$gSum <- growth  
  
  return(oldClu)
  
}

#Analyze a given Graph to establish the difference between the performance of two different models
#Performance is defined as the ability for cluster growth to fit a predictive model.
clusterAnalyze <- function(subG) {
  #@param subG: A subGraph cut based on a threshold distance, expecting a member of the multiGraph set
  #@return: Cluster info from simGrow, but annotated with GAIC (difference in performance between 2 models)
  #         also adds the predictive model formula, and the data that the model was based off of.
  
  #subG <- dFilt(tFilt(g, 2012), 0.015)
  
  #Established here for the generation of date-based data
  times <- as.integer(levels(factor(subG$v$Time)))
  
  #Obtain the frequency of edges, for a series of subgraphs cut to each possible year (excluding the newest year)
  ageDi <- data.frame()
  for (t in rev(tail(times,-2))) {
    ssubG <- clsFilt(tFilt(subG, t))
    ageDi <- rbind(ageDi, bpeFreq(ssubG))
  }
  
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

#Remove edges from some graph that sit above a maximum reporting distance.
dFilt <- function(inG, maxD) {
  #@param inG: The input Graph with all edges present
  #@param maxD: The maximum distance, edges with distance above this are filtered out
  #@return: The input Graph with edges above a maximum reporting distance filtered
  
  inG$e <- inG$e[which(inG$e$Distance<maxD),]
  outG <- inG
  return(outG)
}

#Remove vertices from some graph that sit above a maximum time
tFilt <- function(iG, maxT) {
  #@param iG: The input Graph with all vertices present
  #@param maxT: The maximum distance, edges with Time above this are filtered out
  #@return: The input Graph with vertices above a maximum time filtered out
  
  iG$v <- iG$v[iG$v$Time<=maxT,]
  iG$e <- iG$e[iG$e$tMax<=maxT,]
  
  return(iG)
}

#Filter the edges coming from new cases such that the new cases have no edges to eachother and only one edge leading from them to old cases
#This is a simplification to hep resolve the merging involved in cluster growth and to prevent growth by whole clusters.
clsFilt <- function(inG){
  
  #inG <- dFilt(g, 0.015)
  #inG <- dFilt(tFilt(g, 2012), 0.015)
  
  nV <- inG$v$ID[which(inG$v$Time == max(inG$v$Time))]
  fromNew <- which(inG$e$ID1%in%nV)
  toNew <- which(inG$e$ID2%in%nV)
  withinNew <- intersect(fromNew,toNew)
  
  inG$e <- inG$e[-withinNew,]
  fromNew <- which(inG$e$ID1%in%nV)
  toNew <- which(inG$e$ID2%in%nV)
  
  minE <- sapply(nV, function(v) {
    fromV <- which(inG$e$ID1%in%v)
    toV <- which(inG$e$ID2%in%v)
    incV <- c(fromV,toV)
    
    if (length(incV) > 0){
      vE <- inG$e[incV,]
      minE <- which(vE$Distance == min(vE$Distance))[[1]]
      
      index <- c(fromV,toV)[[minE]]
    } else{
      index <- 0
    }
    
    return(index)
  })
    
  remV <- names(minE[minE==0])
  minE <- unname(minE[minE>0])
  remE <- setdiff(as.numeric(c(fromNew, toNew)), minE)

  outG <- inG
  outG$e <- outG$e[-remE,]
  outG$v <- outG$v[-which(outG$v$ID%in%remV),]

  return(outG)
} 

#Run across a set of several graphs from multiGraph, analyzing GAIC at each with clusterAnalyze
gaicRun <- function(inG, cutoffs) {
  #@param gs: A set of graphs, each created with slightly different parameters
  #@param cutoffs: A list of cutoffs which we will build graphs based off of
  #@param threads: To define how many threads to use for parallel functionality
  #@return: A data frame of each runs cluster information (clusterAnalyze output)
  
  #cutoffs <- seq(0, 0.04, 0.0008)
  
  #This script will give warnings because of poor null model fit and lack of convergence
  options(warn=-1)
  
  gs <- lapply(cutoffs, function(d) {dFilt(inG, d)})
  
  #Generate cluster data for each graph in gs
  res <- lapply(gs, function(subG) {
    clusterAnalyze(subG, sing=T)
  })
  gaics <- sapply(res, function(i){i$gaic})
  plot(gaics)
  
  end_time <- Sys.time()
  
  end_time - start_time
  
  return(res)
}


h <- function(){
  print(memV)
  print(nrow(adj))
  print(length(vid))
}


oC$no
head(sort(oC$csize, T), 200)
length(iG$c)
head(sort(unname(iG$c), T), 200)







