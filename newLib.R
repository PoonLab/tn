infile <- "stD.txt"
#infile <- "Data/Seattle/tn93StsubB.txt"
inputFilter <- 0
metData <- NA

#Creates a graph based on some inputted run arguments and potentially patient meta-data
createGraph <- function(infile, inputFilter, metData){
  #@param infile: The name/path of the input file (expecting tn93 output)
  #@param inputFilter: Will drop x of the most recent years from the total data set based on this input
  #@param metData: the filename for a dataframe of associated metadata (optional).
  #@return: A graph with each vertex having a name and an id as well as every edge representing some measure of distance.
  
  #From the input file, a tn93 output file. This will be an edgeList
  input <- read.csv(infile, stringsAsFactors = F)
  temp1 <- sapply(input$ID1, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp2 <- sapply(input$ID1, function(x) (strsplit(x,'_')[[1]])[[2]])
  temp3 <- sapply(input$ID2, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp4 <- sapply(input$ID2, function(x) (strsplit(x,'_')[[1]])[[2]])
  
  #Represents edge data and vertex data as two seperate dataframes
  el <- data.frame(ID1=as.character(temp1), t1=as.numeric(temp2), ID2=as.character(temp3), t2=as.numeric(temp4), 
                   Distance = as.numeric(input$Distance), stringsAsFactors= F)
  el$tDiff <- abs(el$t1 - el$t2)
  el$tMax <- sapply(1:nrow(el), function(i){max(el[i,]$t1, el[i,]$t2)})
  
  vl <- unique(data.frame(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), stringsAsFactors=F))
  
  ##TO-DO: Meta-Data (potentially broken in cluLib)
  
  #Arrange into a named list representing the graph
  g <- list(v=vl[order(vl$Time),], e=el[order(el$tMax),])
  
  return(g)
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

clusters <- function(inG) {
  
  #inG <- dFilt(g, 0.015)
  
  #To set up a column in vector info of cluster membership
  inG$v$Cluster <- vector(mode="numeric", length=nrow(inG$v))
  
  #Simplify the list of vertices (just id's) and edges (just head and tail id's)
  vid <- inG$v[,"ID"]
  adj <- inG$e[,c("ID1","ID2")]
  
  #Initialize the first cluster name and the list of clusters, with c0 reserved for all singletons
  i <- 1
  clu <- list("c0"=vector())
  
  #Initialize the search term, our first vertex to base the clustering off of. 
  #This also becomes the first member of this cluster
  search <- vid[1]
  member <- search
  
  #Because the first search vertex has been sorted into a cluster, we remove it from temp 
  vid <- setdiff(vid, member)
  repeat {

    #Find vertices from the search set and to the search set but not within the search set
    fromSearch <- which(adj$ID1%in%search)
    toSearch <- which(adj$ID2%in%search)
    withinSearch <- intersect(fromSearch, toSearch)
    fromSearch <- setdiff(fromSearch,withinSearch)
    toSearch <- setdiff(toSearch,withinSearch)

    #Find all neighbouring vertices to the search vertex (or search vertices)
    #These are added to the current cluster and removed from temp
    neighbours <- unique(c(adj$ID2[fromSearch],adj$ID1[toSearch]))
    member <- c(member, neighbours) 

    vid <- setdiff(vid, member)

    #If there are no more neigbours to the search vertex, the cluster is completed and we reset the search parameters
    if (length(neighbours)==0) {
      
      #To catch a singleton
      if (length(member)==1) {
        clu[["c0"]] <- c(clu[["c0"]], member)
      }
      else {
        clu[[paste0("c",i)]] <- member
        inG$v$Cluster[which(inG$v$ID%in%member)] <- i
        i <- i+1
      }
      
      #The end condition, catching the event that there are no vertices to assign to clusters
      if (length(vid)==0) {break}
      
      #Reset search parameters
      search <- vid[1]
      member <- search
      vid <- vid[-which(vid%in%search)]

      next
    }
    
    #Remove all edges within the current cluster from the adjacency list
    adj <- adj[-c(withinSearch,fromSearch,toSearch),]
    search <- neighbours
  }

  #Add some summary information regarding clusters
  inG$cSum <- sapply(tail(clu,-1), function(x){length(x)})
  inG$cNo <- length(clu)-1
  
  outG <- inG
  
  return(list(clu, outG))
}

simGrow <- function(inG) {
  
  #inG <- tFilt(g, 2012)
  
  maxT <- max(inG$v$Time)
  inG <- clsFilt(inG)
  
  oldG <- tFilt(inG, (maxT-1))
  newG <- inG
  
  oldClu <- clusters(oldG) 
  oldC <- oldClu[[1]]
  oldG <- oldClu[[2]]
  
  sing <- oldC$c0
  newG$e <- newG$e[-which(newG$e$ID1%in%sing),]
  newG$e <- newG$e[-which(newG$e$ID2%in%sing),]
  
  newG <- clusters(newG)[[2]]
  
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
  clu <- simGrow(subG)
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
  }
  else {clu$gaic <- 0}

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
tFilt <- function(inG, maxT) {
  #@param inG: The input Graph with all vertices present
  #@param maxT: The maximum distance, edges with Time above this are filtered out
  #@return: The input Graph with vertices above a maximum time filtered out
  
  inG$v <- inG$v[inG$v$Time<=maxT,]
  inG$e <- inG$e[inG$e$tMax<=maxT,]
  
  outG <- inG
  
  return(outG)
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
    
  minE <- unname(minE[minE>0])
  remE <- setdiff(as.numeric(c(fromNew, toNew)), minE)
  
  outG <- inG
  outG$e <- outG$e[-remE,]

  return(outG)
} 

#######WIP
#Run across a set of several graphs from multiGraph, analyzing GAIC at each with clusterAnalyze
gaicRun <- function(inG, cutoffs) {
  #@param gs: A set of graphs, each created with slightly different parameters
  #@param cutoffs: A list of cutoffs which we will build graphs based off of
  #@param threads: To define how many threads to use for parallel functionality
  #@return: A data frame of each runs cluster information (clusterAnalyze output)
  
  #cutoffs <- seq(0, 0.03, 0.002)
  
  #This script will give warnings because of poor null model fit and lack of convergence
  options(warn=-1)
  
  gs <- lapply(cutoffs, function(d) {dFilt(inG, d)})
  
  #Generate cluster data for each graph in gs
  res <- lapply(gs, function(subG) {clusterAnalyze(subG)})
  
  return(res)
}


for (i in gs) {print(nrow(i$e))}
gaics <- sapply(res, function(i){i$gaic})




