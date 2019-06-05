#Import Libraries
library(igraph,verbose = FALSE)
library(dplyr,verbose = FALSE)
library(parallel,verbose = FALSE)
library(ggplot2,verbose = FALSE)


#Obtain the frequency of edges in a bipartite Graph between two different years as a function of the difference between those years
bpeFreq <- function(iG) {
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A data frame of Number of positives (edges from one year to the newest year) 
  #         with total possible edges and time difference (in years) between the two years
  
  # Obtain the range of years
  maxY <- max(V(iG)$year)
  minY <- min(V(iG)$year)
  ys <- seq(minY, (maxY-1), 1)
  nV <- V(iG)[year==maxY]  # nodes in more recent year of subgraph
  
  # Obtain the frequency of new cases being connected to each year
  frequency <- sapply(ys, function(x) {
    pV <- V(iG)[year==x]  # nodes in older year
    bE <- E(iG)[pV%--%nV]  # bipartite edges
    pos <- length(bE)
    tot <- length(pV)*length(nV)
    return(c(pos,tot))
  })
  
  #Assign age to every case
  tDiff <- sapply(ys, function(x) maxY-x)
  
  #Create a data frame of case attachment frequency and case age
  df <- data.frame(tDiff = tDiff, Positive = frequency[1,], Total = frequency[2,])
  
  return(df)
}

#Filters the input graph such that all new cases are only linked to old cases by their closest edge to old cases
minFilt <- function(iG, home=F) {
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@param home: To define whether or not to use parallel functionality and risk errors
  #@return: A filtered version of this same graph, with all new cases holding only one edge to old cases
  
  #Obtain the new year
  nY <- max(V(iG)$year)
  
  #Obtain edge id's of all of the shortest edge lengths from new cases (A new case can only be linked to 1 case)
  bE <- E(iG)[V(iG)[year==nY]%--%V(iG)[year<nY]]
  
  #To catch a case where no new cases link to old ones
  if (length(bE) > 0) {
    
    #Resolve Parallel problem
    if (home) {
      #Obtain the closest edges for each new case
      cE <- lapply(V(iG)[year==nY], function(x) {
        xE <- bE[inc(x)]
        
        #To catch a case that is new, but has no linkages to old cases
        ifelse(length(xE)==0, 0, (xE[Distance == min(Distance)])[1] ) 
      }, mc.cores=8)
    
    } else {
      #Obtain the closest edges for each new case
      cE <-mclapply(V(iG)[year==nY], function(x) {
        xE <- bE[inc(x)]
        
        #To catch a case that is new, but has no linkages to old cases
        ifelse(length(xE)==0, 0, (xE[Distance == min(Distance)])[1] ) 
      }, mc.cores=8)
    }
    
    #Remove the entries from new cases that dont connect to  old cases
    cE <- unname(unlist(cE[cE!=0]))
    
    #Filter out all edges except for the closest edges
    if(!is.null(cE)){
      iG <- iG - difference(bE, E(iG)[cE])
    }
  }
  
  return(iG)
}

#Obtains the growth of predefined old clusters based on the addition of new clusters
simGrow <- function(iG) { 
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: The cluster information for that subgraph, annotated with the growth of each cluster
  
  #Split the input graph into the new cases and present clusters
  nV <- V(iG)[year==max(V(iG)$year)]
  pG <- induced_subgraph(iG, V(iG)[year<max(V(iG)$year)])
  clu <- components(pG)
  
  #Assign cluster growth based on number of new cases linked to old cases in clusters 
  temp <- sapply(1:clu$no, function(x) {
    members <- names(clu$membership[unname(clu$membership)==x])
    memV <- V(iG)[name%in%members]
    bE <- E(iG)[memV%--%nV]
    forecast <- sum(memV$freq) 
    growth <- length(bE)
    return(c(growth,forecast))
  })
  
  #Assign growth, forecast (based on diagnostic date), and incidence
  clu$growth <- temp[1,]
  clu$forecast <- temp[2,]
  clu$inc <- length(nV)
  
  return(clu)
}

#Analyze a given Graph to establish the difference between the performance of a null model and 
clusterAnalyze <- function(subG) {
  #@param subG: A subGraph cut based on a threshold distance, expecting a member of the multiGraph set
  #@return: A list of cluster information: 
    #A table to show membership, cluster sizes, number of clusters, a proposed model predicting cluster growth
    #In addition: Actual cluster growth and forecasted growth (from simGrow) and gaic (difference in fit between 2 models)
  
  #Established here for ageDi purposes  
  years <- as.integer(levels(factor(V(subG)$year)))
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case ag
  #This data may contain missing cases, hense the complete cases addition
  ageDi <- bind_rows(lapply(rev(tail(years,-2)), function(y){
    ssubG <- minFilt(induced_subgraph(subG, V(subG)[year<y]))
    bpeFreq(subG)
  }))
  
  mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')
  
  #Assign a predicted growth value to each member of the graph
  V(subG)$freq <- predict(mod, data.frame(tDiff=V(subG)$tDiff), type='response')
  
  #Obtain growth based on two models restricted model
  clu <- simGrow(subG)
  
  #Place growth and forecast data in dfs for fit and full growth
  df1 <- data.frame(Growth = clu$growth, Pred = clu$forecast)
  df2 <- data.frame(Growth = clu$growth, Pred = clu$csize * (sum(clu$growth)/sum(clu$csize)))
  
  #Model growth as a function of forecast for fit and full growth models
  mod1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
  mod2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
  
  #Save, gaic, model and age data
  clu$gaic <- mod1$aic-mod2$aic
  clu$mod <- mod
  clu$ageD <- ageDi
  
  return(clu)
}

#Creates a graph based on some inputted run arguments and potentially patient meta-data
createGraph <- function(infile, inputFilter, metData){
  #@param infile: The name/path of the input file (expecting tn93 output)
  #@param inputFilter: Will drop x of the most recent years from the total data set based on this input
  #@param metData: the filename for a dataframe of associated metadata (optional).
  #@return: A graph with each vertex having a name and an id as well as every edge representing some measure of distance.

  #From the input file.
  input <- read.csv(infile, stringsAsFactors = F)
  
  #This script will give warnings due to the fact that there are low fit rates on the null model
  options(warn=-1)
  
  #Creates a graph based on the inputted data frame. The tn93 Distances become edge4 attributes
  g <- graph_from_data_frame(input, directed=F, vertices=NULL)
  
  #Adds the ID's and Sample collection years as different vertex attributes for each vertex
  temp <- sapply(V(g)$name, function(x) strsplit(x, '_')[[1]])
  V(g)$name <- temp[1,]
  V(g)$year <- as.numeric(temp[2,])
  
  #Obtain metaData if provided
  if (!is.na(metData)) {
    metD <- read.csv(metData, stringsAsFactors = F)
    V(g)$Age <- sapply(V(g)$name, function(x){
      if(identical(as.numeric(metD[metD$ID==x, 2]), numeric(0))){NA}
      else {as.numeric(metD[metD$ID==x, 2])}
    })
    V(g)$Sex <- sapply(V(g)$name, function(x){
      if(identical(as.numeric(metD[metD$ID==x, 3]), numeric(0))){NA}
      else {as.numeric(metD[metD$ID==x, 3])}
    })
    V(g)$Risk <- sapply(V(g)$name, function(x){
      metD[metD$ID==x, 4]
    })
    V(g)$year <- sapply(V(g)$name, function(x){
      if(identical(as.numeric(metD[metD$ID==x, 5]), numeric(0))){NA}
      else {as.numeric(metD[metD$ID==x, 5])}
    })
  }
  
  #Filter out empty year data
  g <- induced.subgraph(g, V(g)[!is.na(year)])
  
  #Obtain the range of years and the maximum input year
  years <- as.integer(levels(factor(V(g)$year)))
  nY <- max(years)
  while (length(V(g)[year==nY])<63 || inputFilter>0) {
    nY <- nY-1
    if (length(V(g)[year==nY])>63){inputFilter <- inputFilter-1}
  }
  
  #Reset the years based on a newly trimmed graph and obtain the time difference
  g <- induced_subgraph(g, V(g)[year<=nY])
  years <- as.integer(levels(factor(V(g)$year)))
  V(g)$tDiff <- sapply(V(g)$year, function(x) nY-x)
  
  return(g)
}
  
#Return a set of subgraphsover somne set of cutoffs (parameters)
multiGraph <- function(g, cutoffs, home = F) {
  #@param g: The complete graph created with createGraph
  #@param home: To define whether or not to use parallel functionality and risk errors
  #@param cutoffs: A list of cutoffs which we will build graphs based off of
  #@return: A list of multiple graph objects filtered to different cutoffs
  
  #Resolve potential for merging clusters
  g <- minFilt(g)
  
  #Create a set of subgraphs at each cutoff
  #Avoid parallel functionality (breaks home computer)
  if (home) {
    gs <- mclapply(cutoffs, function(d) {
      subgraph.edges(g,E(g)[Distance<=d], delete.vertices = F)
    }, mc.cores=8) 
  }
  else{
    gs <- mclapply(cutoffs, function(d) {
      subgraph.edges(g,E(g)[Distance<=d], delete.vertices = F)
    }, mc.cores=8) 
  }

  return(gs)
}

#Do a run across a set of several graphs, analyzing GAIC at each
gaicRun <- function(gs, home=F) {
  #@param gs: A set of graphs, each created with slightly different parameters
  #@param home: To define whether or not to use parallel functionality and risk errors
  #@return: A data frame of each runs cluster information (clusterAnalyze output)
  
  if (home) {
    #Generate Growth data
    res <- lapply(cutoffs, function(d) {
      cat(paste0("\r", "Running Analysis ", d/max(cutoffs)*100, "%"))
      
      #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
      subG <- gs[[as.character(d)]]
      
      #Obtain cluster information for this subgraph
      clusterAnalyze(subG)
      
    })
  } else {
    #Generate Growth data
    res <- mclapply(cutoffs, function(d) {
      cat(paste0("\r", "Running Analysis ", d/max(cutoffs)*100, "%"))
      
      #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
      subG <- gs[[as.character(d)]]
      
      #Obtain cluster information for this subgraph
      clusterAnalyze(subG)
      
    }, mc.cores=8)
  }
  
  return(res)
}