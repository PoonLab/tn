#A source file for the necessary function regarding cluster threshold selection
#Does not include graphing or results output

#Import Libraries
library(igraph,verbose = FALSE)
library(dplyr,verbose = FALSE)
library(parallel,verbose = FALSE)
library(ggplot2,verbose = FALSE)
library(R.utils,verbose = FALSE)

#Obtain some frequency data regarding number of linkages from a given year to the newest year in the input graph.
bpeFreq <- function(iG) {
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A data frame of Number of positives (edges from one year to the newest year) 
  #         with total possible edges and time difference (in years) between the two years
  
  #Obtain the range of years to plug into an edge-counting function
  maxY <- max(V(iG)$year)
  minY <- min(V(iG)$year)
  ys <- seq(minY, (maxY-1), 1)
  nV <- V(iG)[year==maxY]  # nodes in more recent year of subgraph
  
  #Runs through each year, counting the number of edges from that year to the newest year. 
  #Also counts the number of possible edges between those two years
  frequency <- sapply(ys, function(x) {
    pV <- V(iG)[year==x]
    bE <- E(iG)[pV%--%nV]  
    pos <- length(bE)
    tot <- length(pV)*length(nV)
    return(c(pos,tot))
  })
  
  #Assign a time difference to each year (the number of years between this year and the newest year)
  tDiff <- sapply(ys, function(x) maxY-x)
  
  #Create a data frame of case attachment frequency and case age
  df <- data.frame(tDiff = tDiff, Positive = frequency[1,], Total = frequency[2,])
  
  return(df)
}

#Filters edges the input graph such that all new cases are only linked to old cases by their closest edge to old cases
minFilt <- function(iG) {
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A filtered version of this same graph, with all new cases holding only one edge to old cases
  
  #Obtain the new year
  nY <- max(V(iG)$year)
  
  #Obtain edge id's of all of the shortest edge lengths from new cases (A new case can only be linked to 1 case)
  bE <- E(iG)[V(iG)[year==nY]%--%V(iG)[year<nY]]
  
  #To catch a case where no new cases link to old ones
  if (length(bE) > 0) {
    
    #Obtain the closest edges for each new case
    cE <- mclapply(V(iG)[year==nY], function(x) {
      xE <- bE[inc(x)]

      #To catch a case that is new, but has no linkages to old cases
      ret <- ifelse(length(xE)==0, 0, (xE[Distance == min(Distance)])[1] )
    })
    
    #Remove the entries from new cases that dont connect to old cases
    cE <- unname(unlist(cE[cE!=0]))
    
    #Filter out all edges except for the closest edges
    if(!is.null(cE)){ iG <- iG - difference(bE, E(iG)[cE]) }
    
  }
  
  return(iG)
}

#Simulates the growth of clusters based on the inputted graph's newest year.
#Clusters are built based on the years previous to the newest year 
#Growth is determined based on the cluster membership of the newest cases.
simGrow <- function(iG) { 
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: The cluster information for that subgraph, annotated with the forecasted and actual growth of each cluster.
  #         Clusters are defined as graph components and forecast is based on a binomial regression model from bpeFreq Data
  
  #Split the input graph into the new cases and present clusters
  nV <- V(iG)[year==max(V(iG)$year)]
  pG <- induced_subgraph(iG, V(iG)[year<max(V(iG)$year)])
  clu <- components(pG)
  
  #Assign cluster growth based on number of new cases linked to old cases in clusters 
  temp <- sapply(1:clu$no, function(x) {
    members <- names(clu$membership[unname(clu$membership)==x])
    memV <- V(iG)[name%in%members]
    bE <- E(iG)[memV%--%nV]
    
    #Frequency is assigned based on bpeFreq function
    forecast <- sum(memV$freq) 
    growth <- length(bE)
    return(c(growth,forecast))
  })
  
  #Assign growth, forecast (based on node date), and incidence
  clu$growth <- temp[1,]
  clu$forecast <- temp[2,]
  clu$inc <- length(nV)
  
  return(clu)
}

#Analyze a given Graph to establish the difference between the performance of two different models
#Performance is defined as the ability for cluster growth to fit a predictive model.
clusterAnalyze <- function(subG) {
  #@param subG: A subGraph cut based on a threshold distance, expecting a member of the multiGraph set
  #@return: Cluster info from simGrow, but annotated with GAIC (difference in performance between 2 models)
  #         also adds the predictive model formula, and the data that the model was based off of.
  
  #Established here for the generation of date-based data
  years <- as.integer(levels(factor(V(subG)$year)), threads)
  
  #Obtain the frequency of edges, for a series of subgraphs cut to each possible year (excluding the newest year)
  ageDi <- bind_rows(lapply(rev(tail(years,-2)), function(y){
    ssubG <- minFilt(induced_subgraph(subG, V(subG)[year<y]))
    bpeFreq(subG)
  }))
 
  #Obtain a model of case connection frequency to new cases as predicted by individual case age
  mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')
  
  #Assign a predicted growth value to each member of the graph based on the model
  V(subG)$freq <- predict(mod, data.frame(tDiff=V(subG)$tDiff), type='response')
  
  #Simulate growth of this subgraph 
  clu <- simGrow(subG)
  
  #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
  df1 <- data.frame(Growth = clu$growth, Pred = clu$forecast)
  df2 <- data.frame(Growth = clu$growth, Pred = clu$csize * (sum(clu$growth)/sum(clu$csize)))
  
  #Determine fit when
  fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
  fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
  
  #Save, gaic, model and age data as part of the output
  clu$gaic <- fit1$aic-fit2$aic
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

  #From the input file, a tn93 output file
  input <- read.csv(infile, stringsAsFactors = F)
  
  #Creates a graph based on the inputted data frame. The tn93 Distances become edge attributes
  g <- graph_from_data_frame(input, directed=F, vertices=NULL)
  
  #Adds the ID's and Sample collection years as different vertex attributes for each vertex
  temp <- sapply(V(g)$name, function(x) strsplit(x, '_')[[1]])
  V(g)$name <- temp[1,]
  V(g)$year <- as.numeric(temp[2,])
  
  #Obtain metaData if provided, fills in missing data as NA
  if (!is.na(metData)) {
    
    #Read from meta-data file
    metD <- read.csv(metData, stringsAsFactors = F)
    
    #Patient, age, sex, risk factor and diagnostic year (normally cases are annotated with sample collection year)
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
  
  #Filters out the newest year if it contains a small number of cases
  #This may also be done if the user selects a value for input filter.  
  while (length(V(g)[year==nY])<63 || inputFilter>0) {
    if (length(V(g)[year==nY])>63){inputFilter <- inputFilter-1}
    nY <- nY-1
  }
  
  #Reset the years based on a newly trimmed graph and obtain the time difference for each node 
  #Years between that node and the newest year
  g <- induced_subgraph(g, V(g)[year<=nY])
  years <- as.integer(levels(factor(V(g)$year)))
  V(g)$tDiff <- sapply(V(g)$year, function(x) nY-x)
  
  return(g)
}
  
#Return a set of subgraphs over some set of cutoff parameters. 
#Graphs are created with createGraph based on user defined funcions
multiGraph <- function(iG, cutoffs, threads) {
  #@param g: The complete graph created with createGraph
  #@param threads: To define how many threads to use for parallel functionality
  #@param cutoffs: A list of cutoffs which we will build graphs based off of
  #@return: A list of multiple graph objects filtered to different cutoffs
  
  #Resolve potential for merging clusters
  iG <- minFilt(iG)
  
  #Create a set of subgraphs at each cutoff, removing any edge with weight above the cutoff
  gs <- mclapply(cutoffs, function(d) {
    #Progress Tracking
    cat(paste0("\r", "Building Graphs  ", d/max(cutoffs)*100, "%"))
    
    subgraph.edges(iG,E(iG)[Distance<=d], delete.vertices = F)
  }, mc.cores=threads) 
  
  return(gs)
}

#Run across a set of several graphs from multiGraph, analyzing GAIC at each with clusterAnalyze
gaicRun <- function(gs, cutoffs, threads) {
  #@param gs: A set of graphs, each created with slightly different parameters
  #@param cutoffs: A list of cutoffs which we will build graphs based off of
  #@param threads: To define how many threads to use for parallel functionality
  #@return: A data frame of each runs cluster information (clusterAnalyze output)
  
  #This script will give warnings because of poor null model fit and lack of convergence
  options(warn=-1)
  
  #Generate cluster data for each graph in gs
  res <- mclapply(cutoffs, function(d) {
    #Progress Tracking
    cat(paste0("\r", "Running Analysis   ", d/max(cutoffs)*100, "%")) 
    
    #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
    subG <- gs[[as.character(d)]]
    
    #Obtain cluster information for this subgraph
    clusterAnalyze(subG)
    
  }, mc.cores=threads)
  
  return(res)
}