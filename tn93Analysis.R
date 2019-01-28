#A process process for scoring the effectiveness of clustering from tn93 Output.
#USAGE: Rscript tn93Graph.R tn93output.csv

library(igraph)
####- TO-DO: Use ML for variable selection (see lasso) and better fitting (see Caret or e107) -####

#____________________________________________________________________________________________________________________________#

#Obtains the Growth of clusters based on which clusters haold new cases
growG <- function(inG, full=F) {
  #@param inG: A subG gaph cut based on a threshold distance, with the latest casses representing New cases (ie. Upcoming cases)
  #@param full: An option determining whether or not this is the growth estimate for a fully saturated model
  #@return: Cluster information for the Present year (ie. The year before the newest year in inG)
  
  #Define vertices as new cases and present cases
  newV <- V(inG)[year==max(year)]
  presV <- V(inG)[-newV]
  
  #Used to obtain ther fully saturated model. With case growth measured based on no aggregation in present cases
  if (full){
    bridgeE <- E(inG)[newV%--%V(inG)[-newV]]
    inG <- inG - E(inG)[-bridgeE]
    newV <- V(inG)[year==max(year)]
    presV <- V(inG)[-newV]
  }
  
  #Obtain cluster information
  clu <- components(inG)
  
  #Define the incidence and initialize the growth
  clu$inc <- length(newV)/(length(presV) + length(newV))
  clu$growth <- integer(clu$no) 
  
  #Obtain cluster growth (ie. the number of new cases in each cluster, based on clusters defined at the latest year)  
  gcTable <- table(unname(clu$membership[names(clu$membership) %in% newV$name])) 
  gcIds <- as.integer(names(gcTable)) 
  gcMag <- unname(gcTable)
    
  clu$growth[gcIds] <- gcMag  

  return(clu)
}

#Estimates the growth of clusters based on information from recent years
####- TO-DO: Include include meta-data factors -####
forecast <- function(clu) {
  #@param clu: A set of clusters based on the present year.
  #@return: An attribute for clu representing the past growth of clusters relative to their size. (ie. The predicted absolute growth)
  
  #Obtain a past year to compare too (ideally 5 years before the present year), to establish the recent growth of present clusters
  presY <- max(V(presG)$year)
  minY <- min(V(presG)$year)
  
  #To ensure we don't exceed the bottom limit of the years in data.
  if (presY>(minY+5)) {oldY <- (presY-5)} 
  else {oldY <- minY}
  
  #Difference in past year and present year
  diff <- presY-oldY
  
  #Obtain the cluster sizes of present clusters based only on membership from the old year
  oldMem <- clu$membership[names(clu$membership) %in% V(presG)[year<=oldY]$name]
  oldCsize <- sapply(1:clu$no, function(x) length(oldMem[unname(oldMem)==x]))
  
  #Obtain the Relative, Recent Growth of clusters
  rrG <- (clu$csize-oldCsize) / (diff*sqrt(clu$csize))  
  
  return(rrG)
}

#Simulates growth of current clusters by adding clusters from the upcoming year, with the newest year provided considered "the upcoming" year
growthSim <- function(inG, alt=F, full=F) {
  #@param inG: A subG gaph cut based on a threshold distance, with the latest casses representing New cases (ie. Upcoming cases)
  #@param alt: An option determining whether or not to resolve cluster growth as weighted nodes or closest nodes
  #@param full: An option determining whether or not this is the growth estimate for a fully saturated model
  #@return: Cluster information for the Present year (ie. The year before the newest year in inG)

  #Obtain the present graph (ie. the subgraph contain all cases before the the latest/"upcoming" year) and it's cluster info
  newV <- V(inG)[year==max(year)]
  presG <- inG - newV
  
  #If the full option is selection, we de-aggregate the present year. All vertices are considered clusters of size 1
  if (full) {presG <- presG - E(presG)}
  
  clu <- components(presG)  #Obtain initial cluster info (cluster membership, cluster sizes, number of clusters)
  clu$growth <- integer(clu$no)  #Initialize cluster growth (ie. The appearence of new cases in each cluster) at 0 for each cluster
  clu$inc <- 0   #Initialize cluster incidence (ie. The total number of new cases added over the total number of present cases) at 0
  
  #Establishes the connections between present and new years to measure growth
  bridgeE <- E(inG)[newV%--%V(inG)[-newV]]
  
  #If there are no bridging Edges between present and new cases, there is no growth and the initializeed values represent the cluster growth
  if (length(bridgeE)>0) {
    
    #Establishes a bipartite graph with only the edges that link new cases to present cases
    bridgeG <- subgraph.edges(inG, bridgeE, delete.vertices=T)
    newBridgeV <- V(bridgeG)[year==max(year)]
    presBridgeV <- V(bridgeG)[-newBridgeV]
    
    #Calculate incidence based on the new cases linked to present cases
    clu$inc <- length(newBridgeV) / length(V(presG))
    
    #Branches here to differ between 2 different methods of resolving new cases which bridge 2 clusters
    if(alt){
      
      #Obtain edge id's of all of the shortest edge lengths (A case can only be linked to 1 case)
      closeE <- unname(sapply(newBridgeV, function(x) {
        xE <- E(bridgeG)[inc(x)]
        closest <- xE[Distance == min(Distance)]
        return (closest[[1]])
      }))
      
      #Filter bridge again to only include the shortest edges from new cases
      bridgeG <- subgraph.edges(bridgeG, E(bridgeG)[closeE], delete.vertices=T)
      newBridgeV <- V(bridgeG)[year==max(year)]
      
      #Obtain a table of values, representing the growth of given cluster indices
      growth <- unname(sapply(newBridgeV, function(x) {
        n <- neighbors(bridgeG, x, "all")
        cluId <- unname(clu$membership[names(clu$membership) %in% n$name])
        return(c(cluId, 1))
      }))
        
    } else {
      
      #Obtain the weight of a future case based on the number of past cases it's added to
      weight <- sapply(newBridgeV, function(x) 1/length(E(bridgeG)[inc(x)])) 
      
      #Obtain a table of every link from past to future cases and the weight that those links are worth
      growth <- unname(sapply(presBridgeV, function(x){
        point <- sum(unname(weight[neighbors(bridgeG, x, "all")$name])) 
        cluId <- unname(clu$membership[V(bridgeG)[x]$name])
        return (c(cluId, point))
      }))
    }
    
    #Assign growth as an attribute of the set of clusters based off of the weighted information from growth 
    clu$growth <- sapply(seq(1,clu$no), function(x) round(sum(growth[2,][growth[1,]==x])))
  }
  
  return(clu)
}

#Obtains a filtered subgraph of the full graph. Vertices are removed beyond a given year and edges are removed below a cutoff
subGraph <- function(inG, y, d, plot=F) {
  #@param y: The year that represents the latest year. We forward-censor everything past this.
  #@param d: The distance that represents the cutoff threshold. We remove all edges above this.
  #@param plot: Creates an easy to view subgraph for the purposes of plotting and overview, but not for statistical analysis.
  #@return: The filtered graph (forward censored and cut by a given distance)
  
  #Removes vertices beyond a current year
  outV <- V(inG)[V(inG)$year>y]
  outG <- inG - outV
  
  #Removes edges with distances above a certain cutoff
  outE <- E(outG)[E(outG)$Distance>=d]
  outG <- outG - outE
  
  #Plot option ignores clusters of size 1 and provides a graph (for ease of overview, not for calculations)
  if(plot){ 
    outG <- subgraph.edges(outG, E(outG), delete.vertices = T)
    plot(outG, vertex.size = 2, vertex.label = NA, edge.width = 0.65, edge.color = 'black')
    print(components(outG))
  } 
  
  return(outG)
}

#Obtains several different statistics from a variety of different modelling options
stats <- function(clu) {
  #@param clu: A set of cluster information annotated with incidence and cluster growth at a given year and cutoff
  #@return: A list of statistics based on the provided set of clusters
  
  #Obtain an expectation based off of overall cluster growth
  exp <- clu$inc*(clu$csize)
  
  #zscore can be calculated based off of this diff
  diff <- clu$growth - exp
  zscore <- diff / sqrt(exp)
  
  #Obtain a prediction based on the function forecasting cluster growth 
  pred <- forecast(clu)
  
  #Initialize the list of statistics
  stats <- list()
  
  #Generate data frame at level of clusters
  df <- data.frame(Expectation=exp, Size=clu$csize, Growth=clu$growth, RelGrowth=clu$growth/clu$csize, Prediction=pred)
  
  stats$acc <- glm(RelGrowth~Prediction, data=df, family="poisson") #Models growth against prediction 

  stats$zscore <- mean(zscore) #The mean zscore (representing probability based on expection)
  
  stats$prob <- glm(Growth~Expectation, data=df, family = "poisson") #Models growth against calculated expectation
  
  stats$GbyS <- glm(Growth~Size, data=df, family="poisson") #Models growth by size
  
  stats$Sby1 <- glm(Size~1, data=df, family="poisson") #Models how well size follows a poisson distribution
  
  stats$Gby1 <- glm(Growth~1, data=df, family="poisson") #Models how well growth follows a poisson distribution
  
  stats$RGby1 <- glm(RelGrowth~1, data=df, family="poisson") #Models how well relative growth follows a poisson distribution
  
  stats$var <- var(clu$growth/clu$csize) #The variance of relative growth
  
  return(stats)
}

#Obtains statistics at a set of cutoffs for the same input graph (initially fully connected)
####-TO-DO: Confirm comparisons between Seattle results to select actual methods -####
analyzeG <- function(cutoffs,year=max(years), inG=g) {
  #@param cutoffs: The cutoffs sequence, this will become the independant variable for all results
  #@param year: The year representing the newest (upcoming) year. Ideally the latest year
  #@param inG: The graph to be worked upon. Ideally the complete graph with all connections
  #@return: A dataframe of all calculated results
  
  #Initialize the dataframe of results
  res <- {} 
  
  for (d in cutoffs) {
    print(d) #Important for progress tracking
    
    #Obtain and analyze growth statistics for a given cutoff
    subG <- subGraph(inG,year,d)
    growth <- growG(subG) ##COMPARE: growth <- growthSim(subG) OR growthSim(subG, alt=T)
    fit <- stats(growth)
    
    #Obtain a fully saturated model for the same growth. This will have no cluster aggregation, instead tracking individual new case links
    growth <- growG(subG, full=T) ##COMPARE: growth <- growthSim(subG, full=T) OR growthSim(subG, alt=T, full=T)
    full <- stats(growth) 
    
    #Obtain the Forecast Accuracy
    #This is a measure of how well the relative growth is predicted by the forecast function
    ####- TO-DO: Confirm that this does not need a fully saturated comparison -####
    Accuracy <- 1-(fit$acc)$deviance/(fit$acc)$null.deviance
    
    #Obtain the mean zscore for growths 
    #This is a measure of how unlikely growth rates are based on a theoretical distribution
    zscoreMean <- fit$zscore
    
    #Obtain Deviance and GAICC for a Poisson Probability Map
    #This is a measure of how well the probabilities of growth rates follow a poisson distribution
    DeviancePPMap <- (full$prob)$deviance - (fit$prob)$deviance ##COMPARE: DeviancePPMap <- (fit$prob)$null.deviance - (fit$prob)$deviance
    GAICPPMap <- (fit$prob)$aic- (full$prob)$aic
    
    #This is a measure of how well the growth is predicted by size
    DevianceGbyS <- (full$GbyS)$deviance - (fit$GbyS)$deviance ##COMPARE: DevianceGbyS <- (fit$GbyS)$null.deviance - (fit$GbyS)$deviance
    GAICGbyS <- (fit$GbyS)$aic- (full$prob)$aic
    
    #This is a measure of how well size follows a poisson distribution
    DevianceSby1 <- (full$Sby1)$deviance - (fit$Sby1)$deviance ##COMPARE: DevianceSby1 <- (fit$GbyS)$null.deviance - (fit$GbyS)$deviance
    GAICSby1 <- (fit$Sby1)$aic- (full$Sby1)$aic

    #This is a measure of how well growth follows a poisson distribution    
    DevianceGby1 <- (full$Gby1)$deviance - (fit$Gby1)$deviance ##COMPARE: DevianceGby1 <- (fit$GbyS)$null.deviance - (fit$GbyS)$deviance
    GAICGby1 <- (fit$Gby1)$aic- (full$Gby1)$aic

    #This is a measure of how well Relative growth follows a poisson distribution
    DevianceRGby1 <- (full$RGby1)$deviance - (fit$RGby1)$deviance ##COMPARE: DevianceRGby1 <- (fit$GbyS)$null.deviance - (fit$GbyS)$deviance
    GAICRGby1 <- (fit$RGby1)$aic- (full$RGby1)$aic
    
    #Obtain the Variance partition coefficient
    #This is a measure of how much variance can be attributed to this particular clustering level
    VPC <- fit$var / full$var
    
    stats <- c(Accuracy, zscoreMean, DeviancePPMap, GAICPPMap, DevianceGbyS, GAICGbyS, DevianceSby1, GAICSby1, DevianceGby1, GAICGby1, DevianceRGby1, GAICRGby1, VPC)
    
    #Expand the result df with this data
    res <- cbind(res, stats)
  }
  
  #Apply labels for the results table
  colnames(res) <- cutoffs
  rownames(res) <- c("Accuracy", "zscoreMean", "DeviancePPMap", "GAICPPMap", "DevianceGbyS", "GAICGbyS", "DevianceSby1", "GAICSby1", "DevianceGby1", "GAICGby1", "DevianceRGby1", "GAICRGby1", "VPC")
  
  return(res)
}

#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
args = commandArgs(trailingOnly = T)
input <- read.csv(args[1], stringsAsFactors = F)

#Creates a graph based on the inputted data frame. The tn93 Distances become edge4 attributes
g <- graph_from_data_frame(input, directed=FALSE, vertices=NULL)

#Adds the ID's and Sample collection years as different vertex attributes for each vertex
temp <- sapply(V(g)$name, function(x) strsplit(x, '_')[[1]])
V(g)$name <- temp[1,]
V(g)$year <- as.numeric(temp[2,])

#Obtain the range of years
years <- levels(factor(V(g)$year))

#Test
##########################################################
#Define Resolution of plotted data
cutoffs <- seq(0, 0.05, 0.002)
res <- analyzeG(cutoffs)

#Create Plots
yLabs<- c("Accuracy", "zscoreMean", "DeviancePPMap", "GAICPPMap", "DevianceGbyS", "GAICGbyS", "DevianceSby1", "GAICSby1", "DevianceGby1", "GAICGby1", "DevianceRGby1", "GAICRGby1", "VPC")
for (i in 1:1){
  print(i)
  yLab <- yLabs[[i]]
  plot(cutoffs, res[i,], xlab= "tn93 Cutoffs", ylab= yLab)
}
##########################################################