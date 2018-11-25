#A process process for scoring the effectiveness of clustering from tn93 Output.
#USAGE: Rscript tn93Graph.R tn93output.csv

library(igraph)
####- TO-DO: Use ML for optimization. see Caret or e107 packages -####

#__________________________________________________________________________________________________________________________#

#Estimates a cutoff distance. A subgraph based on edges below that cutoff should produce ideal clusters.
#A list of ideal clusters will have high variance across their growths
prettyGoodDist <- function(inG) {
  #@param inG: The input graph (usually a subgraph forward-censored by year). The latest year represents the Upcoming year. 
  #@return: A cutoff producing "ideal" clusters
  
  #This objective function defines the relationship between cutoff distance and cluster growth variance
  objF <- function(d) {
    cutG <- inG - E(inG)[E(inG)$Distance>=d]
    var <- var(getGrowth(cutG))
    return(var)
  }
  
  #Obtains the interface between vertices from the upcoming year and vertices from past years
  newV <- V(inG)[V(inG)$year==max(V(inG)$year)]
  presG <- V(inG)[V(inG)$year<max(V(inG)$year)]
  bridgeE <- E(inG)[newV%--%presG]
  
  #Initializes some data for modelling (variance~dist). The starting distance (x) is set to the shortest distance between a new vertex and an old vertex. 
  #The step size is based upon how short the starting distance is.
  variance <- c()
  dist <- c()
  x <- min(bridgeE$Distance)
  step <- 0.005
  
  #Generates the data for modelling. Stops generating data once the Maximum informative distance is reached. 
  #(Maximum informative distance is the point at which all old vertices cluster in 1 cluster)
  while (!is.na(objF(x))) {
    dist <- c(dist, x)
    variance <- c(variance, objF(x))
    x<-x+step
    print(x)
  }
  
  #Obtain the distance that maximized variance. 
  out <- dist[variance==max(variance)][[1]]
  
  #Create a poisson family glm based upon the recently generated data.
  ####- TO DO: This is currently unused -####
  distMod <- glm(variance~dist, family=poisson)
  
  return(out)
}

#Estimates the growth of clusters based on information from all years before the latest year.
####- TO DO: Add MANY more predictor variables for this estimation and improve upon the current ones -####
estimateGrowth <- function(inG) {
  #@param inG: The input graph (usually a subgraph forward-censored by year). The latest year represents the Upcoming year. 
  #@return: An attribute for the clusters representing ther predicted growth of each one.
  
  #Obtain the present graph (ie. the subgraph contain all cases before the the latest/"upcoming" year) and it's cluster info
  presV <- V(inG)[V(inG)$year<max(V(inG)$year)]
  presG <- induced_subgraph(inG, presV, impl = "auto")
  clu <- components(presG)
  
  #Obtain the size of current clusters after excluding all vertices from the current year. (ie. The clusters without their recent growth)
  ####- TO DO: Confirm this is a mathematically valid method when taking into account cluster merging -####
  oldClu <- clu$membership[attr(clu$membership, "names") %in% V(presG)[V(presG)$year<max(V(presG)$year)]$name]
  
  #Obtain the recent cluster growth for each cluster (ie. the difference between old cluster sizes and new ones)
  ####- TO DO: Potentially better style and speed with lapply?
  for (i in 1:clu$no) {
    presCsize <- clu$csize[[i]]
    oldCsize <- length(oldClu[unname(oldClu)==i])
    clu$growth[[i]] <- presCsize-oldCsize
  }
  
  return(clu$growth)
}

#Determines the growth of clusters based on the new addition of cases in the latest year (ie. The Upcoming year)
getGrowth <- function(inG) {
  #@param inG: The input graph (usually a subgraph forward-censored by year). The latest year represents the Upcoming year. 
  #@return: An attribute for the clusters representing the growth of each cluster
  
  newV <- V(inG)[V(inG)$year==max(V(inG)$year)]
  
  #Obtain the present graph (ie. the subgraph contain all cases before the the latest/"upcoming" year) and it's cluster info
  presG <- inG - newV
  clu <- components(presG)
  
  #Initialize cluster growth at 0 for each cluster
  clu$growth <- integer(clu$no)
  
  #Obtains the interface between vertices from the upcoming year and vertices from the present year.  
  bridgeE <- E(inG)[newV%--%V(inG)[-newV]]
  bridgeG <- inG - V(inG)[-ends(inG, bridgeE, names=F)]
  
  #Redifine the upcoming and present vertices as only those in the interface
  newV <- V(bridgeG)[V(bridgeG)$year==max(V(bridgeG)$year)]
  presV <- V(bridgeG)[V(bridgeG)$year<max(V(bridgeG)$year)]
  
  #Compare and obtain the growth of each cluster based on the placement of each new vertex
  #A vertex added to multiple clusters is considered to increase each cluster's size by a fraction
  if (length(newV)!=0){
    for (i in 1:length(newV)) {
      v <- newV[[i]]
      es <- E(bridgeG)[inc(v)]
      cluIndex <- unname(clu$membership[attr(clu$membership, "names") %in% ends(bridgeG, es, names = T)])
      vWeight <- 1/length(levels(factor(cluIndex))) 
      clu$growth[cluIndex] <- clu$growth[cluIndex]+vWeight 
    }
  }
  
  return(clu$growth)
}

getGrowthSimp <- function(inG) {
  
  inG <- subgraph.edges(inG, E(inG), delete.vertices = T)
  
  clu <- components(clu)
  newV <- V(inG)[V(inG)$year == max(V(ing)$year)]
  
  newC <- clu$membership[attr(clu$membership, "names") %in% newV$name]
  
  growthTable <- table(unname(newC))
  
  clu$growth <- growthTable[2,]
  
  return(clu)
}

#Obtains a filtered subgraph of the full graph.
subGraph <- function(inG, y, d) {
  #@param y: The year that represents the latest year. We forward censor everything past this.
  #@param d: The distance that represents the optimal cutoff distance for the whole graph
  #@return: The filtered graph (forward censored and cut by an optimal distance)
  
  #Removes vertices beyond a current year
  outV <- V(inG)[V(inG)$year>y]
  outG <- inG - outV
  
  #Removes edges with distances above a certain cutoff
  outE <- E(outG)[E(outG)$Distance>=d]
  outG <- outG - outE
  
  return(outG)
}

#__________________________________________________________________________________________________________________________#


#Warnings are currently generated by maximum informative distance being reached in the objective function. 
#This is handled by the functions that use objF
options(warn=-1)

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
years <- levels(factor(V(g)$year))

#Initialize a dataframe for the output of Accuracy by year
output <- data.frame(Year = as.integer(years), 
                     Accuracy = double(length(years)), 
                     Distance = double(length(years)), 
                     Time = double(length(years)))

#Obtain the optimum cutoff distance for that year
####- TO DO: Have this be based upon optDist function (glm-based estimate) and then erase the outdated function above -####
dist <- prettyGoodDist(g)

#Generate the Output data
for (y in years) {
  #If the case has no previous years to compare against
  if (length(V(g)[V(g)$year<y])==0) next
  
  startT <- proc.time() 
  
  #Filter the input graph to the loop year (y). The loop year will represent the "Upcoming" Year
  filtrG <- subGraph(g, y, dist)
  
  #Obtain the accuracy (ie. The actual - expected cluster growth)
  growthDiff <- estimateGrowth(filtrG) - abs(getGrowth(filtrG))
  acc <- mean(growthDiff)
  
  #Populate the output data with accuracy and ideal cutoff distance
  output$Accuracy[output$Year==y] <- acc
  output$Distance[output$Year==y] <- dist
  
  #Populate the output data with the speed of this iteration
  endT <- proc.time()
  output$Time[output$Year==y] <- (endT-startT)
}

#Test
print(proc.time())
cat("\n")
print(output)