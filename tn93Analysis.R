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
    cutG <- subgraph.edges(inG, E(inG)[E(inG)$Distance<=d], delete.vertices=F)
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
  x <- min(E(inG)[bridgeE]$Distance)
  step <- 10^(floor(log10(x)))
  
  #Generates the data for modelling. Stops generating data once the Maximum informative distance is reached. 
  #(Maximum informative distance is the point at which all old vertices cluster in 1 cluster)
  while (!is.na(objF(x))) {
    dist <- c(dist, x)
    variance <- c(variance, objF(x))
    x<-x+step
  }
  
  #Obtain the distance that maximized variance. 
  out <- dist[variance==max(variance)][[1]]

  return(out)
}

#Obtain the Generalized Linear Model of the relationship between cutoff distance (at which edges appear between clusters) and variation in cluster growth.
optDist <- function(inG) {
  #@param inG: The input graph (usually a subgraph forward-censored by year). The latest year represents the Upcoming year. 
  #@return: The Poisson Family glm, representing the relationship between cutoff distance and cluster growth variance
  
  #This objective function defines the relationship between cutoff distance and cluster growth variance
  objF <- function(d) {
    cutG <- subgraph.edges(inG, E(inG)[E(inG)$Distance<=d], delete.vertices=F)
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
  x <- min(E(inG)[bridgeE]$Distance)
  step <- 10^(floor(log10(x)))
  
  #Generates the data for modelling. Stops generating data once the Maximum informative distance is reached. 
  #(Maximum informative distance is the point at which all old vertices cluster in 1 cluster)
  while (!is.na(objF(x))) {
    dist <- c(dist, x)
    variance <- c(variance, objF(x))
    x<-x+step
  }

  #Create a poisson family glm based upon the recently generated data.
  distMod <- glm(scores~dist, family=poisson)
  
  return(distMod)
}

#Estimates the growth of clusters based on information from all years before the latest year.
####- TO DO: Add MANY more predictor variables for this estimation and improve upon the current ones -####
estimateGrowth <- function(inG) {
  #@param inG: The input graph (usually a subgraph forward-censored by year). The latest year represents the Upcoming year. 
  #@return: An estimate of cluster growth based on information before the upcoming year
  
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
  #@return:
  
  #Represents the upcoming vertices to be added to the present cluster
  newV <- V(inG)[V(inG)$year==max(V(inG)$year)]
  
  #Obtain the present graph (ie. the subgraph contain all cases before the the latest/"upcoming" year) and it's cluster info
  presV <- V(inG)[V(inG)$year<max(V(inG)$year)]
  presG <- induced_subgraph(inG, presV, impl = "auto")
  clu <- components(presG)
  
  #Initialize cluster growth at 0 for each cluster
  clu$growth <- integer(clu$no)
  
  #Obtains the interface between vertices from the upcoming year and vertices from the present year.  
  bridgeE <- E(inG)[newV%--%presV]
  bridgeG <- subgraph.edges(inG, bridgeE, delete.vertices=T)
  
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

#Generate the Output data
for (y in years) {
  #If the case has no previous years to compare against
  if (length(V(g)[V(g)$year<y])==0) next
  
  startT <- proc.time() 
  
  #Filter the input graph to the loop year (y). The loop year will represent the "Upcoming" Year
  filtrV <- V(g)[V(g)$year<=y]
  filtrG <- induced_subgraph(g, filtrV, impl="auto")
  
  #Obtain the optimum cutoff distance for that year
  ####- TO DO: Have this be based upon optDist function (glm-based estimate) and then erase the outdated function above -####
  dist <- prettyGoodDist(filtrG)
  
  #If the case has no new vertices that will be added before the maximum informative distance
  if (is.null(dist)) next
  
  #Removes edges based on the optimum cutoff for cluster growth variance
  filtrE <- E(filtrG)[E(filtrG)$Distance<=dist]
  cutG <- subgraph.edges(filtrG, filtrE, delete.vertices=F)
  
  #Obtain the accuracy (ie. The actual - expected cluster growth)
  growthDiff <- estimateGrowth(cutG) - abs(getGrowth(cutG))
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