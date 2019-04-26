#A process which generates cluster growth data as a function of tn93 cutoff threshold
#Creates an external .RData file of paired cluster info sets

### USAGE: Rscript ~/git/tn/tn93GD.R tn93__.txt dateFormat ###
#input Date Format, specified with % (ie. %d-%b-%y for day, written month, 2-digit year or  %Y for simple, 4-digit year)

#Import Libraries
library(igraph)
library(dplyr)
library(MASS)

require(parallel)

## Helper Functions
#____________________________________________________________________________________________________________________________#

#Models frequency of new cases being linked to old cases based on how old the old cases are
linkFreq <- function(inG) {
  #@param inG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A model new case attachment frequency as a function of case age
  
  #Obtain the range of years
  maxY <- max(V(inG)$year)
  minY <- min(V(inG)$year)
  years <- seq(minY, (maxY-1), 1)
  newV <- V(inG)[year==maxY]
  
  #Obtain the frequency of new cases being connected to each year
  frequency <- sapply(years, function(x) {
    presV <- V(inG)[year==x]
    bridgeE <- E(inG)[presV%--%newV]
    pos <- length(bridgeE)
    tot <- length(presV)*length(newV)
    return(c(pos,tot))
  })
  
  #Assign age to every case
  age <- sapply(years, function(x) maxY-x)

  #Create a data frame of case attachment frequency and case age
  df <- data.frame(Age = age, Positive = frequency[1,], Total = frequency[2,])

  return(df)
}

#Obtains a filtered subgraph of the full graph. Vertices are removed beyond a given year and edges are removed below a cutoff
subGraph <- function(inG, y, d) {
  #@param y: The year that represents the latest year. We forward-censor everything past this.
  #@param d: The distance that represents the cutoff threshold. We remove all edges above this.
  #@return: The filtered graph (forward censored and cut by a given distance)
  
  #Removes vertices beyond a current year
  outV <- V(inG)[V(inG)$year>y]
  outG <- inG - outV
  
  #Removes edges with distances above a certain cutoff
  outE <- E(outG)[E(outG)$Distance>=d]
  outG <- outG - outE
  
  return(outG)
}

#Filters the input graph such that all new cases are only linked to old cases by their closest edge to old cases
closeFilter <- function(inG) {
  #@param inG: A subG gaph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A filtered version of this same graph
  
  #Obtain the maximum year
  y <- max(V(inG)$year)
  
  #Obtain edge id's of all of the shortest edge lengths from new cases (A new case can only be linked to 1 case)
  bridgeE <- E(inG)[V(inG)[year==y]%--%V(inG)[year<y]]
  
  #To catch a case where no new cases link to old ones
  if (length(bridgeE) > 0) {
    
    #Obtain the closest edges for each new case
    closeE <- unname(sapply(V(inG)[year==y], function(x) {
      xE <- bridgeE[inc(x)]
      
      #To catch a case that is new, but has no linkages to old cases
      if(length(xE)==0) {
        return(NULL)
      }
      else {    
        closest <- xE[Distance == min(Distance)]
        return (closest[1])
      }
    }))
    closeE <- unlist(closeE[!vapply(closeE, is.null, logical(1))])
    farE <- difference(bridgeE, E(inG)[closeE])
    
    #Filter out all edges except for the closest edges
    if(!is.null(closeE)){
      outG <- inG - farE
    }
  }
  
  #To return a value in the case that no new cases link to old ones
  else {
    outG <- inG
  }
  
  return(outG)
}

#Simulates the growth of clusters
simGrow <- function(inG, full=F) { 
  
  #Obtain the newest date
  newV <- V(inG)[year==newY]

  #Obtain a forward-Censored Graph
  oldG <- inG - newV
  
  if (full) {
    oldG <- oldG-E(oldG)
  }
  
  #Obtain cluster information
  clu <- components(oldG)

  #Assign cluster growth based on number of new cases linked to old cases in clusters 
  temp <- sapply(1:clu$no, function(x) {
    members <- names(clu$membership[unname(clu$membership)==x])
    memV <- V(inG)[name%in%members]
    bridgeE <- E(inG)[memV%--%newV]
    forecast <- sum(memV$freq) 
    growth <- length(bridgeE)
    return(c(growth,forecast))
  })
  
  #Assign growth, forecast (based on diagnostic date), and incidence
  clu$growth <- temp[1,]
  clu$forecast <- temp[2,]
  clu$inc <- length(newV)
  
  return(clu)
}

## Importing Case data
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
args = commandArgs(trailingOnly = T)
input <- read.csv(args[1], stringsAsFactors = F)

#This script will give warnings due to the fact that there are low fit rates on the null model
options(warn=-1)

#Creates a graph based on the inputted data frame. The tn93 Distances become edge4 attributes
g <- graph_from_data_frame(input, directed=F, vertices=NULL)

#Adds the ID's and Sample collection years as different vertex attributes for each vertex
temp <- sapply(V(g)$name, function(x) strsplit(x, '_')[[1]])
V(g)$name <- temp[1,]
V(g)$year <- as.numeric(temp[2,])

#Obtain the range of years and the maximum input year
years <- as.integer(levels(factor(V(g)$year)))
newY <- max(years)
V(g)$age <- sapply(V(g)$year, function(x) newY-x)

#Initialize a set of cutoffs to observe
cutoffs <- seq(0.005, 0.05, 0.001)

#Create a set of subgraphs at each cutoff
gs <- mclapply(cutoffs, function(d) {
  print(d)
  subGraph(g,max(years),d)
}, mc.cores=8) 
names(gs) <- cutoffs

## Obtain a set models of case linkage frequency based on age
#__________________________________________________________________________________________________________________________#

##TO-DO: Reformat this output, it is significantly faster
ageD <- mclapply(gs, function(graph){
  bind_rows(mclapply(rev(tail(years,-2)), function(y){
    subG <- closeFilter(induced_subgraph(graph, V(graph)[year<=y]))
    linkFreq(subG)
  }, mc.cores = 8))
}, mc.cores=8) 

#Save data in accessable file
saveRDS(ageD, file = paste0(gsub("\\..*", "", args), "AD.rds"))


## Generate Growth data
#__________________________________________________________________________________________________________________________#

#Obtain a filtered graph and measure growth
g <- closeFilter(g)

res <- mclapply(cutoffs, function(d) {
  #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
  subG <- gs[[as.character(d)]]
  ageDi <- ageD[[as.character(d)]]
  
  mod <- glm(cbind(Positive, Total) ~ Age, data=ageDi, family='binomial')
  
  #Assign a predicted growth value to each member of the graph
  V(subG)$freq <- predict(mod, data.frame(Age=V(subG)$age), type='response')
  
  #Obtain growth based on two models restricted model
  clu <- simGrow(subG)
  
  #Place growth and forecast data in dfs for fit and full growth
  df1 <- data.frame(Growth = clu$growth, Pred = clu$forecast)
  df2 <- data.frame(Growth = clu$growth, Pred = clu$csize*(sum(clu$growth)/sum(clu$csize)))
  
  #Model growth as a function of forecast for fit and full growth models
  mod1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
  mod2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
  
  #Calculate GAIC
  clu$gaic <- mod1$aic-mod2$aic
  
  return(clu)
}, mc.cores=8)

#Label data
names(res) <- cutoffs

#Save data in accessable file
saveRDS(res, file = paste0(gsub("\\..*", "", args), "GD.rds"))