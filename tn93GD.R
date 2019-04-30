#A process which generates cluster growth data as a function of tn93 cutoff threshold
#Creates an external .RData file of paired cluster info sets

### USAGE: Rscript ~/git/tn/tn93GD.R tn93__.txt ###
#input Date Format, specified with % (ie. %d-%b-%y for day, written month, 2-digit year or  %Y for simple, 4-digit year)

#Import Libraries
library(igraph)
library(dplyr)
library(MASS)
require(parallel)

## Helper Functions
#____________________________________________________________________________________________________________________________#

#Obtains a data set of all possible Bipartite Edge Frequencies in a given subgraph, with the bipartite subsets of vertices representing a given year of data
bpeFreq <- function(iG) {
  #@param ig: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A data frame of Number of positives (edges from one year to the newest year) with total possible edges and time difference (in years) between the two years
  
  #Obtain the range of years
  maxY <- max(V(iG)$year)
  minY <- min(V(iG)$year)
  ys <- seq(minY, (maxY-1), 1)
  nV <- V(iG)[year==maxY]
  
  #Obtain the frequency of new cases being connected to each year
  frequency <- sapply(ys, function(x) {
    pV <- V(iG)[year==x]
    bE <- E(iG)[pV%--%nV]
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
minFilt <- function(iG) {
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A filtered version of this same graph
  
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
      ifelse(length(xE)==0, 0, (xE[Distance == min(Distance)])[1] ) 
    }, mc.cores=8)
    
    #Remove the entries from new cases that dont connect to  old cases
    cE <- unname(unlist(cE[cE!=0]))

    #Filter out all edges except for the closest edges
    if(!is.null(cE)){
      iG <- iG - difference(bE, E(iG)[cE])
    }
  }
  
  return(iG)
}

#Simulates the growth of clusters
simGrow <- function(iG) { 
  
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
nY <- max(years)
V(g)$tDiff <- sapply(V(g)$year, function(x) nY-x)

#Initialize a set of cutoffs to observe
cutoffs <- seq(0.005, 0.05, 0.001)

g <- minFilt(g)

#Create a set of subgraphs at each cutoff
gs <- mclapply(cutoffs, function(d) {
  subgraph.edges(g,E(g)[Distance<=d], delete.vertices = F)
}, mc.cores=8) 
names(gs) <- cutoffs

## Obtain a set models of case linkage frequency based on age
ageD <- mclapply(gs, function(x){
  l <- mclapply(rev(tail(years, -1)), function(y){
    print(y)
    subG <- minFilt(induced_subgraph(x, V(x)[year<=y]))
    bpeFreq(subG)
  } , mc.cores=8)
  bind_rows(l) 
}, mc.cores=8) 

#Save data in accessable file
saveRDS(ageD, file = paste0(gsub("\\..*", "", args), "AD.rds"))

## Generate Growth data
#__________________________________________________________________________________________________________________________#

res <- mclapply(cutoffs, function(d) {
  #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
  subG <- gs[[as.character(d)]]
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case ag
  #This data may contain missing cases, hense the complete cases addition
  ageDi <- ageD[[as.character(d)]]
  
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
  
  #Calculate GAIC
  clu$gaic <- mod1$aic-mod2$aic
  
  return(clu)
}, mc.cores=8)

#Label data
names(res) <- cutoffs

#Save data in accessable files
saveRDS(res, file = paste0(gsub("\\..*", "", args), "GD.rds"))