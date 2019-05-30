#Import Libraries
library(igraph,verbose = FALSE)
library(dplyr,verbose = FALSE)
library(parallel,verbose = FALSE)
library(ggplot2,verbose = FALSE)

#Expecting tn93 output as second param
## USAGE: Rscript ~/git/tn/robRe.R ~/Seattle/tn93St.txt ~/Seattle/Robust/RRout/ ##
## args <- c("~/Seattle/tn93St.txt", "~/Seattle/Robust/RRout/") ##

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
    cE <- lapply(V(iG)[year==nY], function(x) {
      xE <- bE[inc(x)]
      
      #To catch a case that is new, but has no linkages to old cases
      ifelse(length(xE)==0, 0, (xE[Distance == min(Distance)])[1] ) 
    })#, mc.cores=8)
    
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

optFind <- function(iG) {
  
  years <- as.integer(levels(factor(V(iG)$year)))
  nY <- max(years)
  V(iG)$tDiff <- sapply(V(iG)$year, function(x) nY-x)
  
  #Initialize a set of cutoffs to observe
  steps <- head(hist(E(iG)$Distance, plot=FALSE)$breaks,-5)
  cutoffs <- seq(0 , max(steps), max(steps)/50)
  
  iG <- minFilt(iG)
  
  #Create a set of subgraphs at each cutoff
  gs <- lapply(cutoffs, function(d) {
    subgraph.edges(iG,E(iG)[Distance<=d], delete.vertices = F)
  })#, mc.cores=8) 
  names(gs) <- cutoffs
  
  res <- lapply(cutoffs, function(d) {
    cat(paste0("\r", "Running Analysis ", d/max(cutoffs)*100, "%"))
    #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
    subG <- gs[[as.character(d)]]
    
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
  })#, mc.cores=8)
  
  #Label data
  names(res) <- cutoffs
  
  
  return(res)
}

## Importing Case data
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
args = commandArgs(trailingOnly = T)
tn93out <- read.csv(args[1], stringsAsFactors = F)

#This script will give warnings due to the fact that there are low fit rates on the null model
options(warn=-1)

#Creates a graph based on the inputted data frame. The tn93 Distances become edge4 attributes
g <- graph_from_data_frame(tn93out, directed=F, vertices=NULL)

#Adds the ID's and Sample collection years as different vertex attributes for each vertex
temp <- sapply(V(g)$name, function(x) strsplit(x, '_')[[1]])
V(g)$name <- temp[1,]
V(g)$year <- as.numeric(temp[2,])

#Obtain the range of years and the maximum input year
years <- as.integer(levels(factor(V(g)$year)))
nY <- max(years)
while (length(V(g)[year==nY])<63) {nY <- nY-1}

g <- induced_subgraph(g, V(g)[year<=nY])

for (i in 1:10){
  fname <- paste0(args[2],"St8s",i,".rds")
  res1 <- induced_subgraph(g,sample(V(g), round(length(V(g))*0.80)))
  saveRDS(optFind(res1), file = fname)
  
  fname <- paste0(args[2],"St6s",i,".rds")
  res2 <- induced_subgraph(g,sample(V(g), round(length(V(g))*0.60)))
  saveRDS(optFind(res2), file = fname)

  fname <- paste0(args[2],"St4s",i,".rds")
  res3 <- induced_subgraph(g,sample(V(g), round(length(V(g))*0.40)))
  saveRDS(optFind(res3), file = fname)
}
