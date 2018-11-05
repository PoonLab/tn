#A process process for scoring the effectiveness of clustering from tn93 Output.
#USAGE: Rscript tn93Graph.R tn93output.csv

library(igraph)
####- TO-DO: Use ML for optimization. see Caret or e107 packages -####

#Takes in a graph and returns a subgraph containing only vertices up to a certain year and edges under a certain distance
cutAndCensor <- function(inG, y, dist) {
  #@param inG: The graph to be processed
  #@param year: The specified year, beyond which all vertices will be censored
  #@param dist: The specified distance, beyond which all edges will be excluded
  #@return The processed subgraph
  
  vertUpToYear <- V(inG)[V(inG)$year<=y]
  outG <- induced_subgraph(inG, vertUpToYear , impl = "copy_and_delete")
  
  eWithinDist <- E(outG)[E(outG)$Distance<=dist]
  outG <- subgraph.edges(outG, eWithinDist, delete.vertices = FALSE)
  
  return(outG)
}

#Generates a difference representing how accurately a set of predictor variablespredict cluster size. 
#This score is based off of the difference between clusters actual growth in a year vs their estimate of their own growth.
#Currently the only predictor of growth is recent cluster growth.
predict <- function(inG, y, dist) {
  #@param inG: The graph representing all of the data
  #@param y: The reference point year (recent growth and clusters will be based off of this year)
  #@param dist: The cutoff genetic distance. Only vertices below this distance may be in the same cluster
  #@return: The Average difference between estimated and actual growth for clusters.
  
  #Sub Graphs representing the total input graph sensored up to the next year and present year (respectively)
  newG <- cutAndCensor(inG, (y+1), dist)
  presG <- cutAndCensor(newG, y, dist)
  
  #Obtaining clusters. In this case, simply connected components within dist parameter
  clu <- components(presG)
  
  #A subset of clusters and their membership from a previous year
  ####- TO-DO: Review the validity of this retrospective method. The subset of members in previous years may not be connected -####
  oldClu <- clu$membership[attr(clu$membership, "names") %in% V(newG)[V(newG)$year <= (y-1)]$name]
  
  #Check the current sizes of clusters against their previous growth to obtain recent growth predictor variable
  ####- TO-DO: Possibly better style with sapply? ie. temp <- sapply(clu$growth, function(x) x/length(oldClu[unname(oldClu)==?CURRENTINDEX?])) -####
  for (i in 1:clu$no) {
    presCsize <- clu$csize[[i]]
    oldCsize <- length(oldClu[unname(oldClu)==i])
    clu$growth[[i]] <- presCsize-oldCsize
  }
  
  #Creates a sub graph, representing only the interface between the current and the next year
  #Only edges between a current year and next year vertex are included
  bridgeEs <- E(newG)[(V(newG)[V(newG)$year==(y+1)]) %--% (V(newG)[V(newG)$year==y])]
  bridgeG <- subgraph.edges(newG, bridgeEs, delete.vertices = T) 
  print(bridgeG)
  newVs <- V(bridgeG)[V(bridgeG)$year==(y+1)]
  
  #Based on the bridge between the present and new year, sees the actual growth of clusters and compares the
  if (length(newVs)!=0){
    for (i in 1:length(newVs)) {
      v <- newVs[[i]]
      es <- E(bridgeG)[inc(v)]
      vWeight <- 1/length(es)
      cluIndex <- unname(clu$membership[attr(clu$membership, "names") %in% ends(bridgeG, es, names = T)])
      temp <- sapply(clu$growth, function(x) x-vWeight) 
      clu$growth <- temp
    }
  }
  
  diff <- mean(sapply(clu$growth, function(x) abs(x)))
  
  return(diff)
}

#_______________________________________________________________________________________________________________#

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
years <- as.numeric(levels(factor(V(g)$year)))

###-TO DO: Optimize distance cutoff input to obtain reasonable clusters-####

#Test
for(i in 5:5) {
  predict(g, (2000+i), 0.05)
}