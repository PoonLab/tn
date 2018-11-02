#A process process for scoring the effectiveness of clustering from tn93 Output.
#USAGE: Rscript tn93Graph.R tn93output.csv

library(igraph)
####- TO-DO: Use ML for optimization. see Caret or e107 packages -####

#Takes in a graph and creates a train and test partition of that graph in the global environment
partition <- function(inG, trnRat) {
  #@param inG: A graph to be randomly split into 2 subgraphs
  #@param divide: Determines the ratio of the first subgraph (ie. trnRat=0.60 would mean 60% of the total graph would be the training subgraph )
  
  
  #Attain a list of indices that represent the training portion of the graph 
  partition <- floor(trnRat*length(V(inG)))
  index <- sample(length(V(inG)), size=partition)
  
  #Create a training partition based on the indices list and a test partition based on all indices not in the list
  train <- induced_subgraph(inG, V(inG)[index], impl = "copy_and_delete")
  test <- induced_subgraph(inG, V(inG)[-index], impl = "copy_and_delete")
  
  output <- list(2)
  output[0] <- train
  output[1] <- test
  
  return(output)
}

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

#In Progress...
predict <- function(inG, y, dist) {
  #@param inG:
  #@param y:
  #@param dist:
  #@return:
  
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
  #!!!BROKEN: Possibly creating empty graphs.
  bridgeEs <- E(newG)[(V(newG)[V(newG)$years==(year+1)]) %--% (V(newG)[V(newG)$years==year])]
  bridgeG <- subgraph.edges(newG, bridgeEs, delete.vertices = T) 
  newVs <- V(bridgeG)[V(bridgeG)$year==(y+1)]
  
  #!!!BROKEN: Recieving 0 input crashes
  #!!!BROKEN: Does not award full score to a cluster that recieves the same vertex twice
  for (i in 1:length(newVs)) {
    print(i)
    v <- newVs[[i]]
    print(v)
    es <- E(bridgeG)[inc(v)]
    vWeight <- 1/length(es)
    print(vWeight)
    cluIndex <- unname(clu$membership[attr(clu$membership, "names") %in% ends(bridgeG, es, names = T)])
    print(cluIndex)
    temp <- sapply(clu$growth, function(x) x+vWeight) 
    clu$growth <- temp
  }
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