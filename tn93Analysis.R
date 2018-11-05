#A process process for scoring the effectiveness of clustering from tn93 Output.
#USAGE: Rscript tn93Graph.R tn93output.csv

library(igraph)
####- TO-DO: Use ML for optimization. see Caret or e107 packages -####

#optimize <- function(inG, y) {}

#_______________________________________________________________________________________________________________#

#Generates a difference representing how accurately a set of predictor variablespredict cluster size. 
#This score is based off of the difference between clusters actual growth in a year vs their estimate of their own growth.
#Currently the only predictor of growth is recent cluster growth.
predict <- function(inG, y) {
  #@param inG: The graph representing all of the data
  #@param y: The reference point year (recent growth and clusters will be based off of this year)
  #@param dist: The cutoff genetic distance. Only vertices below this distance may be in the same cluster
  #@return: The Average difference between estimated and actual growth for clusters.
  
  #Sub Graph representing the total input graph sensored up to the next year
  vertUpToYear <- V(inG)[V(inG)$year<=(y+1)]
  newG <- induced_subgraph(inG, vertUpToYear , impl = "copy_and_delete")
  
  #Optimization function for distance
  f <- function(d) {
    #Obtain a subgraph of all edges below a given distance
    eWithinDist <- E(newG)[E(newG)$Distance<=d]
    cutG <- subgraph.edges(newG, eWithinDist, delete.vertices = F)
    clu <- components(cutG)

    #obtain a subgraph of only the interface between the present and next year
    bridgeEs <- E(cutG)[(V(cutG)[V(cutG)$year==(y+1)]) %--% (V(cutG)[V(cutG)$year<=y])]
    bridgeG <- subgraph.edges(cutG, bridgeEs, delete.vertices=T)
    
    #The percentage of current clusters that are capturing new cases
    presClu <- length(levels(factor(unname(clu$membership[attr(clu$membership, "names") %in% V(newG)[V(newG)$year <= y]$name]))))/clu$no
    #The number of new cases captured, based on the distance
    newCases <- length(V(bridgeG)[V(bridgeG)$year==(y+1)]) 
    
    ratV <- newCases/presClu
    return(ratV)
  }
  
  dist <- optimize(f, lower = 0, upper = 1, maximum = T, tol = 0.001)
  
  #Filter forward censored subgraph to optimal distance
  eWithinDist <- E(newG)[E(newG)$Distance<=d]
  newG <- subgraph.edges(newG, eWithinDist, delete.vertices = F)
  
  #obtains a subgraph, forward censored only up to the current year
  vertUpToYear <- V(inG)[V(inG)$year<=y]
  presG <- induced_subgraph(newG, vertUpToYear, delete.vertices = F)
  
  ####- TO-DO: This could be its own function (ie. everything below could be predict, above is just optimize) -####
#___________________________________________________________________________________________________________________#
  
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
  E(newG)[(V(newG)[V(newG)$year==(y+1)]) %--% (V(newG)[V(newG)$year==y])]
  
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
numYears <- length(levels(factor(V(g)$year)))

df <- data.frame(Year = integer(numYears), Score = double(numYears))

temp <- sapply(df$year, function(x) predict(g, x))
df$Score <- temp
print(df)