#A process process for scoring the effectiveness of clustering from tn93 Output.
#USAGE: Rscript tn93Graph.R tn93output.csv

library(igraph)

#Takes in a graph and creates a train and test partition of that graph in the global environment
partition <- function(inG, trnRat) {
  #@param inG: A graph to be randomly split into 2 subgraphs
  #@param divide: Determines the ratio of the first subgraph (ie. trnRat=0.60 would mean 60% of the total graph would be the training subgraph )
  
  
  #Attain a list of indices that represent the training portion of the graph 
  partition <- floor(trnRat*length(V(inG)))
  index <- sample(length(V(inG)), size=partition)
  
  #Create a training partition based on the indices list and a test partition based on all indices not in the list
  train <<- induced_subgraph(inG, V(inG)[index], impl = "copy_and_delete")
  test <<- induced_subgraph(inG, V(inG)[-index], impl = "copy_and_delete")
}


#Forward Censors a graph, discluding all vertices past a certain year
forwardCensor <- function(inG, year) {
  #@param inG: A graph to be censored after a certain year
  #@param year: The year that defines future censorship
  #@return: The inputted graph excluding all vertices beyond year
  
  gUpToYear <- induced_subgraph(inG, V(inG)[V(inG)$years<=year], impl = "copy_and_delete")
  
  return(gUpToYear)
}

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
args = commandArgs(trailingOnly = T)
input <- read.csv(args[1], stringsAsFactors = F)

#Creates a graph based on the inputted data frame. The tn93 Distances become edge4 attributes
g <- graph_from_data_frame(input, directed=FALSE, vertices=NULL)

#Adds the ID's and Sample collection years as different vertex attributes for each vertex
temp <- sapply(V(g)$name, function(x) strsplit(x, '_')[[1]])
V(g)$name <- temp[1,]
V(g)$years <- as.numeric(temp[2,])