#A process process for scoring the effectiveness of clustering from tn93 Output.
#USAGE: Rscript tn93Graph.R tn93output.csv

library(igraph)
####- TO-DO: Use ML for optimization. see Caret or e107 packages -####

#__________________________________________________________________________________________________________________________#

altGrowth <- function(inG, y, full=F) {
  clu <- components(inG)
  oldV <- V(inG)[V(inG)$year<=y]
  oldClu <- clu$membership[attr(clu$membership, "names") %in% oldV$name]
  
  if (full) {
    clu$growth <- integer(length(oldV)) 
    clu$growth <- sapply(oldV, function(x) length(E(inG)[inc(x)&&inc(E(inG)[-oldV])]))
  } else {
    clu$growth <- sapply(1:clu$no, function(x) clu$csize[[x]]-length(oldClu[unname(oldClu)==x]))
  }
  
  return(clu)
}

#Estimates the growth of clusters based on information from all years before the latest year.
####- TO DO: Decide upon a statistically reasonable way to include the "Full" parameter -####
####- TO DO: Implement use of this function -####
####- TO DO: Consider LASSO functionality to add more predictor variables for this estimation -####
estimateGrowth <- function(inG, full=F) {
  #@param inG: The input graph (usually a subgraph forward-censored by year). The latest year represents the Upcoming year.
  #@param full: A boolean determining whether or not this is the growth estimate for a fully saturated model
  #@return: An attribute for the clusters representing the predicted growth of each one.
  
  #Obtain the present graph (ie. the subgraph contain all cases before the the latest/"upcoming" year) and it's cluster info
  presV <- V(inG)[V(inG)$year<max(V(inG)$year)]
  presG <- induced_subgraph(inG, presV, impl = "auto")
  clu <- components(presG)
  
  #Obtain the size of current clusters after excluding all vertices from the current year. (ie. The clusters without their recent growth)
  ####- TO DO: Confirm this is a mathematically valid method when taking into account cluster merging -####
  oldClu <- clu$membership[attr(clu$membership, "names") %in% V(presG)[V(presG)$year<max(V(presG)$year)]$name]
  
  #Obtain the recent cluster growth for each cluster (ie. the difference between old cluster sizes and new ones)
  ####- TO DO: Potentially better style and speed with lapply? -####
  for (i in 1:clu$no) {
    presCsize <- clu$csize[[i]]
    oldCsize <- length(oldClu[unname(oldClu)==i])
    clu$growth[[i]] <- presCsize-oldCsize
  }
  
  return(clu)
} 

#Determines the growth of clusters based on the new addition of cases in the latest year (ie. The Upcoming year)
####- TO DO: Test and patch up the alt parameter so that the alternative method is usable -####
getGrowth <- function(inG, alt=F, full=F) {
  #@param inG: The input graph (usually a subgraph forward-censored by year). The latest year represents the Upcoming year. 
  #@param alt: A boolean determining whether or not to resolve cluster growth as weighted nodes or closest nodes
  #@param full: A boolean determining whether or not this is the growth estimate for a fully saturated model
  #@return: An attribute for the clusters representing the growth of each cluster
  
  #Defines the new vertices (which determine cluster growth)
  newV <- V(inG)[V(inG)$year==max(V(inG)$year)]
  
  #Obtain the present graph (ie. the subgraph contain all cases before the the latest/"upcoming" year) and it's cluster info
  presG <- inG - newV
  if (full) {
    presG <- presG - E(presG)
  }
  clu <- components(presG)
  
  #Initialize cluster growth at 0 for each cluster
  clu$growth <- integer(clu$no)
  clu$inc <-0
  
  #Obtains the interface between vertices from the upcoming year and vertices from the present year.  
  bridgeE <- E(inG)[newV%--%V(inG)[-newV]]
  if (length(bridgeE > 0)) {
    bridgeG <- subgraph.edges(inG, bridgeE, delete.vertices=T)
    
    #Create a new bipartite graph of future cases linked to previous (clustered) cases
    newV <- V(bridgeG)[V(bridgeG)$year==max(V(bridgeG)$year)]
    clu$inc <- length(newV) /  length(V(presG))
    presV <- V(bridgeG)[-newV]
    
    #### TO-DO: Finish alt-option (currently broken, this option would assign case addition based on closest cluster) -####
    if(alt){
      #Get edge id's of all of the shortest edge lengths (A case can only be linked to 1 case)
      closeE <- lapply(newV, function(x) {
        xE <- E(bridgeG)[inc(x)]
        closest <- xE[xE$Distance == min(xE$Distance)]
        return (closest[1])
      })
      
      bridgeG <- subgraph.edges(bridgeG, closeE, delete.vertices = T)
      weight <- rep(newV, function(x) 1)
      
    } else {
      #Establish the weight of a future case based on the number of past cases it's added to
      weight <- sapply(newV, function(x) 1/length(E(bridgeG)[inc(x)])) 
    }
    
    #Create a table of every link from past to future cases and the weight that those links are worth
    temp <- unname(sapply(presV, function(x){
      point <- sum(unname(weight[neighbors(bridgeG, x, "all")$name])) 
      cluId <- unname(clu$membership[V(bridgeG)[x]$name])
      return (c(cluId, point))
    }))
    
    #Assign growth as an attribute of the set of clusters based off of the weighted information from temp 
    clu$growth <- sapply(seq(1,clu$no), function(x) round(sum(temp[2,][temp[1,]==x])))
  }
  
  return(clu)
}

#Obtains a filtered subgraph of the full graph.
subGraph <- function(inG, y, d, plot=F) {
  #@param y: The year that represents the latest year. We forward censor everything past this.
  #@param d: The distance that represents the optimal cutoff distance for the whole graph
  #@param plot: Creates an easy to view subgraph for the purposes of plotting and overview, but not for statistical analysis.
  #@return: The filtered graph (forward censored and cut by an optimal distance)
  
  #Removes vertices beyond a current year
  outV <- V(inG)[V(inG)$year>y]
  outG <- inG - outV
  
  #Removes edges with distances above a certain cutoff
  outE <- E(outG)[E(outG)$Distance>=d]
  outG <- outG - outE
  if(plot){ #Ignores clusters of size 1 and provides a graph (for ease of overview and no calculations)
    outG <- subgraph.edges(outG, E(outG), delete.vertices = T)
    plot(outG, vertex.size = 2, vertex.label = NA, edge.width = 0.65, edge.color = 'black')
    print(components(outG))
  }   
  else {return(outG)}
}

#Obtains fit measurements for 
stats <- function(subG, full=F) {
  #@param inG: The input graph (usually a subgraph forward-censored by year). The latest year represents the Upcoming year. 
  #@param full: A boolean determining whether or not this is the growth estimate for a fully saturated model
  #@return: The glm, Poisson linked tracking growth relative to size as well as a poisson probability map

  #Obtain the cluster growth and size based upon most recent cases and cutoff parameter for subG
  clu <- getGrowth(subG,full=full)
  exp <- clu$inc*(clu$csize)
  relGrowth <- clu$growth - exp
  zscore <- relGrowth / sqrt(exp)

  #Generate data frame at level of clusters
  df <- data.frame(Expectation=exp, Actual=clu$growth)
  
  #Number of new cases per cluster is count outcome
  fit <- glm(Actual~Expectation, data=df, family = "poisson")
  fit$var <- var(clu$growth) 
  
  return(fit)
}

altStats <- function(subG, full=F)
  

#__________________________________________________________________________________________________________________________#

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

#Plot generation, currently records deviance and GAIC for a MAUP based model as modulated by Cutoff
##########################################################

y <- max(years) 
cutoffs <- seq(0, 0.05, 0.002)
res <- {}

for (d in cutoffs) {
  print(d)
  subG <- subGraph(g,y,d)
  fit <- stats(subG)
  full <- stats(subG, full=T)
  
  Deviance <- full$deviance - fit$deviance
  GAIC <- fit$aic- full$aic
  VPC <- fit$var / full$var
  
  res <- cbind(res, c(Deviance, GAIC, VPC))
}

colnames(res) <- cutoffs
rownames(res) <- c("Deviance", "GAIC", "VPC")

plot(cutoffs, res[1,], xlab= "tn93 Cutoffs", ylab= "Deviance" )
plot(cutoffs, res[2,], xlab= "tn93 Cutoffs", ylab= "GAIC" )
plot(cutoffs, res[3,], xlab= "tn93 Cutoffs", ylab= "VPC")

##########################################################
