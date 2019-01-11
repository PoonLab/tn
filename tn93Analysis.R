#A process process for scoring the effectiveness of clustering from tn93 Output.
#USAGE: Rscript tn93Graph.R tn93output.csv

library(igraph)
####- TO-DO: Use ML for optimization. see Caret or e107 packages -####

#__________________________________________________________________________________________________________________________#

#Estimates the growth of clusters based on information from all years before the latest year.
####- TO DO: Add MANY more predictor variables for this estimation and improve upon the current ones -####
estimateGrowth <- function(inG) {
  #@param inG: The input graph (usually a subgraph forward-censored by year). The latest year represents the Upcoming year. 
  #@return: An attribute for the clusters representing the predicted growth of each one.
  
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
  
  return(clu)
} ####- CURRENTLY UNUSED -####

#Determines the growth of clusters based on the new addition of cases in the latest year (ie. The Upcoming year)
getGrowth <- function(inG, alt=F, full=F) {
  #@param inG: The input graph (usually a subgraph forward-censored by year). The latest year represents the Upcoming year. 
  #@param alt: Checked true 
  #@return: An attribute for the clusters representing the growth of each cluster
  
  newV <- V(inG)[V(inG)$year==max(V(inG)$year)]
  
  #Obtain the present graph (ie. the subgraph contain all cases before the the latest/"upcoming" year) and it's cluster info
  presG <- inG - newV
  if (full) {
    presG <- presG - E(presG)
  }
  clu <- components(presG)
  
  #Initialize cluster growth at 0 for each cluster
  clu$growth <- integer(clu$no)
  
  #Obtains the interface between vertices from the upcoming year and vertices from the present year.  
  bridgeE <- E(inG)[newV%--%V(inG)[-newV]]
  if (length(bridgeE > 0)) {
    bridgeG <- subgraph.edges(inG, bridgeE, delete.vertices=T)
    
    #Create a new bipartite graph of future cases linked to previous (clustered) cases
    newV <- V(bridgeG)[V(bridgeG)$year==max(V(bridgeG)$year)]
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

#Obtains the Deviance and AIC of a comparison between our predicted growth as a measure of  
cutoffStats <- function(subG, full=F) {
  #Da = -2(la - lfull)

  clu <- getGrowth(subG, full=full)
  size <- clu$csize
  growth <- clu$growth
  
  # generate data frame at level of clusters=
  df <- data.frame(Size=size, Growth=growth)
  
  # number of new cases per cluster is count outcome 
  fit <- glm(Growth~Size, data=df, family='poisson')
  
  return(fit)
}

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

##########################################################

y <- max(years) 
cutoffs <- seq(0, 0.050, 0.002)
res <- {}

for (d in cutoffs) {
  print(d)
  subG <- subGraph(g, y, d)
  fit <- cutoffStats(subG) 
  
  deviance <- fit$deviance
  AIC <- fit$aic
  var <- var(fit$data$Growth)
  rel <- var(fitdata$Growth/sqrt(fit$data$Growth))
  no <- length(fit$data$Growth) 
  size <- mean(fit$data$Size)
  
  res <- cbind(res, c(Deviance, AIC, var, rel, no, size))
}

colnames(res) <- cutoffs
rownames(res) <- c('Deviance', 'GAIC', 'Variance in Absolute Growth', "Variance in Relative Growth", "Number of Clusters", "Mean Cluster Size" )

##########################################################