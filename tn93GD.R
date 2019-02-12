#A process which generates cluster growth data as a function of tn93 cutoff threshold
#Creates an external .RData file of paired cluster info sets

### USAGE: Rscript tn93GD.R tn93output.csv ###

library(igraph)

#____________________________________________________________________________________________________________________________#

#Models frequency of new cases being linked to old cases based on how old the old cases are
linkFreq <- function(inG) {
  #@param inG: A subG gaph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A model new case attachment frequency as a function of case age
  
  #Obtain the range of years
  maxY <- max(V(inG)$year)
  minY <- min(V(inG)$year)
  years <- seq(minY, (maxY-1), 1)
  newV <- V(inG)[year==maxY]
  
  #Obtain the frequency of new cases being connected to each year
  frequency <- sapply(years, function(x) {
    presV <- V(inG)[year==x]
    bridgeV <- ends(inG, E(inG)[presV%--%newV])
    freq <- length(bridgeV[bridgeV%in%newV$name])
    return(freq)
  })
  
  #Assign age to every case
  age <- sapply(years, function(x) maxY-x)

  #Create a data frame of case attachment frequency and case age
  df <- data.frame(Age = age, Frequency = frequency)

  #Model new case attachment frequency as a function of case age
  fit <- glm(Frequency~Age, data = df, "poisson")
  
  return(fit)
}

#Obtains an estimate of cluster growth based on the ages of cases within
forecast <- function(inG, full=F) {
  #@param inG: A graph cut based on a threshold distance, with the latest casses representing New cases (ie. Upcoming cases)
  #@param full: An option determining whether or not this is the growth estimate for a fully saturated model
  #@return: An attribute for clu representing the estimation of their growth based on case ages
  
  #Obtain the new year, from which we will measure growth
  newY <- max(V(inG)$year)
  
  #In the case of a full model, we need to remove all edges between prevvious cases
  if (full) {
    inG <- inG - E(inG)[(V(inG)[year<newY])%--%(V(inG)[year<newY])]
  }
  
  #Obtain a model of  case growth estimates - predicting growth via the age of member cases of a cluster
  fit <- linkFreq(inG)
  intercept <- unname(fit$coefficients[1])
  slope <- unname(fit$coefficients[2])
  
  #Establish a function which uses the model to obtain an individual cases propensity to grow
  predMod <- function(x) {
    res <- slope*(newY - x) + intercept
    inc <- length(V(inG)[year==x])
    if (res <= 0) {
      return(0)
    }
    else {
      return(res/inc)
    }
  }
  
  #Assign a predicted growth value to each member of the graph
  V(inG)$freq <- sapply(V(inG)$year, function(x) sum(sample(c(1,0), 1, prob = c(predMod(x), 1-predMod(x)))))
  
  #Obtain cluster information
  clu <- components(inG)

  #Obtain a prediction of growth based on the predicted growth of a cluster's members  
  forecast <- sapply(1:clu$no, function(x) {
    members <- names(clu$membership[unname(clu$membership)==x])
    memV <- V(inG)[name%in%members]
    presMV <- memV[year<newY]
    sumFreq <- sum(presMV$freq)
    return(sumFreq)
  })
  
  return(forecast)
}

#Obtains the Growth of clusters based on which clusters hold new cases
growF <- function(inG) {
  #@param inG: A subG gaph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: Cluster information including growth and estimated growth for the Present year (ie. The year before the newest year in inG)
  
  #Obtain the new year
  newY <- max(V(inG)$year)
  
  #obtain a forecast based off of case age for this sub graph
  forecast <- forecast(inG, full=T)
  
  #Define vertices as new cases and present cases
  presV <- V(inG)[year<newY]
  
  #Filter out all edges between present cases (disaggregation)
  es <- E(inG)[presV %--% presV]
  inG <- inG - es
  newV <- V(inG)[year==newY]
  presV <- V(inG)[-newV]
  
  #Obtain cluster information
  clu <- components(inG)
  
  #Assign the previously established forecast to the cluster information
  clu$forecast <- forecast
  
  #Assign the growth of individual clusters (of size 1), based on the newly formatted graph
  clu$growth <- sapply(1:clu$no, function(x){
    members <- names(clu$membership[unname(clu$membership)==x])
    memV <- V(inG)[name%in%members]
    newV <- memV[year==newY]
    return(length(newV))
  })
  
  return(clu)
}

#Obtains the Growth of clusters based on which clusters hold new cases
grow <- function(inG) {
  #@param inG: A subG gaph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: Cluster information including growth and estimated growth for the Present year (ie. The year before the newest year in inG)
  
  #Obtain the newest date
  newY <- max(V(inG)$year)
  
  #Obtain cluster information
  clu <- components(inG)
  
  #Assign the previously established forecast to the cluster information
  clu$forecast <- forecast(inG)
  
  #obtain the number of new cases
  clu$inc <- length(V(inG)[year==newY])

  #Assign cluster growth based on number of new cases embedded in clusters 
  clu$growth <- sapply(1:clu$no, function(x) {
    members <- names(clu$membership[unname(clu$membership)==x])
    memV <- V(inG)[name%in%members]
    newMV <- memV[year==newY]
    return(length(newMV))
  })
  
  return(clu)
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

#____________________________________________________________________________________________________________________________sub#

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

g <- g-V(g)[year==2013]

#Obtain the range of years and the maximum input year
years <- levels(factor(V(g)$year))
y <- max(years)

#Obtain edge id's of all of the shortest edge lengths from new cases (A new case can only be linked to 1 case)
bridgeE <- E(g)[V(g)[year==y]%--%V(g)[year<y]]
closeE <- unname(sapply(V(g)[year==y], function(x) {
  xE <- bridgeE[inc(x)]
  if(length(xE)==0) {
    return(NULL)
  }
  else {    
    closest <- xE[Distance == min(Distance)]
    return (closest[1])
  }
}))
closeE <- unlist(closeE[!vapply(closeE, is.null, logical(1))])
farE <- difference(bridgeE, E(g)[closeE])


#Filter out all edges except for the closest edges
if(!is.null(closeE)){
  g <- g- farE
}

#Initialize the dataframe of results
res <- {}
cutoffs <- seq(0, 0.05, 0.001)

#Generate growth data for each cutoff in a series of cutoffs
for (d in cutoffs) {
  print(noquote(paste0(as.integer(d/max(cutoffs)*100), "%")))  #Progress tracking
  
  #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
  subG <- subGraph(g,y,d)

  #Obtain growth based on a restricted model
  growth <- grow(subG) 
  
  #Obtain growth based on a full model
  growthF <- growF(subG)

  #Group the full and restricted growth models in a list
  l <- list(growth, growthF)
  
  #Add to growing dataframe of results
  res <- cbind(res, l)  
} 

#Label data
rownames(res) <- c("Restricted", "Full")
colnames(res) <- cutoffs

#Save data in accessable file
save(res, file = paste0(gsub("\\..*", "", args), "GD.RData"))