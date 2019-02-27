#A process which generates cluster growth data as a function of tn93 cutoff threshold
#Creates an external .RData file of paired cluster info sets

### USAGE: Rscript tn93GD.R tn93output.csv ###

####- TO-DO: Cross-Validate with train / test -####

#Import Libraries
library(igraph)
library(dplyr)

####- TO-DO: Use coefficients from a non-linear model to establish age instead of means-####
#library(easynls)

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
    freq <- length(bridgeE) / (length(presV)*length(newV))
    return(freq)
  })
  
  #Assign age to every case
  age <- sapply(years, function(x) maxY-x)

  #Create a data frame of case attachment frequency and case age
  df <- data.frame(Age = age, Frequency = frequency)

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
simGrow <- function(inG, fit, full=F) { 
  
  #Obtain the newest date
  newV <- V(inG)[year==newY]

  #Obtain a forward-Censored Graph
  oldG <- inG - newV
  
  if (full) {
    oldG <- oldG-E(oldG)
  }
  
  #Obtain cluster information
  clu <- components(oldG)

  #Assign cluster growth based on number of new cases embedded in clusters 
  temp <- sapply(1:clu$no, function(x) {
    members <- names(clu$membership[unname(clu$membership)==x])
    memV <- V(inG)[name%in%members]
    bridgeE <- E(inG)[memV%--%newV]
    forecast <- sum(memV$freq) 
    growth <- length(bridgeE)
    return(c(growth,forecast))
  })
  
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

#Creates a graph based on the inputted data frame. The tn93 Distances become edge4 attributes
g <- graph_from_data_frame(input, directed=FALSE, vertices=NULL)

#Adds the ID's and Sample collection years as different vertex attributes for each vertex
####- TO-DO: Assign Dates by Month -####
temp <- sapply(V(g)$name, function(x) strsplit(x, '_')[[1]])
V(g)$name <- temp[1,]
V(g)$year <- as.numeric(temp[2,])

#Obtain the range of years and the maximum input year
years <- as.integer(levels(factor(V(g)$year)))
newY <- max(years)

#Obtain the range of years and the maximum input year
years <- as.integer(levels(factor(V(g)$year)))
newY <- max(years)

## Obtain a set models of case linkage frequency based on age
#__________________________________________________________________________________________________________________________#

#Initialize a list of cases
ldf <- {}

#Initialize a set of cutoffs to observe
cutoffs <- seq(0, 0.05, 0.001)

#Progress tracking
print("Modelling age and cutoff effects on node linkage to new cases...")

#Create a set of frequency data, representing the frequency of new case additions as a function of old case age
for (y in years) {
  #Progress tracking
  print(noquote(paste0(as.integer((1-(newY-y)/length(years))*100) , "%")))
  
  #Obtain a subgraph of cases below a given year and close-filtered edges to the new year
  cfG <- g - V(g)[year>y]
  cfG <- closeFilter(cfG)
  
  #To catch the first year of a set of years (no retrospective years to count to)
  if (y==min(years)) {
    next()
  }
  
  #Obtain a set of link frequencies for this year at various cutoffs
  ldfy <- lapply(cutoffs, function(d) {
    subG <- subGraph(cfG, y, d)
    linkFreq(subG)
  })

  #Add the set of link frequencies to our growing data set
  ldf <- cbind(ldf, ldfy)
}

#Collapse data to a dataframe per cutoff, with case Age and Frequency data for each data frame
ageD <- lapply(1:nrow(ldf), function(x) bind_rows(ldf[x,]))
names(ageD) <- cutoffs

#Save data in accessable file
save(ageD, file = paste0(gsub("\\..*", "", args), "AD.RData"))

## Generate Growth data
#__________________________________________________________________________________________________________________________#

#Obtain a filtered graph for the purposes of measuring growth
g <- closeFilter(g)

#Initialize the dataframe of results
res <- {}

#Progress tracking
print("Modelling cutoff effects on case growth...")

#Generate growth data for each cutoff in a series of cutoffs
for (d in cutoffs) {
  #Progress tracking
  print(noquote(paste0(as.integer(d/max(cutoffs)*100), "%")))
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case age 
  df <- ageD[[as.character(d)]]
  fit <- sapply(levels(factor(df$Age)), function(x) mean(df$Frequency[which(df$Age == x)]))
  
  ####- TO-DO: Actually make fit an exponential model instead of a set of means -####
  #fit <- nlsfit(data=df,model=6, start=c(y0,1))
  
  #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
  subG <- subGraph(g,newY,d)
  
  #Assign a predicted growth value to each member of the graph
  V(subG)$freq <- sapply(V(subG)$year, function(x) {
    
    ####- TO-DO: Use coefficients to establish age instead of means-####
    age <- newY - x
    if (age==0) { return(NaN) }
    else {
      f <- unname(fit[age])
      #weight <- sum(sample(c(1,0), length(V(inG)[year==newY]), replace = T, prob = c(f, 1-f)))
      return(f)
    }
  })

  #Obtain growth based on two models restricted model
  growth <- simGrow(subG, fit) 
  growthF <- simGrow(subG, fit, full=T)
  l <- list(growth, growthF)
  
  #Add to growing dataframe of results
  res <- cbind(res, l)  
} 

#Label data
rownames(res) <- c("Restricted", "Full")
colnames(res) <- cutoffs

#Save data in accessable file
save(res, file = paste0(gsub("\\..*", "", args), "GD.RData"))