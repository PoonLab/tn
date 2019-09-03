library(dplyr,verbose = FALSE)
require(clmp)

#Import Tree Data and annotate with sequence ID and Time
#Sequences must be dated with the date separated from the id by '_'. 
##TO-DO: Currently only accepts year dates. Work to allow more specific dates. 
impTree <-function(iFile, minNS=63){
  #@param iFile: The name/path of the input file (expecting a newick file)
  #@param minNS: The minimum number of acceptible new sequences By default we keep this high. Passed to sizeCheck.
  #@preturn: An ape tree object with associated lists of sequence ID and Time

  #Creating an ape tree object from the newick file
  t <- read.tree(iFile)
  
  #Obtain lists of sequence ID and Time
  temp <- sapply(t$tip.label, function(x) strsplit(x, '_')[[1]])
  ids <- temp[1,]
  times <- as.numeric(temp[2,])
  
  #Assign those lists to the tree
  t$Time <- times
  t$ID <- ids
  
  #If the newest timepoint contains a small number of sequences, we remove that timepoint from consideration
  while(length(t$Time[t$Time==max(t$Time)]) <= minNS) {
    t <- tFilt(t, max(t$Time)-1)
  }
  
  return(t)
}

#Import TN93 Distance Data for the original sequence data which generated the tree
impTN93Dist <- function(iFile) {
  #@param iFile: The name/path of the input file (expecting tn93 output csv)
  #@return: All the distances between sequences formatted as an edgelist
  
  #From the input file, a tn93 output file. This
  idf <- read.csv(iFile, stringsAsFactors = F)
  temp1 <- sapply(idf$ID1, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp3 <- sapply(idf$ID2, function(x) (strsplit(x,'_')[[1]])[[1]])
  
  #Create data frame of edges (ie. Vertex interactions)
  el <- data.frame(ID1=as.character(temp1), ID2=as.character(temp3),
                   Distance = as.numeric(idf$Distance), stringsAsFactors= F)
  
  return(el)
}

#A simple function, removing tips that sit above a maximum time point (@param: maxT)
tFilt <- function(iT, maxT) {
  maxT <- max(iT$Time)
  maxTi <- which(iT$Time==maxT)
  iT <- drop.tip(iT, maxTi)
  iT$Time <- iT$Time[-maxTi]
  iT$ID <- iT$ID[-maxTi]
  return(iT)
}

#Filter the edges coming from new cases such that the new cases have no edges to eachother and only one edge leading from them to old cases
#This is a simplification to help resolve the merging involved in cluster growth and to prevent growth by whole clusters.
clsFind <- function(nC, oC, Dist) {
  #@param nC: A list of cases from the newest time point with cluster info obtained through clmp() 
  #@param oC: A list of all cases excluding the newest time point with clusters info through clmp() 
  #@param Dist: Genetic Distance information from impTN93Dist()
  #return: A list of the closest non-new neighbour for each new cases
  
  #Remove unclustered cases from consideration
  nC <- subset(nC, Cluster>0)
  oC <- subset(oC, Cluster>0) 
  
  #Find the closest non-new neighbour for each new case
  #Neighbours are defined as cases that share a cluster together
  closeNeighbs <- sapply(1:nrow(nC), function(i){

    #Obtain the case and all of it's "old" neighbours
    x <- nC[i,]
    ioNeighb <- subset(oC, Cluster==x$Cluster & Time<max(nC$Time))
    
    #Obtain the distances of all old members of case x's cluster
    iDist <- subset(Dist, ID1%in%x$ID | ID2%in%x$ID)
    iDist <- subset(iDist, ID1%in%ioNeighb$ID | ID2%in%ioNeighb$ID)
    iDist$ID <- c(iDist$ID1,iDist$ID2)[c(iDist$ID1,iDist$ID2)%in%ioNeighb$ID]
    iDist <- iDist[,c("Distance", "ID")]
    
    #Catching the instance where a given new case has no old neighbours
    if(nrow(iDist)>0) {
      iMin <- subset(iDist, Distance==min(Distance))[1,]$ID
    } else {
      iMin <- ""
    }
    
    return(iMin)
  })
  
  #Remove instances where a given new case has no close neighbours
  closeNeighbs <- closeNeighbs[!closeNeighbs%in%""]
  
  return(closeNeighbs)
}

#Simulate the growth of clusters, showing the difference in cluster size between the newest and the penultimate time point
#The frame of reference for clusters is the penultimate year, simulating one making forcasting decisions based on one time point and validating them with the next
simGrow <- function(nClu, oClu, Dist) {
  #@param nC: The results of a clmp() run on a the inputted tree 
  #@param oC: The results of a clmp() run on a the inputted tree with the tips from the newest timepoint removed 
  #@param Dist: Genetic Distance information from impTN93Dist(). Passed to clsFind()
  #@return: The actual growth and cluster information
  
  #Extract a data frame of individual cases from the clmp() runs, separating new cases and non-new cases
  nC <- data.frame(ID=nClu$ID, Time=nClu$Time, Cluster=as.numeric(head(nClu$clusters, (length(nClu$clusters)+1)/2)), stringsAsFactors=F)
  nC <- subset(nC, Time==max(nC$Time))
  oC <- data.frame(ID=oClu$ID, Time=oClu$Time, Cluster=as.numeric(head(oClu$clusters, (length(oClu$clusters)+1)/2)), stringsAsFactors=F)
  
  #Re-difine the clusters such that cases which would be singletons are considered to have their own cluster
  oC[oC$Cluster==0,]$Cluster <- seq((max(oC$Cluster)+1), nrow(oC[oC$Cluster==0,])+(max(oC$Cluster)))
  
  #A filter to resolve the idea of merging clusters (ie. Cases in two different clusters at one time point, but the same cluster at a more recent one) 
  closeNeighbs <- clsFind(nC, oC, Dist)
  
  #Calculate cluster Growth as a table, by counting the number of close, new neighbours per cluster
  #This is meant to respond to the question "Which old cluster would these new cases be in?"
  growth <- table(oC$Cluster)
  growth[names(growth)] <- rep(0,length(growth))
  posGrowth <- table(sapply(closeNeighbs, function(id){subset(oC, ID%in%id)$Cluster}))
  growth[names(posGrowth)] <- unname(posGrowth)
  
  return(growth)
}

#Obtains some likelihood data in order to weight cases based on their recency  
likData <- function(oClu, Dist) {
  #@param iClu: The results of a clmp() run on a the inputted tree with the tips from the newest timepoint removed 
  #@param Dist: Genetic Distance information from impTN93Dist(). Passed to clsFind()
  #@return: A data frame of "Positives" (related cases) between one time point and another, annotated with the number of possible total positives.
  #         Each positive is also annotated with it's time point difference (between time points)
  
  #Extract a data frame of individual cases from the clmp() run as well as the time points
  c <- data.frame(ID=oClu$ID, Time=oClu$Time, Cluster=as.numeric(head(oClu$clusters, (length(oClu$clusters)+1)/2)), stringsAsFactors=F)
  times <- as.numeric(levels(as.factor(c$Time)))
  
  #Obtain subsets of "Positives" (related cases) between one time point and another, annotated with the number of possible total positives.
  #Positives are only considered if, from the perspective of the newest time point, there are no old cases with a closer TN93 distance measurement
  ageD <- bind_rows(lapply(tail(times,-1), function(x) {
    
    #Separating new cases and non-new cases, paying special attention to the time difference between cases
    nC <- subset(c, Time==x)
    oC <- subset(c, Time<x) 
    tDiffs <- max(nC$Time) - as.numeric(levels(as.factor(oC$Time)))
    
    #The number of total possible positives after clsFilter() is equal to the length of nC, thus we assign that as the total
    nTot <- nrow(nC)
    clsIDs <- clsFind(nC, oC, Dist)
    
    #Sort time difference, total possible positives and total actual positives (at a given time difference) into a data frame
    pos <- table(subset(oC, ID %in% clsIDs)$Time)
    totPos <- rep(0,length(tDiffs))
    totPos[max(nC$Time)-as.numeric(names(pos))] <- unname(pos)
    tdD <- data.frame(Positive=totPos, Total=nTot, tDiff=tDiffs)
    
    return(tdD)
  }))
  
  
  return(ageD)
}

#Analyze a given clustered Tree to establish the difference between the performance of two different models
#Performance is defined as the ability for cluster growth to fit a predictive model.
clmpAnalyze <- function(iT, Dist, nrates=2) {
  #@param iT: An inputted tree representing the total dataset of interest (after some necessary filtering)
  #@param Dist: Genetic Distance information from impTN93Dist() (to pass to likData())
  #@param nrates: The number of rates for clmps MMPP model (to pass to clmp())
  #@return: A summary of cluster information, case information and predictive model performance information for clmp at a given parameter

  #Run clmp() on the tree up to the newest time poiunt and the penultimate time point
  oT <- tFilt(iT, max(iT$Time)-1)
  nClu <- clmp(iT, nrates)
  oClu <- clmp(oT, nrates)
  
  #Obtain likelihood data to inform a weighted model from the penultimate time point
  ageD <- likData(oClu, Dist)
  
  #Extract a data frame of individual cases from the clmp() run as well as the time points
  c <- data.frame(ID=oClu$ID, Time=oClu$Time, Cluster=as.numeric(head(oClu$clusters, (length(oClu$clusters)+1)/2)), stringsAsFactors=F)
  
  #Obtain a model of case clustering frequency as predicted by individual case age
  mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageD, family='binomial')
  
  #Use this to weight cases for a time-informed model
  c$Weight <- predict(mod, data.frame(tDiff=max(iT$Time)-c$Time), type='response')
  
  #Re-difine the clusters such that cases which would be singletons are considered to have
  c[c$Cluster==0,]$Cluster <- seq((max(c$Cluster)+1), nrow(c[c$Cluster==0,])+(max(c$Cluster)))
  
  #Calculate growth from the two clmp runs and cluster size from the penultimate clmp run
  growth <- simGrow(nClu, oClu, Dist)
  csize <- table(c$Cluster)

  
  #Create two data frames from two predictive models, one based on absolute size (NULL) and our time-informed model
  df1 <- data.frame(Growth = as.numeric(growth), Pred = sapply(names(csize), function(x) { sum(subset(c, Cluster==as.numeric(x))$Weight) }))
  df2 <- data.frame(Growth = as.numeric(growth), Pred =  as.numeric(csize) * (sum(growth)/sum(csize)))
  fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
  fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
  
  #Summarize the GAIC, cluster info summary and package it all into a final result
  cSum <- data.frame(csize=as.numeric(csize), NullPred=df2$Pred, PropPred=df1$Pred, Growth=as.numeric(growth))
  gaic <- fit2$aic - fit1$aic 
  res <- list(c=c, cSum=cSum, GAIC=gaic, PropMod=fit1, NullMod=fit2)
  
  return(res)
}

###### Testing ######
#Import Data
TN93File <- "~/Data/Seattle/tn93StsubB.txt" 
treeFile <- "~/Data/Seattle/analysis/FTStsubB.nwk"

#Run Code
Dist <- impTN93Dist(TN93File)
t <- impTree(treeFile)
clmpAnalyze(t, Dist)