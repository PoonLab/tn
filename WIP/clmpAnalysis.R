#A process which interprets clmp cluster data such that it is comparable to tn93 cluster data for scoring the effectiveness of clustering from tn93 Output.
#USAGE: Rscript clmpAnalysis FastTreeOutput.nwk
library(dplyr,verbose = FALSE)
require(clmp)

#Import Distance Data
impTN93Dist <- function(iFile) {
  #From the input file, a tn93 output file. This
  idf <- read.csv(iFile, stringsAsFactors = F)
  temp1 <- sapply(idf$ID1, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp3 <- sapply(idf$ID2, function(x) (strsplit(x,'_')[[1]])[[1]])
  
  #Create data frame of edges (ie. Vertex interactions)
  el <- data.frame(ID1=as.character(temp1), ID2=as.character(temp3),
                   Distance = as.numeric(idf$Distance), stringsAsFactors= F)
  
  return(el)
}

#Check the size to make sure the new year is big enough
sizeCheck <- function(iT) {

  maxT <- max(iT$Time)
  maxTi <- which(iT$Time==maxT)
  
  while(length(maxTi)<63) {
    iT <- drop.tip(iT, maxTi)
    iT$Time <- iT$Time[-maxTi]
    iT$ID <- iT$ID[-maxTi]
    
    maxT <- max(iT$Time)
    maxTi <- which(iT$Time==maxT)
  }
  
  return(iT)
}

#Cut the Newest Year
cutTime <- function(iT) {
  
  iT <- sizeCheck(iT)
  
  maxT <- max(iT$Time)
  maxTi <- which(iT$Time==maxT)
  
  iT <- drop.tip(iT, maxTi)
  iT$Time <- iT$Time[-maxTi]
  iT$ID <- iT$ID[-maxTi]
  
  iT <- sizeCheck(iT)
  
  return(iT)
}

#Import Tree Data and annotate with ID and Time
impTree <-function(iFile){
  #args = commandArgs(trailingOnly = T)
  t <- read.tree(iFile)
  
  #Establish a set of node ids coupled with collection dates
  temp <- sapply(t$tip.label, function(x) strsplit(x, '_')[[1]])
  ids <- temp[1,]
  times <- as.numeric(temp[2,])
  
  t$Time <- times
  t$ID <- ids
  
  t <- sizeCheck(t)
  
  return(t)
}

simGrow <- function(iT, Dist=Dist, nrates=2) {

  nT <- iT
  oT <- cutTime(iT)
  
  #Establish new and old clusters
  nRes <- clmp(nT, nrates)
  oRes <- clmp(oT, nrates)
  
  nC <- data.frame(ID=nRes$ID, Time=nRes$Time, Cluster=as.numeric(head(nRes$clusters, (length(nRes$clusters)+1)/2)), stringsAsFactors=F)
  oC <- data.frame(ID=oRes$ID, Time=oRes$Time, Cluster=as.numeric(head(oRes$clusters, (length(oRes$clusters)+1)/2)), stringsAsFactors=F)
  
  oC[oC$Cluster==0,]$Cluster <- seq((max(oC$Cluster)+1), nrow(oC[oC$Cluster==0,])+(max(oC$Cluster)))
  
  #Calculate growth (through closest membership)
  niC <- subset(nC, Time==max(nC$Time))
  oiC <- subset(nC, Time<max(nC$Time))
  closeNeighbs <- clsFilt(niC, oiC, Dist)
  
  growth <- table(oC$Cluster)
  growth[names(growth)] <- rep(0,length(growth))
  posGrowth <- table(sapply(closeNeighbs, function(id){subset(oC, ID%in%id)$Cluster}))
  growth[names(posGrowth)] <- unname(posGrowth)
  
  return(growth)
}

bpeFreq <- function(iClu, Dist=Dist) {
  
  clu <- data.frame(ID=iClu$ID, Time=iClu$Time, Cluster=as.numeric(head(iClu$clusters, (length(iClu$clusters)+1)/2)), stringsAsFactors=F)
  times <- as.numeric(levels(as.factor(clu$Time)))
  
  eCounts <- bind_rows(lapply(tail(times,-1), function(x) {

    nC <- subset(clu, Time==x)
    oC <- subset(clu, Time<x) 
    nTot <- nrow(nC)
    tDiffs <- max(nC$Time) - as.numeric(levels(as.factor(oC$Time)))
    clsIDs <- clsFilt(nC, oC, Dist)
    pos <- table(subset(oC, ID %in% clsIDs)$Time)
    totPos <- rep(0,length(tDiffs))
    totPos[max(nC$Time)-as.numeric(names(pos))] <- unname(pos)
    data.frame(Positive=totPos, Total=nTot, tDiff=tDiffs)
    
  }))
  
  
  return(eCounts)
}

#TO-DO: This is a Slow Function (thus slowing down eCounts and simGrow).
clsFilt <- function(nC, oC, Dist=Dist) {
  
  nC <- subset(nC, Cluster>0)
  oC <- subset(oC, Cluster>0) 
  
  #For every new case, find it's closest neighbour from the old cases
  closeNeighbs <- sapply(1:nrow(nC), function(i){
    x <- nC[i,]
    
    #The list of old neighbours to new case "x"
    ioNeighb <- subset(oC, Cluster==x$Cluster & Time<max(nC$Time))
    
    #Obtain the distances of all old members of case x's cluster
    iDist <- subset(Dist, ID1%in%x$ID | ID2%in%x$ID)
    iDist <- subset(iDist, ID1%in%ioNeighb$ID | ID2%in%ioNeighb$ID)
    iDist$ID <- c(iDist$ID1,iDist$ID2)[c(iDist$ID1,iDist$ID2)%in%ioNeighb$ID]
    iDist <- iDist[,c("Distance", "ID")]
    
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

clmpAnalysis <- function(iT, Dist=Dist, nrates=2) {

  oT <- cutTime(iT)
  clu <- clmp(oT, nrates)
  ageD <- bpeFreq(clu, Dist)
  
  clu <- data.frame(ID=clu$ID, Time=clu$Time, Cluster=as.numeric(head(clu$clusters, (length(clu$clusters)+1)/2)), stringsAsFactors=F)
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case age
  #Use this to weight cases by age
  mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageD, family='binomial')
  clu$Weight <- predict(mod, data.frame(tDiff=max(iT$Time)-clu$Time), type='response')
  
  growth <- simGrow(iT,Dist, nrates)
  clu[clu$Cluster==0,]$Cluster <- seq((max(clu$Cluster)+1), nrow(clu[clu$Cluster==0,])+(max(clu$Cluster)))
  csize <- table(clu$Cluster)

  #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
  df1 <- data.frame(Growth = as.numeric(growth), Pred = sapply(names(csize), function(x) { sum(subset(clu, Cluster==as.numeric(x))$Weight) }))
  df2 <- data.frame(Growth = as.numeric(growth), Pred =  as.numeric(csize) * (sum(growth)/sum(csize)))
  fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
  fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
  
  gaic <- fit2$aic - fit1$aic 
  
  return(gaic)
}

###### MAIN ######
#Import Data
TN93File <- "~/Data/Seattle/tn93StsubB.txt" 
treeFile <- "~/Data/Seattle/analysis/FTStsubB.nwk"

Dist <- impTN93Dist(TN93File)
t <- impTree(treeFile)

clmpAnalysis(t, Dist, nrates = 4)
