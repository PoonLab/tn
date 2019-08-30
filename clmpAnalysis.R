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

simGrow <- function(nClu, oClu, Dist=Dist, nrates=2) {
  
  nV <- data.frame(ID=nClu$ID, Time=nClu$Time, Cluster=as.numeric(head(nClu$clusters, (length(nClu$clusters)+1)/2)), stringsAsFactors=F)
  oV <- data.frame(ID=oClu$ID, Time=oClu$Time, Cluster=as.numeric(head(oClu$clusters, (length(oClu$clusters)+1)/2)), stringsAsFactors=F)
  
  oV[oV$Cluster==0,]$Cluster <- seq((max(oV$Cluster)+1), nrow(oV[oV$Cluster==0,])+(max(oV$Cluster)))
  
  #Calculate growth (through closest membership)
  niV <- subset(nV, Time==max(nV$Time))
  oiV <- subset(nV, Time<max(nV$Time))
  closeNeighbs <- clsFilt(niV, oiV, Dist)
  
  growth <- table(oV$Cluster)
  growth[names(growth)] <- rep(0,length(growth))
  posGrowth <- table(sapply(closeNeighbs, function(id){subset(oV, ID%in%id)$Cluster}))
  growth[names(posGrowth)] <- unname(posGrowth)
  
  return(growth)
}

bpeFreq <- function(iClu, Dist=Dist) {
  
  v <- data.frame(ID=iClu$ID, Time=iClu$Time, Cluster=as.numeric(head(iClu$clusters, (length(iClu$clusters)+1)/2)), stringsAsFactors=F)
  times <- as.numeric(levels(as.factor(v$Time)))
  
  eCounts <- bind_rows(lapply(tail(times,-1), function(x) {

    nV <- subset(v, Time==x)
    oV <- subset(v, Time<x) 
    nTot <- nrow(nV)
    tDiffs <- max(nV$Time) - as.numeric(levels(as.factor(oV$Time)))
    clsIDs <- clsFilt(nV, oV, Dist)
    pos <- table(subset(oV, ID %in% clsIDs)$Time)
    totPos <- rep(0,length(tDiffs))
    totPos[max(nV$Time)-as.numeric(names(pos))] <- unname(pos)
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
  nClu <- clmp(iT, nrates)
  oClu <- clmp(oT, nrates)
  ageD <- bpeFreq(oClu, Dist)
  
  v <- data.frame(ID=oClu$ID, Time=oClu$Time, Cluster=as.numeric(head(oClu$clusters, (length(oClu$clusters)+1)/2)), stringsAsFactors=F)
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case age
  #Use this to weight cases by age
  mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageD, family='binomial')
  v$Weight <- predict(mod, data.frame(tDiff=max(iT$Time)-v$Time), type='response')
  growth <- simGrow(nClu, oClu, Dist, nrates)
  v[v$Cluster==0,]$Cluster <- seq((max(v$Cluster)+1), nrow(v[v$Cluster==0,])+(max(v$Cluster)))
  csize <- table(v$Cluster)

  #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
  df1 <- data.frame(Growth = as.numeric(growth), Pred = sapply(names(csize), function(x) { sum(subset(v, Cluster==as.numeric(x))$Weight) }))
  df2 <- data.frame(Growth = as.numeric(growth), Pred =  as.numeric(csize) * (sum(growth)/sum(csize)))
  fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
  fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
  
  cSum <- data.frame(csize=as.numeric(csize), NullPred=df2$Pred, PropPred=df1$Pred, Growth=as.numeric(growth))
  
  gaic <- fit2$aic - fit1$aic 
  
  g <- list(v=v, cSum=cSum, GAIC=gaic, PropMod=fit1, NullMod=fit2)
  
  return(g)
}

###### Testing ######
#Import Data
TN93File <- "~/Data/Seattle/tn93StsubB.txt" 
treeFile <- "~/Data/Seattle/analysis/FTStsubB.nwk"

Dist <- impTN93Dist(TN93File)
t <- impTree(treeFile)

#Workaround for a bug ensuring that R nrates can not be set through a variable in a loop
res2 <- clmpAnalysis(t, Dist, nrates = 2)
res3 <- clmpAnalysis(t, Dist, nrates = 3)
res4 <- clmpAnalysis(t, Dist, nrates = 4)
res5 <- clmpAnalysis(t, Dist, nrates = 5)
res6 <- clmpAnalysis(t, Dist, nrates = 6)
res7 <- clmpAnalysis(t, Dist, nrates = 7)
res8 <- clmpAnalysis(t, Dist, nrates = 8)
res9 <- clmpAnalysis(t, Dist, nrates = 9)
res10 <- clmpAnalysis(t, Dist, nrates = 10)


