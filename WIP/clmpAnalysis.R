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

simGrow <- function(iT, Dist=Dist) {

  nT <- iT
  oT <- cutTime(iT)
  
  #Establish new and old clusters
  nRes <- clmp(nT, nrates = 2)
  oRes <- clmp(oT, nrates = 2)
  
  nC <- data.frame(ID=nT$ID, Time=nT$Time, Cluster=head(nRes$clusters, (length(nRes$clusters)+1)/2))
  oC <- data.frame(ID=oT$ID, Time=oT$Time, Cluster=head(oRes$clusters, (length(oRes$clusters)+1)/2))
  
  oC[oC$Cluster==0,]$Cluster <- seq((max(oC$Cluster)+1), nrow(oC[oC$Cluster==0,])+(max(oC$Cluster)))
  
  #Calculate growth (through closest membership)
  niC <- subset(nC, Time==max(nC$Time) & Cluster>0)
  oiC <- subset(nC, Time<max(nC$Time) & Cluster>0)
  closeNeighbs <- clsFilt(niC, oiC, Dist)
  
  growth <- table(oC$Cluster)
  growth[names(growth)] <- rep(0,length(growth))
  posGrowth <- table(sapply(closeNeighbs, function(id){subset(oC, ID%in%id)$Cluster}))
  growth[names(posGrowth)] <- unname(posGrowth)
  
  return(growth)
}

bpeFreq <- function(iClu, Dist=Dist) {
  
  clu <- data.frame(ID=iClu$ID, Time=iClu$Time, Cluster=head(iClu$clusters, (length(iClu$clusters)+1)/2))
  times <- as.numeric(levels(as.factor(clu$Time)))
  
  eCounts <- bind_rows(lapply(tail(times,-1), function(x) {
    nC <- subset(clu, Time==x & Cluster>0)
    oC <- subset(clu, Time<x & Cluster>0)
    nTot <- nrow(nC) 
    clsIDs <- clsFilt(nC, oC, Dist)
    
    tDiffs <- max(nC$Time) - as.numeric(levels(as.factor(oC$Time)))
    pos <- table(subset(oC, ID %in% clsIDs)$Time)
    totPos <- rep(0,length(tDiffs))
    totPos[max(nC$Time)-as.numeric(names(pos))] <- unname(pos)
    data.frame(Positive=totPos, Total=nTot, tDiff=tDiffs)
  }))
  
  return(eCounts)
}

#TO-DO: This is a Slow Function (thus slowing down eCounts and simGrow).
clsFilt <- function(nC, oC, Dist=Dist) {
  
  #For every new case, find it's closest neighbour from the old cases
  closeNeighbs <- sapply(1:nrow(nC), function(i){
    x <- nC[i,]
    
    #The list of old neighbours to new case "x"
    ioNeighb <- subset(oC, Cluster==x$Cluster & Time<max(iC$Time))
    
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
  closeNeighbs <-  
  
  return(closeNeighbs)
}

clmpAnalysis <- function(iT, Dist=Dist) {

  oT <- cutTime(iT)
  clu <- clmp(oT, nrates = 2)
  ageD <- bpeFreq(clu, Dist)
  
  clu <- data.frame(ID=clu$ID, Time=clu$Time, Cluster=head(clu$clusters, (length(clu$clusters)+1)/2))
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case age
  #Use this to weight cases by age
  mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageD, family='binomial')
  clu$Weight <- predict(mod, data.frame(tDiff=max(iT$Time)-clu$Time), type='response')
  
  growth <- simGrow(iT,Dist)
  clu[clu$Cluster==0,]$Cluster <- seq((max(clu$Cluster)+1), nrow(clu[clu$Cluster==0,])+(max(clu$Cluster)))
  csize <- table(clu$Cluster)
  
  #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
  df1 <- data.frame(Growth = growth, Pred = sapply(names(csize), function(x) { sum(subset(clu, Cluster==as.numeric(x))$Weight) }))
  df2 <- data.frame(Growth = growth, Pred =  csize * (sum(growth)/sum(csize)))
  fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
  fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
  
}

#Import Data
TN93File <- "~/Data/Seattle/tn93StsubB.txt" 
treeFile <- "~/Data/Seattle/analysis/FTStsubB.nwk"
Dist <- impTN93Dist(TN93File)
t <- impTree(treeFile)






cc <- clmp(t,nrates = 2)



saveRDS(res, file="BPE.rds")



######################Playing with nRates

l <- lapply(1:10, function(nr){
  ####- TO-DO: Modulate parameters of clump to produce different sets of cluster data
  table(clusters$Cluster)
  return(table(clusters$Cluster))
})

sings <- sapply(l, function(lsub) {unname(lsub["0"])})
maxs <- sapply(l, function(lsub) {unname(max(tail(lsub,-1)))})
means <-  sapply(l, function(lsub) {
  cluA <- tail(lsub,-1)
  cluB <- sum(cluA)-max(cluA)
  cluB/(length(cluA)-1)
})

dfOld <- data.frame(sings,maxs,means)

####- TO-DO: Modulate parameters of clump to produce different sets of cluster data
res <- clmp(t)

#Establish a set of node ids coupled with collection dates
temp <- sapply(res$tip.label, function(x) strsplit(x, '_')[[1]])
ids <- temp[1,]
times <- as.numeric(temp[2,])

clusters <- data.frame(ID=ids, Time=times, Cluster=head(res$clusters, (length(res$clusters)+1)/2))
return(table(clusters$Cluster))

#To make clmp cluster data comparable to tn93 cluster data
####- TO-DO: Pass to tn93 analysis software
clu <- list()
clu$membership <- res$clusters
clu$csize <- unname(table(res$clusters)) #Cluster number 0 is reserved for all singletons 
clu$no <- clu1$csize[[1]] + (length(clu1$csize) - 1)