#A process which interprets clmp cluster data such that it is comparable to tn93 cluster data for scoring the effectiveness of clustering from tn93 Output.
#USAGE: Rscript clmpAnalysis FastTreeOutput.nwk

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
  oT <- cutTime(oT)
  
  #Establish new and old clusters
  nRes <- clmp(nT, nrates = 2)
  oRes <- clmp(oT, nrates = 2)
  
  nC <- data.frame(ID=nT$ID, Time=nT$Time, Cluster=head(nRes$clusters, (length(nRes$clusters)+1)/2))
  oC <- data.frame(ID=oT$ID, Time=oT$Time, Cluster=head(oRes$clusters, (length(oRes$clusters)+1)/2))
  
  oC[oC$Cluster==0,]$Cluster <- seq((max(oC$Cluster)+1), nrow(oC[oC$Cluster==0,])+(max(oC$Cluster)))
  
  #Calculate growth (through closest membership)
  iC <- nC
  niC <- subset(iC, Time==max(iC$Time) & Cluster>0)
  closeNeighbs <- sapply(1:nrow(niC), function(i){
    x <- niC[i,]
    
    
    ioNeighb <- subset(iC, Cluster==x$Cluster & Time<max(iC$Time))
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
  
  growth <- table(oC$Cluster)
  growth[names(growth)] <- rep(0,length(noGrowth))
  posGrowth <- table(subset(oC, ID%in%closeNeighbs)$Cluster)
  growth[names(posGrowth)] <- unname(posGrowth)
  
  return(growth)
}

bpeFreq <- function(iT, Dist=Dist) {
  iT <- cutTime(iT)
  res <- clmp(iT, nrates = 2)
  
  c <- data.frame(ID=iT$ID, Time=iT$Time, Cluster=head(res$clusters, (length(res$clusters)+1)/2))
  
  for (i in )
}

#Import Data
TN93File <- "~/Data/Seattle/tn93StsubB.txt" 
treeFile <- "~/Data/Seattle/analysis/FTStsubB.nwk"
Dist <- impTN93Dist(TN93File)
t <- impTree(treeFile)






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