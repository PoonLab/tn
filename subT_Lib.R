require("ape") # Possibly Unneeded
require("phangorn")
library("dplyr",verbose = FALSE)

#Import Tree Data and annotate with sequence ID and Time
#Sequences must be dated with the date separated from the id by '_'. 
##TO-DO: Currently only accepts year dates. Work to allow more specific dates. 
impTree <-function(tFile){
  #@param iFile: The name/path of the input file (expecting a newick file)
  #@preturn: An ape tree object with associated lists of sequence ID and Time
  
  #Creating an ape phylogeny object from the newick file, store in a greater list "t" as "p" for phylogeny
  t <- list()
  t$p <- read.tree(tFile)
  
  #Obtain lists of sequence ID and Time
  temp <- sapply(t$p$tip.label, function(x) strsplit(x, '_')[[1]])
  ids <- unname(temp[1,])
  times <- as.numeric(temp[2,])
  t$v <- data.frame(ID=ids, Time=times, stringsAsFactors = F)
  tips <- 1:nrow(t$v)
  
  #Summarize internal branch length information
  t$e <- list()
  t$e$dist <- dist.nodes(t$p)
  t$e$tDiff <- sapply(tips, function(i){
    t1 <- t$v$Time[i]
    sapply(tips, function(j){ t1 - t$v$Time[j] })
  })
  
  #Summarize information by node, representing the potential to cluster by subtree
  #Obtain the mean branch length under each node
  ##TO-DO: Add bootstrap here. Requires bootstrap to be added in tree creation.
  nodes <- (max(tips)+2):(max(tips)*2-2) 
  des <- lapply(nodes, function(x){Descendants(t$p,x,"all")})
  meanDist <- sapply(des, function(x) {
    x <- sort(x[which(x%in%tips)]) 
    m <- mean(unlist(sapply(head(x,-1), function(tip){t$e$dist[tip,x[which(x>tip)]]})))
    return(m)
  })
  
  #The minimum retrospective edge to each tip (except those at the earliest time point)
  t$v <- cbind(t$v, bind_rows(lapply(tips, function(i){
    
    if(t$v$Time[i]==min(t$v$Time)){
      return(data.frame(ClosestV=NA, Distance=NA, tDiff=NA))
    }
    else {
      inc <- t$e$dist[i,-i]
      iTD <- t$e$tDiff[i,-i]
      ret <- inc[which(iTD<0)]
      
      minE <- which(ret==min(ret))[[1]] 
      minD <- ret[minE]
      minTD <- (iTD[which(iTD<0)])[minE]
      
      df <- data.frame(ClosestV=minE, Distance=minD, tDiff=minTD)
      return(df)
    }

  })))
  
  #Set up a node-summary data frame
  t$n <- list()
  t$n$des <- des
  t$n$agg <- data.frame(ID=nodes, mDist=meanDist)
  names(t$n$des) <- t$n$agg$ID
  
  return(t)
}

#After simulating the growth of trees by placing recent sequences as tips on a fixed ML tree
growthSim <- function(iT, gFile) {
  
  #Obtain a set of trees with tips added
  ts <- read.tree(gFile)
  
  #The nodes in iT associated with new cases in t
  iT$g<- bind_rows(lapply(ts, function(t) {
    
    #Find the terminal branch length of the new tip and the node it is adjacent to
    temp <- sapply(t$tip.label, function(x) strsplit(x, '_')[[1]])
    times <- as.numeric(temp[2,])
    ids <- temp[1,]
    nTip <- which(times>max(iT$v$Time))
    nAdj <- Siblings(t, nTip)
    dist <- t$edge.length[which(t$edge[,2]==nTip)]
    
    #If the sibling is a tip, returns the sibling's parent in iT
    if(nAdj<=length(t$tip.label)) {
      adj <- which(iT$v$ID%in%ids[nAdj])
      parent <- iT$p$edge[which(iT$p$edge[,2]==adj),1]

      return(data.frame(ID=ids[nTip], node=parent, dist=dist)) 
    }
    
    else{
      return(data.frame(ID=ids[nTip], node=nAdj, dist=dist))
    }
  }))
  
  return(iT)
}

#Cluster using a modified subtree-based method
STClu <- function(iT, maxD) {

  #Identify potential clusters by some criterion as well as tips within those potential clusters
  cluNames <- as.character(subset(iT$n$agg, mDist<maxD)$ID)
  cluTips <- unlist(iT$n$des[cluNames]) 
  cluTips <- unique(as.numeric(cluTips[cluTips<length(iT$v$ID)]))
  
  #A loop to remove any clusters which are subclusters of a larger cluster
  temp <- vector()
  
  while(length(cluTips)>0) {
    
    tip <- cluTips[1]
    tClu <- cluNames[sapply(iT$n$des[cluNames], function(iClu) {tip%in%iClu})]
    tCluS <- sapply(tClu, function(cluName) {length(iT$n$des[[cluName]])})
    tClu <- tClu[which(tCluS==max(tCluS))[[1]]]
    temp <- c(temp,tClu)
    cluTips <- cluTips[-which(cluTips%in%iT$n$des[[tClu]])]
    
  }
  
  cluNames <- temp

  #Remove internal nodes for summary. Only cases are of interest
  clu <- lapply(cluNames, function(i){

    #Obtain Tips already in Cluster and the mean distance of this cluster
    iDes <- iT$n$des[[i]]
    iTipDes <- iT$v$ID[iDes[which(iDes<length(iT$v$ID))]]
    mDist <- subset(iT$n$agg, ID==as.numeric(i))$mDist
    n <- length(iTipDes)
    
    #New cases associated with an internal node within the potential cluster
    gTips <- subset(iT$g, iT$g$node%in%c(iDes,as.numeric(i)))
    
    if(nrow(gTips)>0) {
      #Obtaining only new cases satisfying distance requirement
      breakCon <- sapply(1:nrow(gTips), function(j) {
        gTip <- gTips[j,]
        gDist <- iT$e$dist[gTip$node, iDes[which(iDes<length(iT$v$ID))]] + gTip$dist
        bDist <- (sum(gDist) + mDist*factorial(n)) / factorial(n+1)
        
        return(bDist>maxD)
      })
      
      gTips <- gTips[!breakCon,]
    }
    
    if(nrow(gTips)==0) {
      gTips <- NULL
    }
  
    return(c(iTipDes, gTips$ID))
  })
  
  iT$c <- list()
  
  #Add singletons and summerize clusters
  sing <- as.list(setdiff(c(iT$v$ID, iT$g$ID), unlist(clu)))
  iT$c$membership <- c(clu,sing)
  temp <- sapply(iT$c$membership, function(x){
    newCases <- length(which(x%in%iT$g$ID))
    oldCases <- length(x) - newCases
    
    return(c(oldCases,newCases))
  }) 
  
  #For growth, record the number of old and new cases in the clusters which contain old cases
  colnames(temp) <- 1:length(temp[1,])
  iT$c$growth <- temp[,which(temp[1,]>0)]
    
  return(iT)
}

#Obtain GAIC at several different cutoffs
GAICRun <- function(iT, cutoffs) {
  
  #Clustering Test
  res <- lapply(cutoffs, function(x){
    print(x)
    subT <- STClu(iT, x)
    
    times <- table(subT$v$Time)
    
    #Get a rating, accounting for the potential biases created by old outbreak
    bias <- table(iT$v$Time+iT$v$tDiff) / 
      table(iT$v$Time)[which(names(table(iT$v$Time))%in%names(table(iT$v$Time+iT$v$tDiff)))]
    bias <- sapply(1:length(bias), function(i){bias[i]/ mean(bias[-i])})
 
    
    #Obtain the age-based data in order to weight vertices
    ageD <- bind_rows(lapply(2:length(times), function(i){
      
      t <- times[i]
      tV <- subset(subT$v, Time==as.numeric(names(t)))
      
      #Obtain a set of tDiffs (retrospective time lag) and positives at those time lags
      tDiffs <- as.numeric(names(times))-as.numeric(names(t))
      otBias <- bias[names(times[which(tDiffs<0)])]
      pos <-sapply(tDiffs[tDiffs<0], function(td){nrow((subset(tV, (tDiff==td)&(Distance<x))))})
      
      #Total Possible positives (attempts)
      tot <- rep(unname(t), length(pos))
      df <- data.frame(Total=tot, Positives=pos, tDiff=abs(tDiffs[tDiffs<0]), otBias=as.numeric(otBias))

      return(df)
    }))
    
    #Weighting cases using a model based on case recency
    mod <- glm(cbind(Positives, Total) ~ tDiff + otBias, data=ageD, family='binomial')
    subT$v$Weight <- predict(mod, type='response', 
                             data.frame(tDiff=max(subT$v$Time)+1-subT$v$Time, 
                                        otBias=sapply(subT$v$Time, function(x){bias[as.character(x)]})))
    
    #Weight clusters using their member case weights (0 represents a cluster of all new cases, and is ignored)
    cluWeights <- sapply(subT$c$membership, function(c){sum(subset(subT$v, ID%in%c)$Weight)})
    cluWeights <- cluWeights[cluWeights>0]
    
    #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
    df1 <- data.frame(Growth = subT$c$growth[2,], Pred = cluWeights) 
    df2 <- data.frame(Growth = subT$c$growth[2,], Pred = subT$c$growth[2,]*(subT$c$growth[1,]/sum(subT$c$growth[1,])))
    fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
    fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
    
    subT$a <- list()
    
    #Save, gaic, model and age data as part of the output
    subT$a$gaic <- fit1$aic-fit2$aic
    subT$a$propFit <- fit1
    subT$a$nullFit <- fit2
    subT$a$mod <- mod
    subT$a$ageD <- ageD

    return(subT)
  })
  
  return(res)
  
}

###############Testing

tFile <- "~/Data/Seattle/SeattleB_RAxML/RAxML_result.Tree" 
gFile <- "~/Data/Seattle/SeattleB_pplacer/st.tre"

#tFile <- "~/Data/Tennessee/TennesseeB_Diag_RAxML/RAxML_result.Tree" 
#gFile <- "~/Data/Tennessee/TennesseeB_Diag_pplacer/tn.tre"

oT <- impTree(tFile) 
oT <- growthSim(oT, gFile)

cutoffs <- seq(0,0.25,0.002)
res <- GAICRun(oT, cutoffs)

gaics <- sapply(res, function(x){x$a$gaic})
aic <- sapply(res, function(x){x$a$propFit$aic})
ageDs <- sapply(res, function(x){x$ageD})


#saveRDS(res, "~/Data/Tennessee/STtnsubB_GD.rds")
