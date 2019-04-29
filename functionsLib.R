##Functions Library 

#Obtains the Growth of clusters based on which clusters hold new cases
embGrow <- function(inG) {
  #@param inG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: Cluster information including growth and estimated growth for the Present year (ie. The year before the newest year in inG)
  
  #Obtain the new year
  newY <- max(V(inG)$year)
  
  #Filter out all edges between present cases (disaggregation)
  presV <- V(inG)[year<newY]
  es <- E(inG)[presV %--% presV]
  inG <- inG - es
  
  #Redifine new cases and present cases
  newV <- V(inG)[year==newY]
  presV <- V(inG)[-newV]
  
  #Obtain cluster information
  clu <- components(inG)
  
  #Assign the growth of individual clusters (of size 1), based on the newly formatted graph
  clu$growth <- sapply(1:clu$no, function(x){
    members <- names(clu$membership[unname(clu$membership)==x])
    memV <- V(inG)[name%in%members]
    newV <- memV[year==newY]
    return(length(newV))
  })
  
  return(clu)
}

#Estimates the growth of clusters based on information from recent years
####- TO-DO: Include include meta-data factors -####
forecastRG <- function(inG, clu, full=F) {
  #@param inG: A subG gaph cut based on a threshold distance, with the latest casses representing New cases (ie. Upcoming cases)
  #@param clu: A set of clusters based on the present year. Obtained from growG()
  #@param full: An option determining whether or not this is the growth estimate for a fully saturated model
  #@return: An attribute for clu representing the past growth of clusters relative to their size. (ie. The predicted absolute growth)
  
  #In the case that this is the fully saturated model
  if (full) {
    inG <- inG - E(inG)
  }
  
  #Obtain a past year to compare too (ideally 5 years before the present year), to establish the recent growth of present clusters
  presY <- max(V(inG)$year)
  minY <- min(V(inG)$year)
  
  #To ensure we don't exceed the bottom limit of the years in data.
  if (presY>(minY+5)) {oldY <- (presY-5)} 
  else {oldY <- minY}
  
  #Difference in past year and present year
  diff <- presY-oldY
  
  #Obtain the cluster sizes of present clusters based only on membership from the old year
  oldMem <- clu$membership[names(clu$membership) %in% V(inG)[year<=oldY]$name]
  oldCsize <- sapply(1:clu$no, function(x) length(oldMem[unname(oldMem)==x]))
  
  #Obtain the Relative, Recent Growth of clusters
  rrG <- (clu$csize-oldCsize) / (diff*sqrt(clu$csize))  
  
  return(rrG)
}

#To handle input as dates instead of years (Pre-Processing for the NA dataset)
if (dates == T) {
  y <- as.Date(temp[2,])
  #Handling day-month-year common format
  yDMY <- as.Date(temp[2,], format="%d-%b-%y")
  y[is.na(y)] <- yDMY[!is.na(yDMY)]
  V(g)$year <- as.integer(as.integer(y) / 120) #Binned into 120 day blocks
}

#For the case of handling missing data
temp <- {}
temp$Frequency <- ageDi$Frequency[!is.nan(ageDi$Frequency)]
temp$Age <- ageDi$Age[!is.nan(ageDi$Frequency)]
ageDi <- temp

## Skeleton Code For plotting From GD results
stat <- sapply(1:(ncol(res)-1), function(x) {
  fit <- res[[1,x]]
  mean(fit$growth)
})
plot(head(colnames(res), -1), stat, ylab = "Mean Cluster Growth" , xlab= "Cutoff Threshold", cex.lab = 1.65, cex.axis = 1.3, cex = 1.5, pch = 20)

#For GAIC edits on already processed Data
LastMin <- function(res,args) {
  UnW <- gaic(res)
  saveRDS(UnW, file = paste0(gsub("\\..*", "", args), "UnW.rds"))
  DisA <- gaic(res, agg=T)
  saveRDS(UnW, file = paste0(gsub("\\..*", "", args), "DisA.rds"))
}

denZero <- function(inG){
  stat <- sapply(as.integer(levels(as.factor(V(in$year))), function(i) {
    gsub <- induced_subgraph(g, V(g)[year==i])
    edge_density(gsub)
  })
  return(stat)
}

LastMin(GDtn, "ColDateData/tnDGD.rds")
LastMin(GDst, "ColDateData/stDGD.rds")
LastMin(GDna, "ColDateData/naDGD.rds")