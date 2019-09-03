##TO-DO: Specify source-file Location
source("~/git/tn/newLib.R")

g <- readRDS("~/Data/Tennessee/analysis/tn93TnsubB_nomet_G.rds")

#remV <- rev(g$v$ID)[1:98]
#g$v <- subset(g$v, !ID%in%remV)
#g$e <- subset(g$e, !ID1%in%remV & !ID2%in%remV)
#g$f <- bpeFreq(g)


#Create a set of subgraphs with sequentially smaller new-years filtered subgraphs
gs  <- lapply(1:(unname(tail(table(g$v$Time),1))-1), function(x){
  iG <- g
  remV <- rev(iG$v$ID)[1:x]
  iG$v <- subset(iG$v, !ID%in%remV)
  iG$e <- subset(iG$e, !ID1%in%remV & !ID2%in%remV)
  iG$f <- bpeFreq(iG)
  
  print(nrow(subset(iG$v, Time==max(iG$v$Time))))
  
  return(iG)
}) 

#Run the GAIC algorithm for all subgraphs
runs <- lapply(1:length(gs) , function(i) {
  
  print(paste0(round(i/length(gs)*100),"%"))
  iG<-gs[[i]]
  
  gaicRun(iG)
})

#Save all growth data in accessable files
saveRDS(runs, file = paste0(oFile, "_SD.rds"))

#Obtain a list of vectors of GAICs for each filtered run
cutoffs <- as.numeric(names(runs[[1]])) 
gaics <- lapply(runs, function(run){sapply(run, function(x) {x$gaic})})

minsLoc <- sapply(gaics, function(x){step*(which(x==min(x))[[1]]-1)}) 
mins <- sapply(gaics, function(x){min(x)}) 

plot(79:1, mins, xlab = "Number of New Cases", ylab="Absolute Minimum GAIC Obtained")

rMet <- readRDS("~/Data/Tennessee/analysis/tn93TnsubB_met_SD.rds")
rNM <- readRDS("~/Data/Tennessee/analysis/tn93TnsubB_nomet_SD.rds")

mins <- sapply(gaics, function(x){min(x)})
