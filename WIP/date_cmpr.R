# Originally used in paper to make comparable sets of data (with diagnostic vs. collection dates). Currently in the process of generalizability.
## TO-DO: Unstable, un-tested code. ##
source("comp_Lib.R")
require(R.utils)

#Take in tn93 Distances
iFile <- ""

idf <- read.csv(iFile, stringsAsFactors = F)
temp1 <- sapply(idf$ID1, function(x) (strsplit(x,'_')[[1]])[[1]])
temp2 <- sapply(idf$ID1, function(x) (strsplit(x,'_')[[1]])[[2]])
temp3 <- sapply(idf$ID2, function(x) (strsplit(x,'_')[[1]])[[1]])
temp4 <- sapply(idf$ID2, function(x) (strsplit(x,'_')[[1]])[[2]])

el <- data.frame(ID1=as.character(temp1), t1=as.numeric(temp2), ID2=as.character(temp3), t2=as.numeric(temp4), 
                 Distance = as.numeric(idf$Distance), stringsAsFactors=F)
el$tMax <- pmax(el$t1,el$t2)
el$tDiff <- abs(el$t1-el$t2)
vl <- unique(data.frame(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), stringsAsFactors=F))

#Attatch meta data file
mFile <- "MD.csv"

df <- read.table(mFile,sep='\t', header=T)
vm <- data.frame(ID=as.character(df$CFAR_PID),Time=df$YEAR_OF_HIV_DX, stringsAsFactors=F) #Likely more members than vl
vm <- subset(vm, (ID%in%vl$ID)) #Likely just as many members as vl, but contains missing data
vm <- subset(vm, !is.na(Time)) #Likely just as many members as vl
elm <- subset(el, ID1%in%vm$ID & ID2%in%vm$ID ) #Could be fewer members than el

temp5 <- sapply(elm$ID1, function(x){which(vm$ID %in% x)})
temp6 <- sapply(elm$ID2, function(x){which(vm$ID %in% x)})

elm$t1 <- vm$Time[temp5]
elm$t2 <- vm$Time[temp6]
elm$tMax <- pmax(elm$t1,elm$t2)
elm$tDiff <- abs(elm$t1-elm$t2)
vlm <-unique(data.frame(ID = c(elm$ID1, elm$ID2), Time = c(elm$t1, elm$t2), stringsAsFactors=F))

vl <- subset(vl, ID%in%vm$ID)
el <- subset(el, ID1%in%vm$ID & ID2%in%vm$ID )

#Order both list elements by time point
g1 <- list(v=vl[order(vl$Time),], e=el[order(el$tMax),], f=el[order(el$tMax),])
g2 <- list(v=vlm[order(vlm$Time),], e=elm[order(elm$tMax),], f=elm[order(elm$tMax),])

#Filter out newest years for the sake of sample size
while(nrow(subset(g1$v,Time==max(Time)))<=63) {
  g1 <- tFilt(g1, max(g1$v$Time)-1)
}
 
#Ensure these two are comparable with eachother (equal number of cases per year) 
t1 <- table(g1$v$Time)
t2 <- table(g2$v$Time)

m1 <- as.numeric(max(names(t1[t1>100])))
m2 <- as.numeric(max(names(t2[t1>100])))

g1 <- tFilt(g1, m1)
g2 <- tFilt(g2, m2)

t1 <- table(g1$v$Time)
t2 <- table(g2$v$Time)

maxT1 <- as.numeric(tail(table(t1),1))
maxT2 <- as.numeric(tail(table(t2),1))
maxT <- max(maxT1, maxT2)
minT <- min(maxT1, maxT2)

maxG <- ifelse(maxT1>maxT2, g1, g2)
                                       
remV <- rev(maxG$v$ID)[1:(maxT-minT)]
maxG$v <- subset(maxG$v, !ID%in%remV)
maxG$e <- subset(maxG$e, !ID1%in%remV & !ID2%in%remV)

#
if(maxT1>maxT2) {
  g1 <- maxG
}else{
  g2 <- maxG
}

minE <- function(g) {
  #Permanently remove edges from the new year such that only the closest edge between new vertices and old vertices remains
  #Obtain new vertices and remove any internal edges within the new vertices 
  #We are not interested in completely new clustering
  nV <- subset(g$v, Time==max(Time))
  
  #The subset of vertices excluding those at the oldest year
  subV <- subset(g$v, Time>min(Time))
  
  #Obtain the closest retrospective edge of every vertex beyond the oldest year
  minE <- bind_rows(lapply(1:nrow(subV), function(i){
    v <- subV[i,]
    incE <- subset(g$e, (ID1%in%v$ID)|(ID2%in%v$ID))
    retE <- subset(incE, (tMax==v$Time)&(tDiff>0))
    retE[which(retE$Distance==min(retE$Distance))[[1]],]
  }))
  
  #Only closest retrospective edges are kept for edges from new cases. 
  g$e <- subset(g$e, tMax!=max(tMax))
  g$e <- rbind(g$e, subset(minE, tMax==max(tMax)))
  g$f <- subset(minE, tMax < max(tMax))
  
  return(g)
}

g1 <- minE(g1)
g2 <- minE(g2)

g1$v <- subset(g2$v, ID%in%g1$v$ID)
g1$e <- subset(g2$e, ID1%in%g1$v$ID & ID2%in%g1$v$ID)
g2$v <- subset(g2$v, ID%in%g1$v$ID)
g2$e <- subset(g2$e, ID1%in%g1$v$ID & ID2%in%g1$v$ID)

saveRDS(g1, "nomet_G.rds")
saveRDS(g2, "met_G.rds")
