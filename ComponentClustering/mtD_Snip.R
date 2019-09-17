#MetData Runtest
source("~/git/tn/comp_An.R")

iFile <- "~/Data/Tennessee/analysis_PRO/tn93TnsubB.txt"

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

if (T) {
  mFile <- "~/Data/Tennessee/sourceData/TnMetD/tnMD.csv"
  
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
}

vl <- subset(vl, ID%in%vm$ID)
el <- subset(el, ID1%in%vm$ID & ID2%in%vm$ID )

#Order both list elements by time point
g1 <- list(v=vl[order(vl$Time),], e=el[order(el$tMax),], f=el[order(el$tMax),])
g2 <- list(v=vlm[order(vlm$Time),], e=elm[order(elm$tMax),], f=elm[order(elm$tMax),])

#Filter out newest years for the sake of sample size
while(nrow(subset(g1$v,Time==max(Time)))<=63) {
  g1 <- tFilt(g1, max(g1$v$Time)-1)
}

g2$v <- subset(g2$v, ID%in%g1$v$ID)
g2$e <- subset(g2$e, ID1%in%g1$v$ID & ID2%in%g1$v$ID)


#Close Filter the overall graph at this point to save future time complexity
g1 <- clsFilt(g1)
g2 <- clsFilt(g2)

#Save a copy of the complete list of minimum edges
g1$f <- likData(g1)
g2$f <- likData(g2)

#TOFIX - Edges exist referencing Vertices that do not

saveRDS(g1, "tn93TnsubB_nomet_G.rds")
saveRDS(g2, "tn93TnsubB_met_G.rds")

g1NM <- readRDS("~/Data/Tennessee/analysis_PRO/tn93TnsubB_nomet_G.rds")
g2met <- readRDS("~/Data/Tennessee/analysis_PRO/tn93TnsubB_met_G.rds")
