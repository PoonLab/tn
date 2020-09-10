library(dplyr,verbose = FALSE)
library(data.table)

#Obtain tn93 edge list from file
iFile <- "Data/gisaid/sc2tn93.csv"
reVars='/|\\|'
varInd=c(5,6,2)
varMan=NA
dateFormat="%Y-%m-%d"
partQ=0.95


idt <- fread(iFile)

#Reformat edge list as data table object with predictors extracted from sequence header
temp1 <- sapply(idt$ID1, function(x) (strsplit(x,'/|\\|')[[1]]))
temp2 <- sapply(idt$ID2, function(x) (strsplit(x,'/|\\|')[[1]]))
el <- data.table(ID1=as.character(temp1[5,]), t1=as.Date(temp1[6,]), l1=as.character(temp1[2,]),
                 ID2=as.character(temp2[5,]), t2=as.Date(temp2[6,]), l2=as.character(temp2[2,]),
                 Distance = as.numeric(idt$Distance), stringsAsFactors= F)


#Obtain the maximum time and time difference between the head and tail of each edge
el[,"tMax" := pmax(el$t1,el$t2)] 
el[,"tDiff" := (el$t1-el$t2)]

#Obtain list of unique sequences (also data table)
vl <- data.table(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), Location = c(el$l1, el$l2), stringsAsFactors=F)
vl <- vl[match(unique(vl$ID), vl$ID),]

#Vertices and edges lists together start a graph object
g <- list(v=vl[order(vl$Time),], e=el[order(el$tMax),])

#Get time info on the scale of months
g$v[, "Month" := as.numeric(round(julian(g$v$Time, origin = min(g$v$Time))/30))]

#Split into testing and training partitions
nV <- g$v[Month>=7]
subV <- g$v[Month<7]

#Remove edges from edge list if they are between two new vertices
g$e <- g$e[!((ID1%in%nV$ID) & (ID2%in%nV$ID))]

#Minimum retrospective edges (saved as "f" component of graph object)
start <- Sys.time()
temp <- g$e
temp[,1:3] <- g$e[,4:6]
temp[,4:6] <- g$e[,1:3]
temp$tDiff <- -(g$e$tDiff)
dt <- rbindlist(list(g$e,temp))

g$f <- dt[, list(Distance=min(.SD[tDiff<0]$Distance), 
            tDiff=.SD$tDiff[which.min(.SD[tDiff<0]$Distance)],
            ID2=.SD$ID2[which.min(.SD[tDiff<0]$Distance)]), .(ID1)]
Sys.time()-start

g$f




