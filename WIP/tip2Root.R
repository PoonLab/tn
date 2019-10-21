library(dplyr,verbose = FALSE)
library(ape)

iFile <- "~/Data/Seattle/analysis_PRO/FTStsubB.nwk"

#Creating an ape tree object from the newick file
t <- read.tree(iFile)

#Obtain lists of sequence ID and Time
temp <- sapply(t$tip.label, function(x) strsplit(x, '_')[[1]])
t$ID <- temp[1,]
t$Time  <- as.numeric(temp[2,])

inds <-  1:length(t$ID)
res <- lapply(inds, function(a) {
  ta <- t$Time[a]
  ida <- t$ID[a]
  
  bs <- lapply(inds[inds>a], function(b){
    tb <- t$Time[b]
    idb <- t$ID[b]
    p <- nodepath(t, a,b)
    p <- tail(p, -1)
    p <- head(p, -2)
    dist <- sum(t$edge.length[which(t$edge[,2]%in%p)])
    return(list(idb, tb, dist))
  })
  
  return(list(ta, ida, bs))
})

saveRDS(res, "~/distances.rds")

#If the newest timepoint contains a small number of sequences, we remove that timepoint from consideration
while(length(t$Time[t$Time==max(t$Time)]) <= minNS) {
  t <- tFilt(t, max(t$Time)-1)
}

t$tip.label


el <- data.frame(ID1=as.character(NULL), t1=as.numeric(NULL), 
                 ID2=as.character(NULL), t2=as.numeric(NULL), 
                 Distance = as.numeric(NULL), stringsAsFactors= F)

ps <- 1:length(t$ID) #Children
ps <- t$edge[which(t$edge[,2]%in%ps),1] #To parent nodes
ds <- rep(0, length(subTs)) #Distance travelled (init at 0)

tab <- table(ps)
cN <- as.numeric(names(which(tab>1)))

opairs <- matrix(nrow=2, ncol=1)
npairs <- lapply(cN, function(x) {
  inds <- combn(which(ps==x),2)
  return(inds)
})

opairs <- cbind(opairs, npairs)


while(length(ps>1)){
  ind <- sapply(ps, function(x){which(t$edge[,2]==x)})
  ps <- t$edge[ind,1]
  ds <- ds+t$edge.length[ind]
  length(levels(as.factor(ps)))
  length(ps)
}

l <- c(1,2,3,4,5,6,7,8,9,10)
