#A process which interprets clmp cluster data such that it is comparable to tn93 cluster data for scoring the effectiveness of clustering from tn93 Output.
#USAGE: Rscript clmpAnalysis FastTreeOutput.nwk

require(clmp)

#Expecting an ML tree from an alignmnent of HIV sequence data
#Expecting patient information in the format ID_Date
args = commandArgs(trailingOnly = T)
t <- read.tree(args[1])

####- TO-DO: Modulate parameters of clump to produce different sets of cluster data
res <- clmp(t)

#Establish a set of node ids coupled with collection dates
temp <- sapply(res$tip.label, function(x) strsplit(x, '_')[[1]])
nodes <- temp
attr(nodes, "year") <- as.numeric(temp[2,])

#To make clmp cluster data comparable to tn93 cluster data
####- TO-DO: Pass to tn93 analysis software
clu <- list()
clu$membership <- res$clusters
clu$csize <- unname(table(res$clusters)) #Cluster number 0 is reserved for all singletons 
clu$no <- clu1$csize[[1]] + (length(clu1$csize) - 1)