#Import Libraries
library(igraph,verbose = FALSE)
library(dplyr,verbose = FALSE)
library(parallel,verbose = FALSE)
library(ggplot2,verbose = FALSE)
source("~/git/tn/CluLib.R")

#Expecting tn93 output as second param
## USAGE: Rscript ~/git/tn/OpClusters.R ___D.txt ##
# runArgs <- c("~/Seattle/tn93St.txt", NA, "0")

## Generating Analysis
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
runArgs <- commandArgs(trailingOnly = T, defaults = c("stdin",NA,0))
gs <- createGraphSet(runArgs)
cutoffs <- as.numeric(names(gs))

# Generate Growth data
res <- mclapply(cutoffs, function(d) {
  cat(paste0("\r", "Running Analysis ", d/max(cutoffs)*100, "%"))
  
  #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
  subG <- gs[[as.character(d)]]
  
  #Obtain cluster information for this subgraph
  clusterAnalyze(subG)
  
}, mc.cores=8)

#Label data
names(res) <- cutoffs

## Generate Pictures and output
#__________________________________________________________________________________________________________________________#

#Obtain Minimum GAIC estemating cutoff threshold and the network associated with it
gaics <- sapply(res, function(x) {x$gaic})
do <- names(which(gaics==min(gaics))[1])
opt <- gs[[do]]

#Plot option ignores clusters of size 1 and provides a graph (for ease of overview, not for calculations)
optPG <- subgraph.edges(opt, E(opt), delete.vertices = T)

#Create output pdf
pdf(file = paste0(gsub("\\..*", "", runArgs[2]), "VS.pdf"))

#Plot GAIC
gaicPlot(res)

#Plot Network
plot(optPG, vertex.size = 2, vertex.label = NA, vertex.color= "orange",
     edge.width = 0.65, edge.color = 'black', 
     margin = c(0,0,0,0))

dev.off()

#Obtain the information from opt cluster and print it to stOut
optClu <- components(opt)
optClu$years <- table(V(g)$year)
optClu$no <- NULL
optClu$csize <- sort(table(optClu$membership)[table(optClu$membership)>1], decreasing =T)
print(optClu)

#Save all growth data in accessable files
saveRDS(res, file = paste0(gsub("\\..*", "", runArgs[2]), "GD.rds"))

cat(paste0("\n","Done" ))