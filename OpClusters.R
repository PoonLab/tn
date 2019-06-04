##TO-DO: Specify source-file Location
source("~/git/tn/CluLib.R")

#Expecting tn93 output as second param
## USAGE: Rscript ~/git/tn/OpClusters.R ___D.txt ##
#EX: runArgs <- c("~/Seattle/tn93St.txt", NA, "0", "F")

## Generating Analysis
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
#The name/path of the output file, will both a pdf summary, a set of all clustering data, and a complete version of the graph in question 
runArgs <- commandArgs(trailingOnly = T, defaults = c("stdin",NA,0,F))
infile <- runArgs[1]
outfile <- ifelse(exists(runArgs[2]), runArgs[2], infile)
inputFilter <- as.numeric(runArgs[3])
home <- as.logical(runArgs[4])

#Save all growth data in accessable files
g <- createGraph(infile, inputFilter)
saveRDS(g, file = paste0(outfile, "G.rds"))

#Initialize a set of cutoffs to observe
steps <- head(hist(E(g)$Distance, plot=FALSE)$breaks,-5)
cutoffs <- seq(0 , max(steps), max(steps)/50)

#Create a set of subgraphs
gs <- multiGraph(g)
names(gs) <- cutoffs

if (home) {
  #Generate Growth data
  res <- lapply(cutoffs, function(d) {
    cat(paste0("\r", "Running Analysis ", d/max(cutoffs)*100, "%"))
    
    #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
    subG <- gs[[as.character(d)]]
    
    #Obtain cluster information for this subgraph
    clusterAnalyze(subG)
    
  })
} else {
  #Generate Growth data
  res <- mclapply(cutoffs, function(d) {
    cat(paste0("\r", "Running Analysis ", d/max(cutoffs)*100, "%"))
    
    #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
    subG <- gs[[as.character(d)]]
    
    #Obtain cluster information for this subgraph
    clusterAnalyze(subG)
    
  }, mc.cores=8)
}

#Label data
names(res) <- cutoffs

#Save all growth data in accessable files
saveRDS(res, file = paste0(outfile, "GD.rds"))

## Generate Pictures and output summary
#__________________________________________________________________________________________________________________________#

#Obtain Minimum GAIC estemating cutoff threshold and the network associated with it
gaics <- sapply(res, function(x) {x$gaic})
do <- names(which(gaics==min(gaics))[1])
opt <- gs[[do]]

#Plot option ignores clusters of size 1 and provides a graph (for ease of overview, not for calculations)
optPG <- subgraph.edges(opt, E(opt), delete.vertices = T)

#Create visual output pdf
pdf(file = paste0(outfile, "VS.pdf"))

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

cat(paste0("\n","Done" ))