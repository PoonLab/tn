##TO-DO: Specify source-file Location
source("~/git/tn/CluLib.R")

#Expecting tn93 output as second param
## USAGE: Rscript ~/git/tn/OpClusters.R ___D.txt ##
#EX: runArgs <- c("~/Seattle/tn93St.txt", NA, "4", "F")

## Generating Analysis
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
#The name/path of the output file, will both a pdf summary, a set of all clustering data, and a complete version of the graph in question 
runArgs <- commandArgs(trailingOnly = T, defaults = c("stdin",NA,"4","F"))
infile <- runArgs[1]
outfile <- ifelse(exists(runArgs[2]), runArgs[2], infile)
filterRange <- seq(0,as.numeric(runArgs[3]),1)
home <- as.logical(runArgs[4])
metData <- runArgs[5]

#Create Multiple Runs at various longitudinal cuts
runs <- mclapply(filterRange, function(x) {
  
  #Save all growth data in accessable files
  g <- createGraph(infile, x, metData)
  saveRDS(g, file = paste0(outfile, "G.rds"))
  
  #Initialize a set of cutoffs to observe
  steps <- head(hist(E(g)$Distance, plot=FALSE)$breaks,-5)
  cutoffs <- seq(0 , max(steps), max(steps)/50)
  
  #Setting Parameters for future graphing of results
  if (x==0) {
    maxY <- max(V(g)$year)
    step <- max(cutoffs) / (length(cutoffs)-1)
  }
  if (x==max(filterRange)) {minY <- max(V(g)$year)}
  
  #Create a set of subgraphs based off of differing cluster parameter
  gs <- multiGraph(g)
  names(gs) <- cutoffs
  
  #Obtain cluster info for all subgraphs
  res <- gaicRun(gs)
  
  #Label data
  names(res) <- cutoffs
})

#Save all growth data in accessable files
saveRDS(runs, file = paste0(outfile, "GD.rds"))

## Generate Pictures and output summary
#__________________________________________________________________________________________________________________________#

#Obtain a list of vectors of GAICs for each filtered run
gaics <- lapply(rev(runs), function(run){sapply(run, function(x) {x$gaic})})

#The cutoff values which aquire the minimum GAIC. Also called the Minimum GAIC Estimator (MGAICE).
minsLoc <- sapply(gaics, function(x){step*(which(x==min(x))[[1]]-1)}) 

#Create output pdf
pdf(file = args[2])

#Plot Generation
par(mfrow=c(2, length(runs)/2))

for (i in 1:length(gaics)) {
  GAIC <- gaics[[i]]
  plot(Cutoffs, GAIC, main = paste0(minY, "-", (maxY-length(runs))+i))
  lines(Cutoffs, GAIC)
  abline(v=minsLoc[i], lty=2, lwd=1.5)
  
  #Represents the location of the past run's MGAICE, Loss Ratio = minGAIC / Past minGAIC
  if (i>1){
    abline(v=minsLoc[i-1], lty=2, col="orangered")
    lossRat <- GAIC[as.character(minsLoc[i-1])]/GAIC[as.character(minsLoc[i])]
    legend("bottomright", legend = paste0("Loss Ratio: ", round(lossRat,2)))
  }
}

dev.off()
