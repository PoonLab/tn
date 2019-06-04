##TO-DO: Specify source-file Location
source("~/git/tn/CluLib.R")

#Expecting tn93 output as second param
## USAGE: Rscript ~/git/tn/OpClusters.R ___D.txt ##
#EX: runArgs <- list("~/Seattle/tn93St.txt", NA, "0", "20")

## Generating Analysis
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
#The name/path of the output file, will both a pdf summary, a set of all clustering data, and a complete version of the graph in question 
runArgs <- commandArgs(trailingOnly = T, defaults = c("stdin",NA,0,F,20))
infile <- runArgs[1]
outfile <- ifelse(exists(runArgs[2]), runArgs[2], infile)
inputFilter <- as.numeric(runArgs[3])
repeats <- as.numeric(runArgs[4])

#Save all growth data in accessable files
g <- createGraph(infile, inputFilter)
saveRDS(g, file = paste0(outfile, "G.rds"))

#Initialize a set of cutoffs to observe
steps <- head(hist(E(g)$Distance, plot=FALSE)$breaks,-5)
cutoffs <- seq(0 , max(steps), max(steps)/50)

#Create a list of runs and an index for naming purposes
runlist <- rep(c(0.8,0.6,0.4), repeats)

#Create Multiple Runs at various sub-samples
runs <- mclapply(runlist, function(run) {
  
  #Create a set of subgraphs
  gs <- multiGraph(g)
  names(gs) <- cutoffs
                 
  #Generate Growth data
  res <- mclapply(cutoffs, function(d) {
    cat(paste0("\r", "Running Analysis ", d/max(cutoffs)*100, "%"))
    
    #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
    subG <- gs[[as.character(d)]]
    
    #Obtain cluster information for this subgraph
    clusterAnalyze(subG)
    
  }, mc.cores=8)
  
  #Label data
  names(res) <- cutoffs
}, mc.cores=8)

#Save all growth data in accessable files
saveRDS(runs, file = paste0(outfile, "RD.rds"))

## Generate Pictures and output
#__________________________________________________________________________________________________________________________#

#Obtain a list of vectors of GAICs for each filtered run
gaics <- lapply(rev(runs), function(run){sapply(run, function(x) {x$gaic})})

#The cutoff values which aquire the minimum GAIC. Also called the Minimum GAIC Estimator (MGAICE).
minsLoc <- sapply(gaics, function(x){step*(which(x==min(x))[[1]]-1)}) 


#Create a dataframe of cutoff to GAIC, and also fund the absolut minimumm and the absolute maximum for scale
df <- data.frame(Cutoff=as.numeric(names(unlist(gaics))), GAIC=unname(unlist(gaics)))
mins <- sapply(gaics, function(x){min(x)}) 
maxs <- sapply(gaics, function(x){max(x)})
minmin <- min(mins)
maxmax <- max(maxs)
Cutoffs <- as.numeric(names(gaics[[1]]))
step <- max(Cutoffs) / (length(Cutoffs)-1)

#Create output pdf
pdf(file = paste0(outfile,"RVS.pdf"))

#Plot Generation
par(mfrow=c(1,2))

## TO-DO: Work on Scope here (at of axis)
plot.new()
title(xlab= "Cutoff values used to construct models and measure growth", 
      ylab = "GAIC: A null model's AIC subtracted from a proposed model AIC")
plot.window(xlim = c(min(Cutoffs), max(Cutoffs)), ylim = c(minmin, maxmax))
axis(2, at= seq(-210,50,10), labels = seq(-210, 50, 10), las=2, pos = 0)
axis(1, at=cutoffs, labels=Cutoffs)
points(df, col="grey")
smooth <- smooth.spline(df)
smthmin <- predict(smooth)$x[predict(smooth)$y == min(predict(smooth)$y)]
lines(smooth, lwd=2)
abline(v=smthmin, lty=2)
axis(3, smthmin)
range <- c((smthmin+sd(minsLoc)),(smthmin-sd(minsLoc)))
axis(3, at=range, pos=20, labels=F)

#Density Plot generation
d <- density(minsLoc)
d1 <- density(minsLoc[1:10])
d2 <- density(minsLoc[11:30])
d3 <- density(minsLoc[21:20])

plot(d,ylim=(c(0,500)), col="white", main = "Kernal Density of MGAICE (Bandwidth = 0.0007)", xlab = "Cutoffs")
polygon(d2, col=alpha("orange",0.6))
polygon(d3, col=alpha("yellow",0.4))
polygon(d1,col=alpha("red",0.5))
abline(v=smthmin, lty=2)
legend("topleft", legend=c("Resample 40% of cases", "Resample 60% of cases", "Resample 80% of cases"), fill=c("red", "orange", "yellow"), title = paste0("Resamples Groups (N=",length(runs)/3))

dev.off()