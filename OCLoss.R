##TO-DO: Specify source-file Location
source("~/git/tn/CluLib.R")

## USAGE: Rscript ~/git/tn/OpClusters.R.R __tn93output.txt__##
#Options...
# -f: The full file name (with path) to a tn93 output file, by default, this takes the standard input.
# -o: The output file without an extension (with path). Several files outputs will be made distinguished by letters.
# -y: A number for the purposes of filtering by year. Removes n years from the end of a data set (defaults to 0)
# -t: Threads - how many parallel processes will be run at once (defaults to 1).
# -m: Takes the name and path of a meta-data csv. Containing Age, sex, risk and Diagnostic Year (overwrites collection year)
# -r: How many yeasrs beyond the most recent year will be examined

#EX: runArgs <- list(f="~/Seattle/tn93St.txt", o=NA, t=8, m=NA, r=4)
#EX2: runArgs <-list(f="~/Tennessee/tn93Tn.txt", o="~/Tennessee/tn93TnMet", t=8, m="~/Tennessee/TnMetD/tnMD.csv", r=4)

## Generating Analysis
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
#The name/path of the output file, will both a pdf summary, a set of all clustering data, and a complete version of the graph in question 
runArgs <- commandArgs(trailingOnly=T, asValues=T, defaults = list(f="stdin",o=NA,t=8,m=NA,r=4))
infile <- runArgs$f
outfile <- ifelse(is.na(runArgs$o), gsub(".txt$", "", infile), runArgs$o)
threads <- as.numeric(runArgs$t)
metData <- runArgs$m
filterRange <- 0:runArgs$r

#Create Multiple Runs at various longitudinal cuts (with different amounts of new years truncated)
runs <- lapply(1:length(filterRange), function(i) {
  
  #To cut n years off of the graph
  x <- filterRange[[i]]

  #Save all growth data in accessable files
  g <- createGraph(infile, x, metData)
  
  #Initialize a set of cutoffs to observe
  steps <- head(hist(E(g)$Distance, plot=FALSE)$breaks,-5)
  cutoffs <- seq(0 , max(steps), max(steps)/50)
  
  #Setting Global Parameters for future graphing of results
  if (x==0) {
    saveRDS(g, file = paste0(outfile, "G.rds"))
    assign("maxY", max(V(g)$year),  envir = .GlobalEnv)
    assign("minY", min(V(g)$year),  envir = .GlobalEnv)
    assign("step", max(cutoffs) / (length(cutoffs)-1),  envir = .GlobalEnv)
    assign("cutoffs",  cutoffs,  envir = .GlobalEnv)
  }
  
  #Create a set of subgraphs based off of differing cluster parameter
  gs <- multiGraph(g, cutoffs, threads)
  names(gs) <- cutoffs
  
  #Obtain cluster info for all subgraphs
  res <- gaicRun(gs, cutoffs, threads)
  
  #Label data
  names(res) <- cutoffs
  
  #Progress Tracking
  cat(paste0("\r", "                 ", "     - Total Progress ", round(i/length(filterRange)*100,1), "%")) 
  
  return(res)
})

#To elimnate potential zombie processes
wait()

#Save all growth data in accessable files
saveRDS(runs, file = paste0(outfile, "LD.rds"))

# runs <- readRDS("~/Seattle/tn93StLD.rds")
# runs <- readRDS("~/Tennessee/tn93TnMetLD.rds")

## Generate Pictures and output summary
#__________________________________________________________________________________________________________________________#

#Obtain a list of vectors of GAICs for each filtered run
gaics <- lapply(rev(runs), function(run){sapply(run, function(x) {x$gaic})})

#The step distance between cutoff points
step <- max(cutoffs) / (length(cutoffs)-1)

#The cutoff values which aquire the minimum GAIC. Also called the Minimum GAIC Estimator (MGAICE).
minsLoc <- sapply(gaics, function(x){step*(which(x==min(x))[[1]]-1)}) 
mins <- sapply(gaics, function(x){min(x)}) 
maxs <- sapply(gaics, function(x){max(x)})

#For defining range on a plot
minmin <- min(mins)
maxmax <- max(maxs)

#Create output pdf
pdf(file = paste0(outfile,"LVS.pdf"), width=6, height=12)

#Plot Generation
par(mfrow=c(length(filterRange), 1), mar = c(1,4,1,2), oma=c(5,4,1,2), cex.lab=1.2)

#Make multiple plots for each run of GAICs with minimum labelled
for (i in 1:length(gaics)) {
  
  #Initialize plot and background
  GAIC <- gaics[[i]]
  plot(cutoffs, GAIC, ylim = c(minmin+(0.2*minmin),maxmax), xlab="", ylab = "GAIC")
  bg <- par('usr')
  rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='linen', border=NA)
  abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
  abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
  abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
  abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
  
  #Plot GAIC
  lines(cutoffs, GAIC, lwd=1.6, col="blue")
  legend("bottomright", legend = paste0("Years ", minY, "-", (maxY-length(runs))+i), cex = 1,bg ="white")
  points(x=c(minsLoc[i]), y=c(mins[i]), cex=1)
  
  #Represents the location of the past run's MGAICE, Loss Ratio = minGAIC / Past minGAIC
  if (i>1){abline(v=minsLoc[i-1], lty=2)}
  
  #Draws an arrow to represent the follow through of the previous minimum
  if (i<length(filterRange)){arrows(minsLoc[i], mins[i], minsLoc[i], minmin+(0.2*minmin), length=0.05)}
  
  #To create the Cutoff label
  if (i==length(filterRange)){
    par(xpd=NA)
    title(xlab="Cutoffs")
  }
}
par(xpd=F)

dev.off()
