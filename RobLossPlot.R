#Expecting a folder containing multiple .RDS files of the same data set, longitudinally truncated at different years
## USAGE: Rscript ~/git/tn/robLossPlot.R ______ _______ ##

## args <- c("~/Seattle/Robust/OCLout", "~/Seattle/Robust/OCLVS.pdf", MinimumYear, MaximumYear )

#From arguments, extract the range of years, the range of cutoffs and the set of Clustering Outputs after filtering  
args <- commandArgs(trailingOnly = T)
runs <- lapply(list.files(args[1]), function(x) {readRDS(file=paste0(args[1], "/", x))})
minY <- args[3]
maxY <- args[4]
Cutoffs <- as.numeric(names(gaics[[1]]))
step <- max(Cutoffs) / (length(Cutoffs)-1)

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