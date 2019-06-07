##TO-DO: Specify source-file Location
source("~/git/tn/CluLib.R")

#Expecting tn93 output as second param
## USAGE: Rscript ~/git/tn/OCRobust.R 
#EX: runArgs <- list(f="~/Seattle/tn93St.txt",o=NA,y=0,t=1,m=NA,r=30)

## Generating Analysis
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
#The name/path of the output file, will both a pdf summary, a set of all clustering data, and a complete version of the graph in question 
runArgs <- commandArgs(trailingOnly=T, asValues=T, defaults = list(f="stdin",o=NA,y=0,t=1,m=NA,r=20))
infile <- runArgs$f
outfile <- ifelse(is.na(runArgs$o), runArgs$f, infile)
inputFilter <- as.numeric(runArgs$y)
threads <- as.logical(runArgs$t)
metData <- runArgs$m
repeats <- runArgs$r

#Save all growth data in accessable files
g <- createGraph(infile, inputFilter, metData)
saveRDS(g, file = paste0(outfile, "G.rds"))

#Initialize a set of cutoffs to observe
steps <- head(hist(E(g)$Distance, plot=FALSE)$breaks,-5)
cutoffs <- seq(0 , max(steps), max(steps)/50)

#Create a list of runs and an index for naming purposes
runlist <- rep(c(0.8,0.6,0.4), repeats)

#Create Multiple Runs at various sub-samples
runs <- mclapply(1:length(runlist), function(i) {
  
  #Progress Tracking
  run <- runlist[[i]]
  cat(paste0("\r", "Running Analysis ", round(i/length(runlist)*100), "%")) 
  
  
  #Create a sample subgraph
  subG <- induced_subgraph(g,sample(V(g), round(length(V(g))*run)))
  
  #Create a set of subgraphs
  gs <- multiGraph(subG, cutoffs, threads)
  names(gs) <- cutoffs
                 
  #Generate Growth data
  res <- gaicRun(gs, cutoffs, threads, tracking=F)

  #Label data
  names(res) <- cutoffs
  
  return(res)
}, mc.cores = threads)

#Save all growth data in accessable files
saveRDS(runs, file = paste0(outfile, "RD.rds"))

## Generate Pictures and output
#__________________________________________________________________________________________________________________________#

#Obtain a list of vectors of GAICs for each filtered run
gaics <- lapply(runs, function(run){sapply(run, function(x) {x$gaic})})

#The stepdistance between cutoff points
step <- max(cutoffs) / (length(cutoffs)-1)

#The cutoff values which aquire the minimum GAIC. Also called the Minimum GAIC Estimator (MGAICE).
minsLoc <- sapply(gaics, function(x){step*(which(x==min(x))[[1]]-1)}) 

#Create a dataframe of cutoff to GAIC, and also fund the absolut minimumm and the absolute maximum for scale
df <- data.frame(Cutoff=as.numeric(names(unlist(gaics))), GAIC=unname(unlist(gaics)))
mins <- sapply(gaics, function(x){min(x)}) 
maxs <- sapply(gaics, function(x){max(x)})
minmin <- min(mins)
maxmax <- max(maxs)

#Create output pdf
pdf(file = paste0(outfile,"RVS.pdf"),width=20, height=10)

#Plot Generation
par(mfrow=c(1,2))

#Plot of GAICs over a wide range of robustness
plot.new()
title(xlab= "Cutoff values used to construct models and measure growth", 
      ylab = "GAIC: A null model's AIC subtracted from a proposed model AIC")
plot.window(xlim = c(min(cutoffs), max(cutoffs)), ylim = c(minmin, maxmax))
axis(2, at=round(seq(minmin,maxmax,((maxmax-minmin)/10))), labels = round(seq(minmin,maxmax,((maxmax-minmin)/10))), las=2, pos = 0)
axis(1, at=cutoffs, labels=cutoffs)
for (i in gaics){lines(as.numeric(names(i)), unname(i), col="grey")}
smooth <- smooth.spline(df)
smthmin <- predict(smooth)$x[predict(smooth)$y == min(predict(smooth)$y)]
lines(smooth, lwd=2)
abline(v=smthmin, lty=2)
axis(3, smthmin)
range <- c((smthmin+sd(minsLoc)),(smthmin-sd(minsLoc)))
axis(3, at=range, pos=maxmax, labels=F)

#Density Plot generation
d <- density(minsLoc)
d1 <- density(minsLoc[seq(1,(repeats*3-2),3)])
d2 <- density(minsLoc[seq(2,(repeats*3-1),3)])
d3 <- density(minsLoc[seq(3,(repeats*3),3)])

plot(d,ylim=(c(0,500)), col="white", main = "Kernal Density of MGAICE (Bandwidth = 0.0007)", xlab = "Cutoffs")
polygon(d2, col=alpha("orange",0.6))
polygon(d3, col=alpha("yellow",0.4))
polygon(d1,col=alpha("red",0.5))
abline(v=smthmin, lty=2)
legend("topleft", legend=c("Resample 40% of cases", "Resample 60% of cases", "Resample 80% of cases"), fill=c("red", "orange", "yellow"), title = paste0("Resamples Groups (N=",length(runs)/3))

dev.off()