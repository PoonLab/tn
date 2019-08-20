##TO-DO: Specify source-file Location
source("~/git/tn/newLib.R")

## USAGE: Rscript ~/git/tn/OpClusters.R.R __tn93output.txt__##
#Options...
# -f: The full file name (with path) to a tn93 output file, by default, this takes the standard input.
# -o: The output file without an extension (with path). Several files outputs will be made distinguished by letters.
# -y: A number for the purposes of filtering by year. Removes n years from the end of a data set (defaults to 0)
# -t: Threads - how many parallel processes will be run at once (defaults to 1).
# -m: Takes the name and path of a meta-data csv. Containing Age, sex, risk and Diagnostic Year (overwrites collection year)
# -r: How many repeats of 0.8, 0.6, and 0.4 resamples will be taken. (defaults to 20)

#EX1: runArgs <- list(f="~/Seattle/tn93St.txt", o=NA, y=0, t=8, m=NA, r=1)
#EX2: runArgs <- list(f="~/Seattle/tn93St.txt", o=NA, y=0, t=8, m=NA, r=30)

#Test-EX: runArgs <- list(f=NA, o="~/Data/Tennessee/analysis/tn93TnsubB_nomet", y=0, m=NA, g="~/Data/Tennessee/analysis/tn93TnsubB_nomet_G.rds", r=30)
#Test-Ex: runArgs <- list(f=NA, o="~/Data/Tennessee/analysis/tn93TnsubB_met", y=0, m=NA, g="~/Data/Tennessee/analysis/tn93TnsubB_met_G.rds", r=30)


## Generating Analysis
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
#The name/path of the output file, will both a pdf summary, a set of all clustering data, and a complete version of the graph in question 
runArgs <- commandArgs(trailingOnly=T, asValues=T, defaults = list(f="stdin",o=NA,y=0,t=1,m=NA,r=20))
iFile <- runArgs$f
oFile <- ifelse(is.na(runArgs$o), gsub(".txt$", "", infile), runArgs$o)
mtD <- runArgs$m
gFile <- runArgs$g
repeats <- runArgs$r

#Load or create a graph, saving a newly created graph in an accessible file for later use
if (!is.nan(gFile)) {
  g <- readRDS(gFile)
} else {
  g <- impTN93(iFile, mtD)
  saveRDS(g, file = paste0(oFile, "_G.rds"))
}

#Create a set of sub-sampled graphs
gs  <- lapply(rep(c(0.8, 0.6, 0.4), repeats), function(x){
  
  #Subsample a random set of n cases from the total graph
  iG <- g
  sID <- sample(iG$v$ID, size=round(x*nrow(iG$v)), replace=F)
  iG$v <- subset(iG$v, ID%in%sID)
  iG$e <- subset(iG$e, ID1%in%iG$v$ID & ID2%in%iG$v$ID)
  
  #Filter out newest years for the sake of sample size
  while(nrow(subset(iG$v,Time==max(Time)))<=63) {
    iG <- tFilt(iG, as.numeric(tail(names(table(iG$v$Time)),2))[[1]])
    iG <- clsFilt(iG)
  }
  
  #Save a copy of the complete list of minimum edges
  iG$f <- bpeFreq(iG)
  
  return(iG)
}) 

saveRDS(gs, file = paste0(oFile, "_GS.rds"))
#gs <- readRDS(paste0(oFile, "_GS.rds")))

#Run analysis on all subsampled graphs
runs <- lapply(gs, function(iG) {
  print(rev(iG$v$ID)[[1]])
  gaicRun(iG)
})

#Save all growth data in accessable files
saveRDS(runs, file = paste0(outfile, "_RD.rds"))
#runs <- readRDS(paste0(oFile, "_RD.rds")))

cutoffs <- as.numeric(names(runs[[1]])) 

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

#For defining range on a plot
minmin <- min(mins)
maxmax <- max(maxs)

#Create output pdf
pdf(file = paste0(outfile,"RVS.pdf"), width = 15, height = 10)

#Plot Generation
par(mfrow=c(1,2), cex.lab=1.2, main.font=1, cex.lab=1.2, font.main=1)

#Plot of GAICs over a wide range of robustness
plot.new()
title(xlab= "Cutoffs", 
      ylab = "GAIC")
plot.window(xlim = c(min(cutoffs), max(cutoffs)), ylim = c(minmin, maxmax))
axis(2, at=round(seq(minmin,maxmax,((maxmax-minmin)/10))), labels = round(seq(minmin,maxmax,((maxmax-minmin)/10))), las=2)
axis(1, at=seq(0,0.04,0.005), labels=seq(0,0.040,0.005))

# create a background
#bg <- par('usr')
#rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='linen', border=NA)
#abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
#abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
#abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
#abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
#abline(0,0, lty=2, lwd=3, col='grey50')
#box()

#Add lines and a smooth trend
for (i in gaics[seq(1,(repeats*3-2),3)]){lines(as.numeric(names(i)), unname(i), col=alpha("darkblue",0.4))}
for (i in gaics[seq(2,(repeats*3-1),3)]){lines(as.numeric(names(i)), unname(i), col=alpha("darkorchid",0.6))}
for (i in gaics[seq(3,(repeats*3),3)]){lines(as.numeric(names(i)), unname(i), col=alpha("darkturquoise",0.6))}
legend("topright", bg="white",
       legend=c("Resample 80% of cases", "Resample 60% of cases", "Resample 40% of cases"), 
       fill=c("darkblue", "darkorchid", "darkturquoise"), 
       title = paste0("Resamples Groups (N=",repeats,")"))

smooth <- smooth.spline(df)

#Add and make clear the minimum
smthmin <- predict(smooth)$x[predict(smooth)$y == min(predict(smooth)$y)]
lines(smooth, lwd=2)
abline(v=smthmin, lty=2)
axis(3, smthmin)
range <- c(quantile(minsLoc, 0.25),quantile(minsLoc, 0.75))
axis(3, at=range, pos=maxmax, labels=F, tcl=0.5)
axis(3, at=range, pos=maxmax, labels=F, tcl=-0.5)

#Density data generation
d <- density(minsLoc)
d1 <- density(minsLoc[seq(1,(repeats*3-2),3)], cut=4)
d2 <- density(minsLoc[seq(2,(repeats*3-1),3)])
d3 <- density(minsLoc[seq(3,(repeats*3),3)])

#Density plot generation
plot(d, ylim=c(0,max(d1$y)), col="white", xlab = "Cutoffs", main="Kernal Density of Optimal Cutoff (Bandwidth = 0.0007)")

# create a background
#bg <- par('usr')
#rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='linen', border=NA)
#abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
#abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
#abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
#abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
#abline(h=0, lty=2, lwd=3, col='grey50')
#box()

#Add data to plot
polygon(d2, col=alpha("darkorchid",1))
polygon(d3, col=alpha("darkturquoise",0.8))
polygon(d1,col=alpha("darkblue",0.6))

abline(v=smthmin, lty=2)

dev.off()