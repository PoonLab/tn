##TO-DO: Specify source-file Location
source("~/git/tn/CluLib.R")

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

## Generating Analysis
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
#The name/path of the output file, will both a pdf summary, a set of all clustering data, and a complete version of the graph in question 
runArgs <- commandArgs(trailingOnly=T, asValues=T, defaults = list(f="stdin",o=NA,y=0,t=1,m=NA,r=20))
infile <- runArgs$f
outfile <- ifelse(is.na(runArgs$o), gsub(".txt$", "", infile), runArgs$o)
inputFilter <- as.numeric(runArgs$y)
threads <- as.numeric(runArgs$t)
metData <- runArgs$m
repeats <- runArgs$r

#Save all growth data in accessable files
g <- createGraph(infile, inputFilter, metData)
saveRDS(g, file = paste0(outfile, "G.rds"))

#Initialize a set of cutoffs to observe
steps <- head(hist(E(g)$Distance, plot=FALSE)$breaks,-5)
cutoffs <- seq(0 , max(steps), max(steps)/50)

#Create a list of runs and an index for naming purposes
runProps <- rep(c(0.8,0.6,0.4), repeats)

#Create Multiple Runs at various sub-samples
runs <- lapply(1:length(runProps), function(i) {
  
  #Progress Tracking
  runProp <- runProps[[i]]

  #Create a sample subgraph
  sampleV <- sample(V(g), round(length(V(g))*runProp))
  subG <- induced_subgraph(g,sampleV)
  
  #Create a set of subgraphs
  gs <- multiGraph(subG, cutoffs, threads)
  names(gs) <- cutoffs
                 
  #Generate Growth data
  res <- gaicRun(gs, cutoffs, threads)

  #Label data
  names(res) <- cutoffs
  
  cat(paste0("\r", "                 ", "      - Total Progress ", round(i/length(runProps)*100,1), "%")) 
  
  return(res)
})

#Save all growth data in accessable files
saveRDS(runs, file = paste0(outfile, "RD.rds"))

# runs <- readRDS("~/Seattle/tn93StRD.rds")

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
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='linen', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(0,0, lty=2, lwd=3, col='grey50')
box()

#Add lines and a smooth trend
for (i in gaics[seq(1,(repeats*3-2),3)]){lines(as.numeric(names(i)), unname(i), col=alpha("darkblue",0.4))}
for (i in gaics[seq(2,(repeats*3-1),3)]){lines(as.numeric(names(i)), unname(i), col=alpha("darkcyan",0.4))}
for (i in gaics[seq(3,(repeats*3),3)]){lines(as.numeric(names(i)), unname(i), col=alpha("cadetblue1",0.4))}
legend("topright", bg="white",
       legend=c("Resample 80% of cases", "Resample 60% of cases", "Resample 40% of cases"), 
       fill=c("darkblue", "darkcyan", "cadetblue1"), 
       title = paste0("Resamples Groups (N=",repeats,")"))

smooth <- smooth.spline(df)

#Add and make clear the absolute minimum of the trendline
smthmin <- predict(smooth)$x[predict(smooth)$y == min(predict(smooth)$y)]
lines(smooth, lwd=2)
abline(v=smthmin, lty=2)
axis(3, smthmin)
range <- c((smthmin+sd(minsLoc)),(smthmin-sd(minsLoc)))
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
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='linen', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(h=0, lty=2, lwd=3, col='grey50')
box()

#Add data to plot
polygon(d2, col=alpha("darkcyan",0.6))
polygon(d3, col=alpha("cadetblue1",0.6))
polygon(d1,col=alpha("darkblue",0.6))
abline(v=smthmin, lty=2)

dev.off()