source("comp_Lib.R")
require(scales)
require(optparse)

#Options...
# -f: The full file name (with path) to a tn93 output file, by default, this takes the standard input.
# -o: The output file without an extension (with path). Several files outputs will be made distinguished by letters.
# -y: A number for the purposes of filtering by year. Removes n years from the end of a data set (defaults to 0)
# -t: Threads - how many parallel processes will be run at once (defaults to 1).
# -m: Takes the name and path of a meta-data csv. Containing Age, sex, risk and Diagnostic Year (overwrites collection year)
# -r: How many repeats of 0.8, 0.6, and 0.4 resamples will be taken. 

## Generating Analysis
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
#The name/path of the output file, will both a pdf summary, a set of all clustering data, and a complete version of the graph in question 
option_list <- list( 
  make_option(c("-f", "--file"), default="stdin"),
  make_option(c("-o", "--output"), default=""),
  make_option(c("-g", "--graph"), default=""),
  make_option(c("-r", "--repeats"), default=20),
  make_option(c("-s", "--sampleVar"), default=F),
  make_option(c("-m", "--meta"), default=""))

opt <- parse_args(OptionParser(option_list=option_list))

iFile <- opt$f
oFile <- ifelse(opt$o%in%"", gsub(".txt$", "", iFile), opt$o)
mtD <- opt$m
gFile <- opt$g
repeats <- opt$r
varSamp <- opt$s

print(opt)

#Load or create a graph, saving a newly created graph in an accessible file for later use
if (file.exists(gFile)) {
  g <- readRDS(gFile)
} else {
  g <- impTN93(iFile, mtD)
  saveRDS(g, file = paste0(oFile, "_G.rds"))
}

#Assign a set of values to repeats
if(varSamp) {
  repeats <- rep(c(0.8,0.6,0.4), repeats)
} else {
  repeats <- rep(0.8, repeats*3)
}

#Create a set of sub-sampled graphs
gs  <- lapply(repeats, function(i){
  
  #Subsample a random set of n cases from the total graph
  iG <- g
  sID <- sample(iG$v$ID, size=round(i*nrow(iG$v)), replace=F)
  iG$v <- subset(iG$v, ID%in%sID)
  iG$e <- subset(iG$e, ID1%in%iG$v$ID & ID2%in%iG$v$ID)
  
  #Filter out newest years for the sake of sample size
  if(max(iG$v$Time)<max(g$v$Time)){iG <- clsFilt(iG)}
  
  #Save a copy of the complete list of minimum edges
  iG$f <- likData(iG)
  
  return(iG)
}) 

saveRDS(gs, file = paste0(oFile, "_GS.rds"))
#gs <- readRDS(paste0(oFile, "_GS.rds")))

#Run analysis on all subsampled graphs
cutoffs <- seq(0,0.04,0.0008) 
runs <- lapply(gs, function(iG) {
  print(rev(iG$v$ID)[[1]])
  gaicRun(iG,cutoffs)
})

#Save all growth data in accessable files
saveRDS(runs, file = paste0(oFile, "_RD.rds"))
#runs <- readRDS("~/Data/Tennessee/analysis/tn93TnsubB_met_RobComp_RD.rds")

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
pdf(file = paste0(oFile,"_RVS.pdf"), width = 15, height = 10)

#Plot Generation
par(cex.lab=1.2, main.font=1, cex.lab=1.2, font.main=1)
if (varSamp) {par(mfrow=c(1,2))}

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

if (varSamp) {
  #Add lines and a smooth trend
  for (i in gaics[which(repeats==0.8)]){lines(as.numeric(names(i)), unname(i), col=alpha("darkblue",0.4))}
  for (i in gaics[which(repeats==0.6)]){lines(as.numeric(names(i)), unname(i), col=alpha("darkorchid",0.6))}
  for (i in gaics[which(repeats==0.4)]){lines(as.numeric(names(i)), unname(i), col=alpha("darkturquoise",0.6))}
  legend("topright", bg="white",
         legend=c(paste0("Resample 80% (n=", length(repeats)/3, ")"), 
                  paste0("Resample 60% (n=", length(repeats)/3, ")"), 
                  paste0("Resample 40% (n=", length(repeats)/3, ")")), 
         fill=c("darkblue", "darkorchid", "darkturquoise"))
  
} else {
  for (i in gaics) {lines(as.numeric(names(i)), unname(i), col="darkturquoise")}
}

smooth <- smooth.spline(df)

#Add and make clear the minimum
smthmin <- predict(smooth)$x[predict(smooth)$y == min(predict(smooth)$y)]
lines(smooth, lwd=2)
abline(v=smthmin, lty=2)
axis(3, smthmin)
range <- c(quantile(minsLoc, 0.25),quantile(minsLoc, 0.75))
axis(3, at=range, pos=maxmax, labels=F, tcl=0.5)
axis(3, at=range, pos=maxmax, labels=F, tcl=-0.5)

if (varSamp) {
  
  #Density data generation
  d <- density(minsLoc)
  d1 <- density(minsLoc[which(repeats==0.8)], cut=4)
  d2 <- density(minsLoc[which(repeats==0.6)])
  d3 <- density(minsLoc[which(repeats==0.4)])
  
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
  
}


dev.off()