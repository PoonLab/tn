#Import Libraries
library(igraph,verbose = FALSE)

## args <- c("~/Seattle/Robust/RRout", "~/Seattle/stGraph.rds", _____)

#Input args, expecting a file of rubustness estimates (many resampled runs)
#Also takes a Graph of the total data and an output file name as arg3
args <- commandArgs(trailingOnly = T)

#Reads these resamples into data files
runs <- lapply(list.files(args[1]), function(x) {readRDS(file=paste0(args[1], "/", x))})
g <- readRDS(file=args[2])

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
pdf(file = paste0(args[3]))


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