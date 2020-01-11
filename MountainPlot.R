source("~/git/tn/comp_Lib.R")
require(R.utils)

## Generating Analysis
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
#The name/path of the output file, will both a pdf summary, a set of all clustering data, and a complete version of the graph in question 
runArgs <- commandArgs(trailingOnly=T, asValues=T, defaults = list(f="stdin",o=NA,t=1,g=""))
iFile <- runArgs$f
oFile <- ifelse(is.na(runArgs$o), gsub(".txt$", "", iFile), runArgs$o)
gFile <- runArgs$g

#Load or create a graph, saving a newly created graph in an accessible file for later use
if ((!is.na(gFile))&file.exists(gFile)) {
    g <- readRDS(gFile)
} else {
  g <- impTN93(iFile)
  saveRDS(g, file = paste0(oFile, "_G.rds"))
}

#Obtain cluster info for all subgraphs
res <- gaicRun(g)

#Save all growth data in accessable files
saveRDS(res, file = paste0(oFile, "_GD.rds"))
cutoffs <- as.numeric(names(res))

## Generate Pictures and output summary
#__________________________________________________________________________________________________________________________#

#Extract GAICs, AICs and cutoffs for graphing purposes
gaics <- sapply(res, function(x) {x$gaic})
modAIC <- sapply(res, function(x) {x$ageFit$aic})
nullAIC <- sapply(res, function(x) {x$nullFit$aic})

#Create visual output pdf
#pdf(file = paste0(oFile, "_VS.pdf"), width = 10, height = 12)
par(mfrow=c(2, 1), mar = c(1,4,1,2), oma=c(5,4,1,2), cex.lab=1.2)

#Plot both AIC measurements for comparison
plot(cutoffs, modAIC, ylab="AIC Comparison", xlab="", ylim=c(0,max(c(modAIC, nullAIC))))

#Background
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='blanchedalmond', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(h=0)

#Locate minimum
abline(v=cutoffs[which(gaics==min(gaics))[[1]]], lty=3, lwd=2)

polygon(c(0, cutoffs, max(cutoffs)), c(0, modAIC, 0) , col=rgb(1,0,0,0.4))
polygon(c(0, cutoffs, max(cutoffs)), c(0, nullAIC, 0) , col=rgb(0,1,1,0.4))
legend("topright", bg="white",
       legend=c("Proposed Model", "Null Model", "Overlap"), 
       fill=c(rgb(1,0,0,0.4), rgb(0,1,1,0.4), rgb(0,0.65,0.65,0.7)))


#Plot GAIC
plot(cutoffs, gaics, type = "n", ylim=c(min(gaics),max(gaics)), xlab="Cutoffs", ylab = "Difference (GAIC)")

#Background
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='blanchedalmond', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(h=0)

#Plot Lines/points
lines(cutoffs, gaics, lwd=1.6, col="orangered")
points(cutoffs, gaics)
abline(h=0)

#Locate minimum and annotate Graph
abline(v=cutoffs[which(gaics==min(gaics))[[1]]], lty=3,lwd=2)
text(cutoffs[which(gaics==min(gaics))[[1]]+2.0], min(gaics), labels= round(min(gaics)))
text(cutoffs[which(gaics==min(gaics))[[1]]+2.0], max(gaics)-0.10*(max(gaics)), labels=cutoffs[which(gaics==min(gaics))[[1]]])

#Title
par(xpd=NA)
title(xlab="Cutoffs")

dev.off()

#Save Optimal cutoff information
saveRDS(res[[which(gaics==min(gaics))]], file = paste0(oFile, "_Optimal.rds"))