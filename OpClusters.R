##TO-DO: Specify source-file Location
source("~/git/tn/newLib.R")

## USAGE: Rscript ~/git/tn/OpClusters.R.R __tn93output.txt__##
#Options...
# -f: The full file name (with path) to a tn93 output file, by default, this takes the standard input.
# -o: The output file without an extension (with path). Several files outputs will be made distinguished by letters.
# -m: Takes the name and path of a meta-data csv. Containing Age, sex, risk and Diagnostic date (overwrites collection date)
# -g: The file path to saved graphical info. If one has already saved a graph, this will save time making one

#EX: runArgs <- list(f="~/Data/Seattle/tn93StsubB.txt", o=NA, y=0, m=NA, g="~/Data/Seattle/analysis/tn93StsubB_G.rds")

## Generating Analysis
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
#The name/path of the output file, will both a pdf summary, a set of all clustering data, and a complete version of the graph in question 
runArgs <- commandArgs(trailingOnly=T, asValues=T, defaults = list(f="stdin",o=NA,t=1,m=NA,g=NA))
iFile <- runArgs$f
oFile <- ifelse(is.na(runArgs$o), gsub(".txt$", "", iFile), runArgs$o)
mtD <- runArgs$m
gFile <- runArgs$g

#Load or create a graph, saving a newly created graph in an accessible file for later use
if (!is.nan(gFile)) {
  g <- readRDS(gFile)
} else {
  g <- impTN93(iFile, mtD)
  saveRDS(g, file = paste0(oFile, "_G.rds"))
}

#Obtain cluster info for all subgraphs
res <- gaicRun(g)

#Save all growth data in accessable files
saveRDS(res, file = paste0(oFile, "_GD.rds"))

## Generate Pictures and output summary
#__________________________________________________________________________________________________________________________#

#Extract GAICs and cutoffs for graphing purposes
gaics <- sapply(res, function(x) {x$gaic})
cutoffs <- names(res)  

#Create visual output pdf
pdf(file = paste0(oFile, "_VS.pdf"))

#Plot GAIC
plot(cutoffs, gaics, type = "n", ylim=c(min(gaics),max(gaics)), xlab="Cutoffs", ylab = "GAIC")
lines(cutoffs, gaics, lwd=1.6, col="orangered")
points(cutoffs, gaics)
abline(h=0)
abline(v=cutoffs[which(gaics==min(gaics))[[1]]], lty=3)
text(cutoffs[which(gaics==min(gaics))[[1]]+1.5], min(gaics), labels= round(min(gaics)))
text(cutoffs[which(gaics==min(gaics))[[1]]], max(c(gaics,gaics))-1.5, labels=cutoffs[which(gaics==min(gaics))[[1]]])

dev.off()