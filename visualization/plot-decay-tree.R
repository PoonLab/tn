source("~/git/tn/subT_Lib.R")

getageD <- function(subT) {
  
  #Check positives now that we have params
  subT$f$Positive <- (subT$f$xDist<=10000) #&(subT$f$BootStrap>=minB)
  
  #Layout of Time and time lag information
  tTab <- table(subT$v$Time)
  tDiffs <- sort(unique(abs(as.vector(subT$e$tDiff))))
  tdTab <- rep(nrow(subT$v)-1, length(tDiffs))
  names(tdTab) <- tDiffs
  largeTD <- tDiffs[which(tDiffs>(ceiling(max(tDiffs/2))))]
  leftOut <- sapply(1:length(largeTD), function(i){
    td <- largeTD[[i]]
    ts <- as.numeric(names(tTab))[between(as.numeric(names(tTab)), (max(subT$v$Time)-td+1), (min(subT$v$Time)+td-1))]
    sum(tTab[as.character(unique(ts))]) 
  })
  
  #Obtains the "attempts". Or how many tips find a given tDiff possible
  #For example, it may be impossible for central time points to see the largest time difference in the set
  tdTab[as.character(largeTD)] <- tdTab[as.character(largeTD)]-leftOut
  
  #Sort Data into Age Data
  names(tdTab) <- tDiffs
  posTab <- rep(0,length(tDiffs))
  names(posTab) <- tDiffs
  tempTab <- table(round(subset(subT$f, Positive)$tDiff))
  posTab[names(tempTab)] <- as.numeric(tempTab)
  ageD <- data.frame(tDiff=tDiffs, Positives=as.numeric(posTab), 
                     Total=as.numeric(tdTab))
  
  df <- data.frame(tDiff=ageD[-1,]$tDiff, Rate=ageD[-1,]$Positives/ageD[-1,]$Total) 
  return(df)
}


sT <- readRDS("Data/ThesisData/sT.rds")
nT <- readRDS("Data/ThesisData/nT.rds")
tT <- readRDS("Data/ThesisData/tT.rds")
tT_Diag <- readRDS("Data/ThesisData/tT_Diag.rds")


col1 <- "dodgerblue"
col2 <- "orange2"
col3 <- "orangered"
col4 <- "indianred1"
col5 <- "indianred4"

sD <- getageD(sT)
nD <- getageD(nT)
tD <- getageD(tT)
tD_Diag <- getageD(tT_Diag)

#pdf(file="~/decayPlot-tree.pdf", width=15, height=8)
par(mfrow=c(1,2),mar=c(5,6,4,1)+.1)


plot(x='n', y='n', xlab="Time lag (years)", ylab="Connection rate", xlim=c(1,13), ylim=c(0,0.3),
     cex.lab=1.8, cex.axis = 1.5)

# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(0.1, 0.5, 0.1), col='white', lwd=3, lend=2)
abline(h=seq(0.05, 0.45, 0.1), col='white', lend=3)

points(sD, col=col1, lwd=1.5)
lines(sD, col=col1, lwd=1.5)
points(nD, col=col2, lwd=1.5)
lines(nD, col=col2, lwd=1.5)
points(tD, col=col3, lwd=1.5)
lines(tD, col=col3, lwd=1.5)

legend('topright', legend=c('Seattle', 'Alberta', 'Tennessee'),
       fill=c(col1, col2, col3), bty='n', cex=1.8)


plot(x='n', y='n', xlab="Time lag (years)", ylab="Connection rate", xlim=c(1,33), ylim=c(0,0.2),
     cex.lab=1.8, cex.axis = 1.5)

# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(0.1, 0.5, 0.1), col='white', lwd=3, lend=2)
abline(h=seq(0.05, 0.45, 0.1), col='white', lend=3)

points(tD, col=col4, lwd=1.5)
lines(tD, col=col4, lwd=1.5)
points(tD_Diag, col=col5, lwd=1.5)
lines(tD_Diag, col=col5, lwd=1.5)

legend('topright', legend=c('Tennessee (Collection Year)', 'Tennessee (Diagnostic Year)'),
       fill=c(col4, col5), bty='n', cex=1.8)

dev.off()
