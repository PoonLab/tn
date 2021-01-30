source("~/git/tn/comp_Lib.R")

getageD <- function(subG) {
  tTab <- as.numeric(table(subG$v$Time))
  tDiffs <-(max(subG$f$tMax))-as.numeric(names(table(subG$f$tMax)))
  maxD <- max(subG$e$Distance)
  ageD <- bind_rows(lapply(tDiffs[tDiffs>0] , function(x){
    data.frame(tDiff=x, 
               Positive=nrow(subset(subG$f, tDiff==x & Distance<=maxD)), 
               vTotal=sum(tail(tTab,-x)))
  }))  
  
  df <- data.frame(tDiff=ageD$tDiff, Rate=ageD$Positive/ageD$vTotal)
  
  return(df)
}

sG <- readRDS("Data/ThesisData/sG.rds")
nG <- readRDS("Data/ThesisData/nG.rds")
nG_Dates <- readRDS("Data/ThesisData/nG_Dates.rds")
tG <- readRDS("Data/ThesisData/tG.rds")
tG_Diag <- readRDS("Data/ThesisData/tG_Diag.rds")

sD <- getageD(sG)
nD <- getageD(nG)
tD <- getageD(tG)
tD_Diag <- getageD(tG_Diag)

col1 <- "dodgerblue"
col2 <- "orange2"
col3 <- "orangered"
col4 <- "indianred1"
col5 <- "indianred4"

if(F){
  lm(log(sD$Rate)~tDiff, sD)
  lm(log(nD$Rate)~tDiff, nD)
  lm(log(tD$Rate)~tDiff, tD)
  tD_Diag_Ad <- tD_Diag
  tD_Diag_Ad$Rate <- tD_Diag$Rate+0.0000000000001
  lm(log(tD_Diag_Ad$Rate)~tDiff, tD_Diag)
  
}

#pdf(file="~/decayPlot-new.pdf", width=15, height=8)
par(mfrow=c(1,2),mar=c(5,6,4,1)+.1)

plot(x='n', y='n', xlab="Time lag (years)", ylab="Connection rate", xlim=c(1,12), ylim=c(0,0.5),
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


plot(x='n', y='n', xlab="Time lag (years)", ylab="Connection rate", xlim=c(1,30), ylim=c(0,0.3),
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
