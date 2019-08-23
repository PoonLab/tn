##TO-DO: Specify source-file Location
source("~/git/tn/newLib.R")

gNM <- readRDS("~/Data/Tennessee/analysis/tn93TnsubB_nomet_G.rds")
gMet <- readRDS("~/Data/Tennessee/analysis/tn93TnsubB_met_G.rds")

tComp <- tail(table(gMet$v$Time), length(table(gNM$v$Time))-1)

gMet$v <- subset(gMet$v, Time>=min(as.numeric(names(tComp))))
gMet$e <- subset(gMet$e, t1>=min(as.numeric(names(tComp))) & t2>=min(as.numeric(names(tComp))))
gMet$f <- bpeFreq(gMet)

gNM$v <- subset(gNM$v, Time>=min(as.numeric(names(tComp))))
gNM$e <- subset(gNM$e, t1>=min(as.numeric(names(tComp))) & t2>=min(as.numeric(names(tComp))))
gNM$f <- bpeFreq(gNM)

keepV <- vector()

for (t in as.numeric(names(tComp))) {
  tV <- subset(gNM$v, Time==t)
  keepV <- c(keepV,tV$ID[1:tComp[as.character(t)]])
}
gNM$v <- subset(gNM$v, ID%in%keepV)
gNM$e <- subset(gNM$e, ID1%in%keepV & ID2%in%keepV)
gNM$f <- bpeFreq(gNM)

saveRDS(gMet, file ="Met_RobComp_G")
saveRDS(gNM, file="NM_RobComp_G")

#For Smooths
####################################################################
smoothNM <- smooth
smoothMet <- smooth

minmin <- -10
maxmax <- 1

plot.new()

# create a background
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='linen', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
box()

title(xlab= "Cutoffs", 
      ylab = "Smoothed Trend in GAIC over 90 Runs")
plot.window(xlim = c(min(cutoffs), max(cutoffs)), ylim = c(minmin, maxmax))
axis(2, at=round(seq(minmin,maxmax,((maxmax-minmin)/10))), labels = round(seq(minmin,maxmax,((maxmax-minmin)/10))), las=2)
axis(1, at=seq(0,0.04,0.005), labels=seq(0,0.04,0.005))
lines(smoothNM, col="darkturquoise", lwd=1.6)
lines(smoothMet, col="orangered", lwd=1.6)
abline(h=0, lty=2)


