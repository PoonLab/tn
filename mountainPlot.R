GD_bad <- readRDS("~/Data/Seattle/analysis_PRO/tn93StsubB_BAD_GD.rds") 
GD <- readRDS("~/Data/Seattle/analysis_PRO/tn93StsubB_GD.rds") 

gaics_bad <- sapply(GD_bad, function(x) {x$gaic})
gaics <- sapply(GD, function(x) {x$gaic})

mod_bad <- sapply(GD_bad, function(x) {summary(x$ageFit)$aic})
mod <- sapply(GD, function(x) {summary(x$ageFit)$aic})

cutoffs <- seq(0, 0.04, 0.0008)

plot(cutoffs, mod, ylab="GAIC is measured near the bottom)", xlab="Tuning", ylim=c(-80,400))

bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='blanchedalmond', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(h=0)

lines(mod_bad, lwd=5, lty=2, col="azure2")


x <- cutoffs
y1 <- mod
y2 <- mod_bad

polygon(c(x,rev(x)), c(y1, rev(y2)) , col="darkslategray4")
lines(cutoffs, mod-mod_bad, lwd=3, col="darkslategray4")
points(cutoffs, mod-mod_bad, pch=21, bg="black", col="white")

abline(v= which(abs(mod-mod_bad)== max(abs(mod-mod_bad)))*0.0008 , lty=2, lwd=2)
text(x= which(abs(mod-mod_bad)==max(abs(mod-mod_bad)))*0.0008, y=min(mod-mod_bad), cex=1.2, adj=0, 
     labels= c(paste0(round(min(mod-mod_bad)), " at ", which(abs(mod-mod_bad)==max(abs(mod-mod_bad)))*0.0008*100, "%"))) 
     
