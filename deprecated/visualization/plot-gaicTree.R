sGs <- readRDS("Data/ThesisData/sTs.rds")
nGs <- readRDS("Data/ThesisData/nTs.rds")
tGs <- readRDS("Data/ThesisData/tTs.rds")
tGs_Comp <- readRDS("Data/ThesisData/tTs_Comp.rds")
tGs_Diag <- readRDS("Data/ThesisData/tTs_Diag.rds")

col1 <- "dodgerblue"
col2 <- "orange2"
col3 <- "orangered"
col4 <- "indianred1"
col5 <- "indianred4"

#Number of Clusters
getGAIC <- function(iGs) {sapply(iGs, function(x){x[[1]]$gaic})}

cutoffs <- seq(0,0.15,0.003)

pdf(file="~/gaicPlots-Trees.pdf", width=14, height=7)
par(mfrow=c(1,2),mar=c(5,6,4,1)+.1)

#PlotA
plot(x=0, y=0, type='n', xlab="Patristic distance threshold", ylab="AIC Difference (Proposed - Null)", xlim=c(0,0.15), ylim=c(-50, 20),
     cex.lab=1.8, cex.axis = 1.5)

sD <- getGAIC(sGs)
nD <- getGAIC(nGs)
tD <- getGAIC(tGs)


# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(10, -50, -20), col='white', lwd=3, lend=2)
abline(h=seq(20, -40, -20), col='white', lend=3)
abline(v=seq(0.05, 0.15, 0.05), col='white', lwd=3, lend=2)
abline(v=seq(0, 0.10, 0.05), col='white', lend=3)
abline(h=0, col="black", lwd=2)

points(cutoffs,sD, col=col1, cex=0.5)
lines(cutoffs,sD, lwd=1.5, col=col1)
abline(v=cutoffs[which(sD==min(sD))], col="black", lty=2, lend=3)
text(x=cutoffs[which(sD==min(sD))+2], y=round(min(sD)), labels = cutoffs[which(sD==min(sD))])

points(cutoffs,nD, col=col2, cex=0.5)
lines(cutoffs,nD, lwd=1.5, col=col2)
abline(v=cutoffs[which(nD==min(nD))], col="black", lty=2, lend=3)
text(x=cutoffs[which(nD==min(nD))+2], y=round(min(nD)), labels = cutoffs[which(nD==min(nD))])

points(cutoffs,tD, col=col3, cex=0.5)
lines(cutoffs,tD, lwd=1.5, col=col3)
abline(v=cutoffs[which(tD==min(tD))], col="black", lty=2, lend=3)
text(x=cutoffs[which(tD==min(tD))+2], y=round(min(tD)), labels = cutoffs[which(tD==min(tD))])

legend('topright', legend=c('Seattle', 'Alberta','Tennessee'),
       fill=c(col1, col2,col3), bty='n', cex=1.2)

#PlotB
plot(x=0, y=0, type='n', xlab="Patristic distance threshold", ylab="AIC Difference (Proposed - Null)", xlim=c(0,0.15), ylim=c(-50, 20),
     cex.lab=1.8, cex.axis = 1.5)

tD_Comp <- getGAIC(tGs_Comp)
tD_Diag <- getGAIC(tGs_Diag)



# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(10, -50, -20), col='white', lwd=3, lend=2)
abline(h=seq(20, -40, -20), col='white', lend=3)
abline(v=seq(0.05, 0.15, 0.05), col='white', lwd=3, lend=2)
abline(v=seq(0, 0.10, 0.05), col='white', lend=3)
abline(h=0, col="black", lwd=2)


points(cutoffs,tD_Comp, col=col4, cex=0.5)
lines(cutoffs,tD_Comp, lwd=1.5, col=col4)
abline(v=cutoffs[which(tD_Comp==min(tD_Comp))], col="black", lty=2, lend=3)
text(x=cutoffs[which(tD_Comp==min(tD_Comp))+2], y=round(min(tD_Comp)), labels = cutoffs[which(tD_Comp==min(tD_Comp))])
points(cutoffs,tD_Diag, col=col5, cex=0.5)
lines(cutoffs,tD_Diag, lwd=1.5, col=col5)
abline(v=cutoffs[which(tD_Diag==min(tD_Diag))], col="black", lty=2, lend=3)
text(x=cutoffs[which(tD_Diag==min(tD_Diag))+2], y=round(min(tD_Diag)), labels = cutoffs[which(tD_Diag==min(tD_Diag))])

legend('topright', legend=c('Tennessee (Collection Year)', 'Tennessee (Diagnostic Year)'),
       fill=c(col4, col5), bty='n', cex=1.2)

dev.off()