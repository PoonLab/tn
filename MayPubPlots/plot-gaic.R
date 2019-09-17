# contains the cluster information for graphs at different distance cutoffs
GDst <- readRDS("Data/Paper1/tn93StsubB_GD.rds")
GDna <- readRDS("Data/Paper1/tn93NAsubB_GD.rds")

st <- sapply(GDst, function(x) {x$gaic})
na <- sapply(GDna, function(x) {x$gaic})
x <- seq(0,0.04,0.0008)

pdf(file="Data/Paper1/gaic_tn.pdf", width=7.5, height=7.5)

par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(x, st, xlim=c(0, 0.04), xlab='TN93 distance cutoff',
     ylab='Generalized AIC', type='n', ylim=c(-60,max(st)), cex.lab=1.2)

# create a background
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='ivory2', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(h=0, lty=2, lwd=3, col='grey50')
box()

if(tn){
  GDtn <- readRDS("Data/Paper1/tn93TnsubB_met_GD.rds")
  tn <- sapply(GDtn, function(x) {x$gaic})
  
  lines(x, tn, col='orangered', lwd=2)
  points(x, tn, pch=24, col='white', bg='orangered', cex=1.1)
  
  text(x=x[which.min(tn)]-0.001, y=min(tn)-2, label=x[which.min(tn)], cex=0.8)
  
  legend(x=0.025, y=-30, legend=c('Seattle', 'N.Alberta', "Tennessee"), pch=c(21,22,24), 
         col='white', pt.bg=c('dodgerblue', 'orange2', "orangered"), pt.cex=1.5,
         bty='n')
  
}

lines(x, st, col='dodgerblue', lwd=2)
points(x, st, pch=21, col='white', bg='dodgerblue', cex=1.1)

text(x=x[which.min(st)]+0.0014, y=min(st), label=x[which.min(st)], cex=0.8)

lines(x, na, col='orange2', lwd=2)
points(x, na, pch=22, col='white', bg='orange2', cex=1.1)

text(x=x[which.min(na)]-0.001, y=min(na)-2, label=x[which.min(na)], cex=0.8)


if (!tn) {
  legend(x=0.025, y=-30, legend=c('Seattle', 'N.Alberta'), pch=c(21,22), 
         col='white', pt.bg=c('dodgerblue', 'orange2'), pt.cex=1.5,
         bty='n')
}

dev.off()

#########Gen Compare

pdf(file="Data/Paper1/gaic_total.pdf", width=15, height=7.5)


par(mfrow=c(1,2), mar=c(5,5,1,1))
GD1 <- readRDS("Data/Paper1/tn93StsubB_GD.rds")
GD2 <- readRDS("Data/Paper1/tn93NAsubB_GD.rds")

col1 <- "dodgerblue"
col2 <- "orange2"

r1 <- sapply(GD1, function(x) {x$gaic})
r2 <- sapply(GD2, function(x) {x$gaic})
x <- seq(0,0.04,0.0008)

par(mfrow=c(1,2), mar=c(5,5,1,1))
plot(x, r1, xlim=c(0, 0.04), xlab='TN93 distance cutoff',
     ylab='Generalized AIC', type='n', ylim=c(-80,12), cex.lab=1.5)

# create a background
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='ivory2', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(h=0, lty=2, lwd=3, col='grey50')
box()


lines(x, r1, col=col1, lwd=2)
points(x, r1, pch=22, col='white', bg=col1, cex=1.5)

text(x=x[which.min(r1)]+0.0035, y=min(r1), label="0.0160", cex=1.1)

lines(x, r2, col=col2, lwd=2)
points(x, r2, pch=22, col='white', bg=col2, cex=1.5)

text(x=x[which.min(r2)]-0.001, y=min(r2)-2, label=format(x[which.min(r2)],digits=3), cex=1.1)

legend("right", legend=c("Northern Alberta", "Seattle"), pch=c(22,22), 
       col='white', pt.bg=c(col2,col1), pt.cex=1.5, cex=1.5, bty='n')


######################################################

GD1 <- readRDS("~/Data/compareFig/compareMet_GD.rds")
GD2 <- readRDS("~/Data/compareFig/compareNM_GD.rds")

col1 <- "indianred4"
col2 <- "indianred1"

r1 <- sapply(GD1, function(x) {x$gaic})
r2 <- sapply(GD2, function(x) {x$gaic})
x <- seq(0,0.04,0.0008)

plot(x, r1, xlim=c(0, 0.04), xlab='TN93 distance cutoff',
     ylab='Generalized AIC', type='n', ylim=c(-80,12), cex.lab=1.5)

# create a background
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='ivory2', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(h=0, lty=2, lwd=3, col='grey50')
box()


lines(x, r1, col=col1, lwd=2)
points(x, r1, pch=22, col='white', bg=col1, cex=1.5)

text(x=x[which.min(r1)]+0.0035, y=min(r1), label=format(x[which.min(r1)],digits=4), cex=1.1)

lines(x, r2, col=col2, lwd=2)
points(x, r2, pch=22, col='white', bg=col2, cex=1.5)

text(x=x[which.min(r2)]-0.001, y=min(r2)-2, label="0.0160", cex=1.1)

legend("right", legend=c("Tennessee (Collection)", "Tennessee (Diagnostic)"), pch=c(22,22), 
       col='white', pt.bg=c(col2,col1), pt.cex=1.5, cex=1.5, bty='n')

dev.off()
