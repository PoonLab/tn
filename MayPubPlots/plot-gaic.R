# contains the cluster information for graphs at different distance cutoffs
GDst <- readRDS("stDGD.rds")
GDna <- readRDS("naDGD.rds")


st <- sapply(GDst, function(x) {x$gaic})
na <- sapply(GDna, function(x) {x$gaic})
x <- seq(0.005,0.05,0.001)

pdf(file='~/papers/maup/images/gaic.pdf', width=5, height=5)

par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(x, st, xlim=c(0.005, 0.04), xlab='TN93 distance cutoff',
     ylab='Generalized AIC', type='n', cex.lab=1.2)

# create a background
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='ivory2', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(h=0, lty=2, lwd=3, col='grey50')
box()

lines(x, st, col='dodgerblue', lwd=2)
points(x, st, pch=21, col='white', bg='dodgerblue', cex=1.1)

text(x=x[which.min(st)]+0.003, y=min(st), label=x[which.min(st)], cex=0.8)

lines(x, na, col='orange2', lwd=2)
points(x, na, pch=22, col='white', bg='orange2', cex=1.1)

text(x=x[which.min(na)]-0.001, y=min(na)-2, label=x[which.min(na)], cex=0.8)

legend(x=0.025, y=-30, legend=c('Seattle', 'N.Alberta'), pch=c(21,22), 
       col='white', pt.bg=c('dodgerblue', 'orange2'), pt.cex=1.5,
       bty='n')

dev.off()