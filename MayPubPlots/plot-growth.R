
GDst <- readRDS("stDGD.rds")
GDna <- readRDS("naDGD.rds")

st <- sapply(GDst, function(x) {
  c(sum(x$growth), sum(x$growth>0))
})
st <- data.frame(cutoff=as.numeric(colnames(st)), total.R=st[1,], n.clust=st[2,])

na <- sapply(GDna, function(x) {
  c(sum(x$growth), sum(x$growth>0))
})
na <- data.frame(cutoff=as.numeric(colnames(na)), total.R=na[1,], n.clust=na[2,])


pdf(file='~/papers/maup/images/growth.pdf', width=5, height=5)

par(mar=c(5,5,1,1))
plot(st$cutoff, st$total.R, type='n', ylim=c(0, max(st$total.R)),
     xlab='TN93 distance cutoff', ylab='Cluster growth (R) / Number of active clusters', cex.lab=1.2)

# draw background
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='ivory2', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
box()

lines(st$cutoff, st$total.R, type='s', col='dodgerblue', lwd=2)
lines(na$cutoff, na$total.R, type='s', col='orange2', lwd=2)
points(st$cutoff, st$n.clust, pch=21, col='white', bg='dodgerblue', cex=1.2)
points(na$cutoff, na$n.clust, pch=22, col='white', bg='orange2', cex=1.2)

legend(x=0.005, y=105, legend=c('Seattle', 'N.Alberta'), 
       pch=c(21, 22), col=c('dodgerblue', 'orange2'), 
       pt.bg=c('dodgerblue', 'orange2'), 
       bty='n', lty=1)

dev.off()
