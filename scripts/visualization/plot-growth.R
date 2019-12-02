
GDst <- readRDS("~/Data/Paper1/tn93StsubB_GD.rds")
GDna <- readRDS("~/Data/Paper1/tn93NAsubB_GD.rds")
GDtn <- readRDS("~/Data/Paper1/tn93TnsubB_GD.rds")

GDtn_met <- readRDS("~/Data/Paper1/tn93TnsubB_met_GD.rds")
GDtn_NM <- readRDS("~/Data/Paper1/tn93TnsubB_nomet_GD.rds")


st <- sapply(GDst, function(x) {
  c(sum(x$g), sum(x$g>0))
})
st <- data.frame(cutoff=as.numeric(colnames(st)), total.R=st[1,], n.clust=st[2,])

na <- sapply(GDna, function(x) {
  c(sum(x$g), sum(x$g>0))
})
na <- data.frame(cutoff=as.numeric(colnames(na)), total.R=na[1,], n.clust=na[2,])

tn <- sapply(GDtn, function(x) {
  c(sum(x$g), sum(x$g>0))
})
tn <- data.frame(cutoff=as.numeric(colnames(tn)), total.R=tn[1,], n.clust=tn[2,])

met <- sapply(GDtn_met, function(x) {
  c(sum(x$g), sum(x$g>0))
})
met <- data.frame(cutoff=as.numeric(colnames(met)), total.R=met[1,], n.clust=met[2,])

NM <- sapply(GDtn_NM, function(x) {
  c(sum(x$g), sum(x$g>0))
})
NM <- data.frame(cutoff=as.numeric(colnames(NM)), total.R=NM[1,], n.clust=NM[2,])


pdf(file="~/Data/Paper1/growth.pdf", width=15, height=7.5)

par(mfrow=c(1,2), mar=c(5,5,4,2))
plot(st$cutoff, st$total.R, type='n', ylim=c(0, max(st$total.R)), cex.axis=1.5,
     xlab='TN93 distance cutoff', ylab='Cluster growth (R) | Number of active clusters', cex.lab=1.8)

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
lines(tn$cutoff, tn$total.R, type='s', col='orangered', lwd=2)
points(st$cutoff, st$n.clust, pch=23, col='white', bg='dodgerblue', cex=1.2)
points(na$cutoff, na$n.clust, pch=22, col='white', bg='orange2', cex=1.2)
points(tn$cutoff, tn$n.clust, pch=24, col='white', bg='orangered', cex=1.2)


legend("topleft", legend=c('Seattle', 'N.Alberta', 'Tennessee'), 
       pch=c(23, 22, 24), col=c('dodgerblue', 'orange2', 'orangered'), 
       pt.bg=c('dodgerblue', 'orange2','orangered'), 
       bty='n', lty=1, cex=1.5)

plot(met$cutoff, met$total.R, type='n', ylim=c(0, max(st$total.R)), cex.axis=1.5,
     xlab='TN93 distance cutoff', ylab='Cluster growth (R) | Number of active clusters', cex.lab=1.8)

# draw background
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='ivory2', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
box()

lines(met$cutoff, met$total.R, type='s', col='indianred1', lwd=2)
lines(NM$cutoff, NM$total.R, type='s', col='indianred4', lwd=2)
points(met$cutoff, met$n.clust, pch=24, col='white', bg='indianred1', cex=1.2)
points(NM$cutoff, NM$n.clust, pch=24, col='white', bg='indianred4', cex=1.2)

legend("topleft", legend=c('Tennessee (Collection Date)','Tennessee (Diagnostic Date)' ), 
       pch=c(24, 24), col=c('indianred1', 'indianred4'), 
       pt.bg=c('indianred1', 'indianred4'), 
       bty='n', lty=1, cex = 1.5)

dev.off()
