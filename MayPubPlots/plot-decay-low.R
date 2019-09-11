tn <- T

stD <- readRDS("Data/Paper1/tn93StsubB_GD.rds")
naD <- readRDS("Data/Paper1/tn93NAsubB_GD.rds")
if(tn) {
  tnD <- readRDS("Data/Tennessee/analysis_PRO/tn93TnsubB_nomet_GD.rds")
}

par(mar=c(5,5,1,1))
limits <- par('usr')

ageD <- lapply(stD, function(x) x$f)
cutoffs <- names(ageD)
r2 <- {}
for (cutoff in cutoffs) {
   ageDi <- ageD[[cutoff]]
   mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')
   r2 <- c(r2, 1-mod$deviance/mod$null.deviance)
}


par(mar=c(5,5,1,1))
plot(cutoffs, r2, type='l', col='dodgerblue', ylim=c(0.2, 0.65),
     xlab='TN93 distance cutoff', ylab=expression("R"^2), cex.lab=1.2)
points(cutoffs, r2, pch=21, col='white', bg='dodgerblue', cex=1.2)


ageD <- lapply(stD, function(x) x$f)
cutoffs <- names(ageD)
r2 <- {}
for (cutoff in cutoffs) {
   ageDi <- ageD[[cutoff]]
   mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')
   r2 <- c(r2, 1-mod$deviance/mod$null.deviance)
}

lines(cutoffs, r2, col='orange2')
points(cutoffs, r2, pch=21, col='white', bg='orange2', cex=1.2)

#========================================================#

# log transformed plot  
pdf(file='decaylow.pdf', width=4.5, height=4.5)

par(mar=c(5,5,1,1))
ageD <- lapply(stD, function(x) x$f)
ageDi <- ageD[["0.004"]]
mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')

set.seed(1)  # for reproducible jitter
plot(jitter(ageDi$tDiff), 
     ageDi$Positive/ageDi$Total, 
#     log='y',
     type='n',
     xlab='Time lag (years)', 
     ylab=expression('Bipartite edge density '), 
     cex.lab=1.2, yaxt='n')
 
# draw background
#bg <- par('usr')
#rect(xl=bg[1], yb=10^(bg[3]), xr=bg[2], yt=10^(bg[4]), col='linen', border=NA)
#abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
##abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
#abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
#abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
#abline(h=0, lty=2, lwd=3, col='grey50')
#box()

points(jitter(ageDi$tDiff), ageDi$Positive/ageDi$Total, 
       pch=21, bg='dodgerblue', col='white', cex=1.5)

axis(2, at=axTicks(side=2), 
     labels=axTicks(side=2),
     las=2)

legend(x=7, y=0.08, legend=c('Seattle', 'N.Alberta'), pch=c(21,22), pt.cex=1.5,
       col='white', pt.bg=c('dodgerblue', 'orange2'), bty='n')

# indicate location of zeroes
#y <- ageDi$Positive/ageDi$Total
#ymin <- min(y[y>0], na.rm=T)
#points(jitter(ageDi$tDiff[ageDi$Positive==0]), 
#       rep(ymin, sum(ageDi$Positive==0)), pch=4, lwd=2, col='dodgerblue')
lines(smooth.spline(x=ageDi$tDiff, y=mod$fitted.values), 
      lwd=2, col='dodgerblue')

# show confidence intervals
#pr <- predict(mod, se=T, type='response')
#ll <- pr$fit - 1.96*pr$se.fit
#lines(smooth.spline(ageDi$tDiff, ll), lty=2)


# North Alberta

ageD <- lapply(naD, function(x) x$f)
ageDi <- ageD[["0.004"]]

mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')

#plot(Positive/Total ~ tDiff, data=ageDi, log='y')

set.seed(4)
points(jitter(ageDi$tDiff), mod$y, pch=22, 
       bg='orange2', col='white', cex=1.5)
lines(smooth.spline(ageDi$tDiff, mod$fitted.values),
      lwd=2, col='orange2')

dev.off()
