
stD <- readRDS("tn93StGD.rds")
naD <- readRDS("tn93NAGD.rds")

# log transformed plot  
pdf(file='decay.pdf', width=5, height=5)

par(mar=c(5,5,1,1))
limits <- par('usr')


ageD <- lapply(stD, function(x) x$ageD)
ageDi <- ageD[["0.04"]]
mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')

set.seed(1)  # for reproducible jitter
plot(jitter(ageDi$tDiff), 
     ageDi$Positive/ageDi$Total, 
     log='y',
     type='n',
     xlab='Time lag (years)', 
     ylab=expression('Bipartite edge density'), # %*%10^-4), 
     cex.lab=1.2, yaxt='n')

# draw background
bg <- par('usr')
rect(xl=bg[1], yb=10^(bg[3]), xr=bg[2], yt=10^(bg[4]), col='ivory2', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
#abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(h=0, lty=2, lwd=3, col='grey50')
box()

points(jitter(ageDi$tDiff), ageDi$Positive/ageDi$Total, 
       pch=21, bg='dodgerblue', col='white', cex=1.5)

axis(2, at=axTicks(side=2), 
     labels=axTicks(side=2), # * 10^4,
     las=2)

legend(x=8, y=0.01, legend=c('Seattle', 'N.Alberta'), pch=c(21,22), pt.cex=1.5,
       col='white', pt.bg=c('dodgerblue', 'orange2'), bty='n')

# indicate location of zeroes
y <- ageDi$Positive/ageDi$Total
ymin <- min(y[y>0], na.rm=T)
points(ageDi$tDiff[ageDi$Positive==0], 
       rep(ymin, sum(ageDi$Positive==0)), pch=4, lwd=2, col='dodgerblue')
lines(smooth.spline(x=ageDi$tDiff, y=mod$fitted.values), 
      lwd=2, col='dodgerblue')

# show confidence intervals
#pr <- predict(mod, se=T, type='response')
#ll <- pr$fit - 1.96*pr$se.fit
#lines(smooth.spline(ageDi$tDiff, ll), lty=2)


# North Alberta

ageD <- lapply(naD, function(x) x$ageD)
ageDi <- ageD[["0.04"]]

mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')

#plot(Positive/Total ~ tDiff, data=ageDi, log='y')

set.seed(4)
points(jitter(ageDi$tDiff), mod$y, pch=22, 
       bg='orange2', col='white', cex=1.5)
lines(smooth.spline(ageDi$tDiff, mod$fitted.values),
      lwd=2, col='orange2')

dev.off()
