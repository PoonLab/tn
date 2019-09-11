tn <- T

stD <- readRDS("Data/Paper1/tn93StsubB_GD.rds")
naD <- readRDS("Data/Paper1/tn93NAsubB_GD.rds")
if(tn) {
  tnD <- readRDS("Data/Tennessee/analysis_PRO/tn93TnsubB_met_GD.rds")
}


# log transformed plot  
if(tn){
  pdf(file='Data/Paper1/decay_tn.pdf', width=7.5, height=7.5)
} else {
  pdf(file='Data/Paper1/decay.pdf', width=5, height=5)
}

par(mar=c(5,5,1,1))
limits <- par('usr')

# Tennessee
if(tn) {
  ageD <- lapply(tnD, function(x) x$f)
  ageDi <- ageD[["0.04"]]
  ageDi <- subset(ageDi, tDiff<15)
}else{
  ageD <- lapply(stD, function(x) x$f)
  ageDi <- ageD[["0.04"]]
}


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
       pch=ifelse(tn,24,21), bg=ifelse(tn,'orangered', 'dodgerblue'), col='white', cex=ifelse(tn,0.75,1.5))

axis(2, at=axTicks(side=2), 
     labels=axTicks(side=2), # * 10^4,
     las=2)

if(tn) {
  legend(x=12, y=0.55, legend=c('Tennessee', 'N.Alberta', 'Seattle'), pch=c(24,21,22), pt.cex=1.5,
         col='white', pt.bg=c('orangered', 'orange2', 'dodgerblue'), bty='n')
} else {
  legend(x=8, y=0.5, legend=c('Seattle', 'N.Alberta'), pch=c(21,22), pt.cex=1.5,
         col='white', pt.bg=c('dodgerblue', 'orange2'), bty='n')
}


# indicate location of zeroes
y <- ageDi$Positive/ageDi$Total
ymin <- min(y[y>0], na.rm=T)
points(ageDi$tDiff[ageDi$Positive==0], 
       rep(ymin, sum(ageDi$Positive==0)), pch=4, lwd=2, col=ifelse(tn,'orangered', 'dodgerblue'))
lines(smooth.spline(x=ageDi$tDiff, y=mod$fitted.values), 
      lwd=2, col=ifelse(tn,'orangered', 'dodgerblue'))


# show confidence intervals
#pr <- predict(mod, se=T, type='response')
#ll <- pr$fit - 1.96*pr$se.fit
#lines(smooth.spline(ageDi$tDiff, ll), lty=2)


# North Alberta

ageD <- lapply(naD, function(x) x$f)
ageDi <- ageD[["0.04"]]

mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')

#plot(Positive/Total ~ tDiff, data=ageDi, log='y')

set.seed(4)
points(jitter(ageDi$tDiff), mod$y, pch=22, 
       bg='orange2', col='white', cex=ifelse(tn,0.75,1.5))
lines(smooth.spline(ageDi$tDiff, mod$fitted.values),
      lwd=2, col='orange2')

# Tennessee
if(tn) {
  ageD <- lapply(stD, function(x) x$f)
  ageDi <- ageD[["0.04"]]
  
  mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')
  
  set.seed(4)
  points(jitter(ageDi$tDiff), mod$y, pch=21, 
         bg='dodgerblue', col='white', cex=0.75)
  lines(smooth.spline(ageDi$tDiff, mod$fitted.values),
        lwd=2, col='dodgerblue')
}

dev.off()
