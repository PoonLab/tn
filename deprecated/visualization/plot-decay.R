tn <- T

stD <- readRDS("Data/Paper1/tn93StsubB_GD.rds")
naD <- readRDS("Data/Paper1/tn93NAsubB_GD.rds")
if(tn) {
  tnD <- readRDS("Data/Paper1/tn93TnsubB_GD.rds")
  tnD_met <- readRDS("Data/Paper1/tn93TnsubB_met_GD.rds")
  tnD_NM <- readRDS("Data/Paper1/tn93TnsubB_nomet_GD.rds")
}


# log transformed plot  
if(tn){
  pdf(file='Data/Paper1/decay.pdf', width=15, height=7.5)
} else {
  pdf(file='Data/Paper1/decay.pdf', width=5, height=5)
}

par(mfrow=c(1,2), mar=c(5,5,4,2))
limits <- par('usr')

# Tennessee
if(tn) {
  ageD <- lapply(tnD, function(x) x$f)
  ageDi <- ageD[["0.04"]]
  ageDi <- subset(ageDi, vTotal>63 & tDiff<14 & oeDens<280)
}else{
  ageD <- lapply(stD, function(x) x$f)
  ageDi <- ageD[["0.04"]]
}

mod <- glm(cbind(Positive, vTotal) ~ tDiff, data=ageDi, family='binomial')

set.seed(1)  # for reproducible jitter
plot(jitter(ageDi$tDiff), 
     ageDi$Positive/(ageDi$vTotal), 
     ylim=c(0.005,0.5), 
     log='y',
     type='n',
     xlab='Collection Time lag (In Years)', 
     ylab=expression('Bipartite edge density'), 
     cex.lab=1.8, yaxt='n', cex.axis=1.5, las=0)

# draw background
bg <- par('usr')
rect(xl=bg[1], yb=10^(bg[3]), xr=bg[2], yt=10^(bg[4]), col='ivory2', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
#abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(h=0, lty=2, lwd=3, col='grey50')
box()

points(jitter(ageDi$tDiff), mod$y, 
       pch=ifelse(tn,24,21), bg=ifelse(tn,'orangered', 'dodgerblue'), col='white', cex=ifelse(tn,1,1.5))

axis(2, at=axTicks(side=2), 
     labels=axTicks(side=2), # * 10^4,
     las=0, cex.axis = 1.5)

if(tn) {
  legend("topright", legend=c('Tennessee', 'N.Alberta', 'Seattle'), pch=c(24,22,23), pt.cex=1.5,
         col='white', pt.bg=c('orangered', 'orange2', 'dodgerblue'), bty='n', cex = 1.8)
} else {
  legend("topright", legend=c('Seattle', 'N.Alberta'), pch=c(21,22), pt.cex=1.5,
         col='white', pt.bg=c('dodgerblue', 'orange2'), bty='n', cex = 1.8)
}


# indicate location of zeroes
y <- ageDi$Positive/(ageDi$vTotal)
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

mod <- glm(cbind(Positive, vTotal) ~ tDiff, data=ageDi, family='binomial')

#plot(Positive/Total ~ tDiff, data=ageDi, log='y')

set.seed(4)
points(jitter(ageDi$tDiff), mod$y, pch=22, 
       bg='orange2', col='white', cex=ifelse(tn,1,1.5))
lines(smooth.spline(ageDi$tDiff, mod$fitted.values),
      lwd=2, col='orange2')

# Tennessee
if(tn) {
  ageD <- lapply(stD, function(x) x$f)
  ageDi <- ageD[["0.04"]]
  
  mod <- glm(cbind(Positive, vTotal) ~ tDiff, data=ageDi, family='binomial')
  
  set.seed(4)
  points(jitter(ageDi$tDiff), mod$y, pch=23, 
         bg='dodgerblue', col='white', cex=1)
  lines(smooth.spline(ageDi$tDiff, mod$fitted.values),
        lwd=2, col='dodgerblue')
}

#############################Second Plot for JUST tnMET


# Tennessee_met
if(tn) {
  ageD <- lapply(tnD_met, function(x) x$f)
  ageDi <- ageD[["0.04"]]
  ageDi <- subset(ageDi,  tDiff<25&vTotal>63 & oeDens<280)
}else{
  ageD <- lapply(stD, function(x) x$f)
  ageDi <- ageD[["0.04"]]
}

mod <- glm(cbind(Positive, vTotal) ~ tDiff, data=ageDi, family='binomial')

set.seed(1)  # for reproducible jitter
plot(jitter(ageDi$tDiff), 
     ageDi$Positive/(ageDi$vTotal), 
     ylim=c(0.005,0.5), 
     log='y',
     type='n',
     xlab='Time lag (years)', 
     ylab=expression('Bipartite edge density'), # %*%10^-4), 
     cex.lab=1.8, yaxt='n', cex.axis=1.5)

# draw background
bg <- par('usr')
rect(xl=bg[1], yb=10^(bg[3]), xr=bg[2], yt=10^(bg[4]), col='ivory2', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
#abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
abline(h=0, lty=2, lwd=3, col='grey50')
box()

points(jitter(ageDi$tDiff), mod$y, 
       pch=ifelse(tn,24,23), bg=ifelse(tn,"indianred4", 'dodgerblue'), col='white', cex=ifelse(tn,1,1.5))

axis(2, at=axTicks(side=2), 
     labels=axTicks(side=2), # * 10^4,
     las=0, cex.axis=1.5)

# indicate location of zeroes
y <- ageDi$Positive/(ageDi$vTotal)
ymin <- min(y[y>0], na.rm=T)
#points(ageDi$tDiff[ageDi$Positive==0], 
       #rep(ymin, sum(ageDi$Positive==0)), pch=4, lwd=2, col=ifelse(tn,'orangered', 'dodgerblue'))
lines(smooth.spline(x=ageDi$tDiff, y=mod$fitted.values), 
      lwd=2, col=ifelse(tn,'indianred4', 'dodgerblue'))


ageD <- lapply(tnD_NM, function(x) x$f)
ageDi <- ageD[["0.04"]]
ageDi <- subset(ageDi, vTotal>63 & tDiff<14 & oeDens<280)

mod <- glm(cbind(Positive, vTotal) ~ tDiff, data=ageDi, family='binomial')

#plot(Positive/Total ~ tDiff, data=ageDi, log='y')

set.seed(4)
points(jitter(ageDi$tDiff), mod$y, pch=24, 
       bg='indianred1', col='white', cex=ifelse(tn,1,1.5))
lines(smooth.spline(ageDi$tDiff, mod$fitted.values),
      lwd=2, col='indianred1')

legend("topright", legend=c('Tennessee (Diagnostic Dates)', 'Tennessee (Collection Dates)'), pch=c(24,24), pt.cex=1.5,
       col='white', pt.bg=c('indianred4', 'indianred1'), bty='n', cex = 1.8)

dev.off()
