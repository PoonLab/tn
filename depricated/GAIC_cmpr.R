#Quick comparison of two sets of GAICs in a single plot 

# Plug in 2 sets of GAICs here
gaicsF <- NA
gaicsT <- NA

#Plot and compare two sets of GAICs
plot(cutoffs, gaicsF, type = "n", ylim=c(min(c(gaicsF,gaicsT)),max(c(gaicsF,gaicsT))), xlab="Cutoffs", ylab = "GAIC")
lines(cutoffs, gaicsF, lwd=1.6, col="orange")
lines(cutoffs, gaicsT, lwd=1.6, col="orangered")
points(cutoffs, gaicsF)
points(cutoffs, gaicsT)
abline(h=0, lty=2)
abline(v=cutoffs[which(gaicsF==min(gaicsF))[[1]]], lty=3)
abline(v=cutoffs[which(gaicsT==min(gaicsT))[[1]]], lty=3)
text(cutoffs[which(gaicsF==min(gaicsF))[[1]]+1], min(gaicsF), labels= round(min(gaicsF)))
text(cutoffs[which(gaicsT==min(gaicsT))[[1]]+1], min(gaicsT), labels= round(min(gaicsT)))
text(cutoffs[which(gaicsF==min(gaicsF))[[1]]], max(c(gaicsT,gaicsF)), labels=cutoffs[which(gaicsF==min(gaicsF))[[1]]])
text(cutoffs[which(gaicsT==min(gaicsT))[[1]]], max(c(gaicsT,gaicsF))-1.5, labels=cutoffs[which(gaicsT==min(gaicsT))[[1]]])
