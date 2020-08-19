sGs <- readRDS("Data/ThesisData/sTs.rds")
nGs <- readRDS("Data/ThesisData/nTs.rds")
tGs <- readRDS("Data/ThesisData/tTs.rds")
tGs_Diag <- readRDS("Data/ThesisData/tTs_Diag.rds")

col1 <- "dodgerblue"
col2 <- "orange2"
col3 <- "indianred1"
col4 <- "indianred4"

#Number of Clusters
getNum <- function(iGs) {
  sapply(iGs, function(x){
    (nrow(x[[1]]$growth))
  })
}

#Total cluster growth
getGrowth <- function(iGs) {
  sapply(iGs, function(x){
    sum(x[[1]]$growth$New)
  })
}

#Proportion of growing clusters
getProp <- function(iGs) {
  sapply(iGs, function(x){
    nrow(x[[1]]$growth[which(x[[1]]$growth$New>0),])
  })
}

#Proportion of training set
getTrain <- function(iGs) {
  sapply(iGs, function(x){
    sum(x[[1]]$ageD$Positive) / max(x[[1]]$ageD$Total)
  })
}


cutoffs <- seq(0,0.15,0.003)
pdf(file="~/clusterPlots-trees.pdf", width=12, height=12)
par(mfrow=c(2,2),mar=c(5,6,4,1)+.1)

#PlotA
plot(x=0, y=0, type='n', xlab="Patristic distance threshold", ylab="Number of clusters", xlim=c(0,0.15), ylim=c(1,2800),
     cex.lab=1.8, cex.axis = 1.5)


title(outer=F, adj=0, "A", cex.main=2.5)

sD <- getNum(sGs)
nD <- getNum(nGs)
tD <- getNum(tGs)
tD_Diag <- getNum(tGs_Diag)

# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(500, 2500, 1000), col='white', lwd=3, lend=2)
abline(h=seq(0, 3000, 1000), col='white', lend=3)
abline(v=seq(0.05, 0.15, 0.05), col='white', lwd=3, lend=2)
abline(v=seq(0, 0.10, 0.05), col='white', lend=3)

points(cutoffs,sD, col=col1, cex=0.5)
lines(cutoffs,sD, lwd=1.5, col=col1)
points(cutoffs,nD, col=col2, cex=0.5)
lines(cutoffs,nD, lwd=1.5, col=col2)
points(cutoffs,tD, col=col3, cex=0.5)
lines(cutoffs,tD, lwd=1.5, col=col3)
points(cutoffs,tD_Diag, col=col4, cex=0.5)
lines(cutoffs,tD_Diag, lwd=1.5, col=col4)

legend('topright', legend=c('Seattle', 'Alberta','Tennessee (Collection Year)', 'Tennessee (Diagnostic Year)'),
       fill=c(col1, col2,col3, col4), bty='n', cex=1.8)


#PlotB

plot(x='n', y='n', xlab="Patristic distance threshold", ylab="Proportion of validation set involved", xlim=c(0,0.15), ylim=c(0,1),
     cex.lab=1.8, cex.axis = 1.5)

title(outer=F, adj=0, "B", cex.main=2.5)

sD <- getGrowth(sGs)/110
nD <- getGrowth(nGs)/110
tD <- getGrowth(tGs)/153
tD_Diag <- getGrowth(tGs_Diag)/129

# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(0, 1, 0.5), col='white', lwd=3, lend=2)
abline(h=seq(0.25, 0.75, 0.50), col='white', lend=3)
abline(v=seq(0.05, 0.15, 0.05), col='white', lwd=3, lend=2)
abline(v=seq(0, 0.10, 0.05), col='white', lend=3)


points(cutoffs,sD, col=col1, cex=0.5)
lines(cutoffs,sD, lwd=1.5, col=col1)
points(cutoffs,nD, col=col2, cex=0.5)
lines(cutoffs,nD, lwd=1.5, col=col2)
points(cutoffs,tD, col=col3, cex=0.5)
lines(cutoffs,tD, lwd=1.5, col=col3)
points(cutoffs,tD_Diag, col=col4, cex=0.5)
lines(cutoffs,tD_Diag, lwd=1.5, col=col4)


#PlotC
plot(x='n', y='n', xlab="Patristic distance threshold", ylab="Number of growing clusters", xlim=c(0,0.15), ylim=c(0,100),
     cex.lab=1.8, cex.axis = 1.5)

title(outer=F, adj=0, "C", cex.main=2.5)

sD <- getProp(sGs)
nD <- getProp(nGs)
tD <- getProp(tGs)
tD_Diag <- getProp(tGs_Diag)

# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(0, 80, 40), col='white', lwd=3, lend=2)
abline(h=seq(20, 100, 40), col='white', lend=3)
abline(v=seq(0.05, 0.15, 0.05), col='white', lwd=3, lend=2)
abline(v=seq(0, 0.10, 0.05), col='white', lend=3)


points(cutoffs,sD, col=col1, cex=0.5)
lines(cutoffs,sD, lwd=1.5, col=col1)
points(cutoffs,nD, col=col2, cex=0.5)
lines(cutoffs,nD, lwd=1.5, col=col2)
points(cutoffs,tD, col=col3, cex=0.5)
lines(cutoffs,tD, lwd=1.5, col=col3)
points(cutoffs,tD_Diag, col=col4, cex=0.5)
lines(cutoffs,tD_Diag, lwd=1.5, col=col4)


#PlotD

plot(x='n', y='n', xlab="Patristic distance threshold", ylab="Proportion of training outcomes remaining", xlim=c(0,0.15), ylim=c(0,1),
     cex.lab=1.8, cex.axis = 1.5, cex=1.5)

title(outer=F, adj=0, "D", cex.main=2.5)

sD <- getTrain(sGs)
nD <- getTrain(nGs)
tD <- getTrain(tGs)
tD_Diag <- getTrain(tGs_Diag)

# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(0, 1, 0.5), col='white', lwd=3, lend=2)
abline(h=seq(0.25, 0.75, 0.50), col='white', lend=3)
abline(v=seq(0.05, 0.15, 0.05), col='white', lwd=3, lend=2)
abline(v=seq(0, 0.10, 0.05), col='white', lend=3)


points(cutoffs,sD, col=col1, cex=0.5)
lines(cutoffs,sD, lwd=1.5, col=col1)
points(cutoffs,nD, col=col2, cex=0.5)
lines(cutoffs,nD, lwd=1.5, col=col2)
points(cutoffs,tD, col=col3, cex=0.5)
lines(cutoffs,tD, lwd=1.5, col=col3)
points(cutoffs,tD_Diag, col=col4, cex=0.5)
lines(cutoffs,tD_Diag, lwd=1.5, col=col4)



dev.off()
