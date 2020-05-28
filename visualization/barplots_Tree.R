require(scales)


stD <- readRDS("Data/Seattle/stst_T.rds")
naD <- readRDS("Data/NAlberta/naT.rds")
tnD <- readRDS("Data/Tennessee/tnst_Trim_Diag_T.rds")

pdf(file="~/barplots_Tree.pdf", width=15, height=8)
par(mar=c(5,6,4,1)+.1)
col1 <- "dodgerblue"
col2 <- "orange2"
col3 <- "orangered"

distMat <- naD$e$dist[1:nrow(naD$v),1:nrow(naD$v)]
d2 <- as.numeric(unlist(lapply(1:nrow(iT$v), function(x){distMat[x,x:nrow(iT$v)]})))

h1 <- hist(stD$e$el$Distance[stD$e$el$Distance<0.15], breaks=50, plot=F)
h2 <- hist(d2[which(d2<0.15)], breaks=h1$breaks, plot=F)
h3 <- hist(tnD$e$el$Distance[tnD$e$el$Distance<0.15], breaks=h1$breaks, plot=F)

n1 <- length(stD$v$ID)
n2 <- length(naD$v$ID)
n3 <- length(tnD$v$ID)



barplot((h3$counts / nrow(tnD$e$el)), col=alpha("white", 1), 
        border=rgb(0,0,0,0), space=0, xaxt='n',
        ylab='Frequency', xlab='Patristic distance', cex.lab=1.8, cex.axis = 1.5)

#x <- par('usr')
#rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
#abline(h=seq(0.02, 0.12, 0.02), col='white', lwd=3, lend=2)
#abline(h=seq(0.01, 0.11, 0.02), col='white', lend=3)

axis(side=1, at=seq(0, length(h1$counts), 5),
     labels=seq(0, 0.150, 0.025), cex.axis=1.5)


#h2a <- c((h2$counts / choose(n2,2))[1:30] , rep(NA,20))
#h2b <- c(rep(NA,30), (h2$counts / choose(n2,2))[31:50])

barplot(h1$counts/nrow(stD$e$el), add=T, col=alpha(col1,1), space=0, 
        border=rgb(0,0,0,0), cex.axis = 1.5, cex.names = 1.5)

barplot(h2$counts/length(d2), add=T, col=alpha(col2,0.85), space=0, 
       border=rgb(0,0,0,0), cex.axis = 1.5, cex.names = 1.5)

#barplot(h2a, add=T, col=alpha(col2,0.5), space=0, 
 #       border=rgb(0,0,0,0), cex.axis = 0, cex.names = 1.5)
#barplot(h2b, add=T, col=alpha(col2,1), space=0, 
 #       border=rgb(0,0,0,0), cex.axis = 0, cex.names = 1.5)

barplot(h3$counts/nrow(tnD$e$el), add=T, col=alpha(col3,0.5), space=0, 
        border=rgb(0,0,0,0), cex.axis = 0, cex.names = 1.5)


legend('topleft', legend=c('Seattle', 'N.Alberta', 'Tennessee'),
       fill=c('dodgerblue', 'orange2', 'orangered'), bty='n', cex=1.8)
dev.off()
