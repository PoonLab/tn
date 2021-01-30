require(scales)

sG <- readRDS("Data/ThesisData/sG.rds")
nG <- readRDS("Data/ThesisData/nG.rds")
tG <- readRDS("Data/ThesisData/tG.rds")
tG_Diag <- readRDS("Data/ThesisData/tG_Diag.rds")


pdf(file="~/barplots.pdf", width=15, height=15)
par(mfrow=c(2,1))

#tab1 <- table(sG$v$Time)
#tab2 <- table(nG$v$Time)
#tab3 <- table(tG$v$Time)

col1 <- "dodgerblue"
col2 <- "orange2"
col3 <- "orangered"
col3 <- "indianred1"
col4 <- "indianred4"

h1 <- hist(sG$e$Distance[sG$e$Distance<0.05], breaks=50, plot=F)
h2 <- hist(nG$e$Distance[nG$e$Distance<0.05], breaks=h1$breaks, plot=F)
h3 <- hist(tG$e$Distance[tG$e$Distance<0.05], breaks=h1$breaks, plot=F)

n1 <- length(sG$v$ID)
n2 <- length(nG$v$ID)
n3 <- length(tG$v$ID)
n4 <- length(tG_Diag$v$ID)

barplot(h3$counts / choose(n3,2), col=alpha(col1, 1), 
        border=rgb(0,0,0,0), space=0, xaxt='n',
        ylab='Frequency', xlab='TN93 distance', cex.lab=1.8, cex.axis = 1.5)


axis(side=1, at=seq(0, length(h1$counts), 5),
     labels=seq(0, 0.05, length.out=11), cex.axis=1.5)


h2a <- c((h2$counts / nrow(sG$v$e))[1:30] , rep(NA,20))
h2b <- c(rep(NA,30), (h2$counts / choose(n2,2))[31:50])

barplot(h1$counts / choose(n1,2), add=T, col=alpha(col1,1), space=0, 
        border=rgb(0,0,0,0), cex.axis = 1.5, cex.names = 1.5)

barplot(h2a, add=T, col=alpha(col2,0.5), space=0, 
        border=rgb(0,0,0,0), cex.names=1.5)
barplot(h2b, add=T, col=alpha(col2,1), space=0, 
        border=rgb(0,0,0,0), cex.axis = 1.5, cex.names = 1.5)

barplot(h3$counts / choose(n3,2), add=T, col=alpha(col3,0.75), space=0, 
        border=rgb(0,0,0,0), cex.axis =1.5, cex.names = 1.5)

legend('topleft', legend=c('Seattle', 'Alberta', 'Tennessee'),
       fill=c('dodgerblue', 'orange2', 'orangered'), bty='n', cex=1.8)




h3 <- hist(tG$e$Distance[tG$e$Distance<0.05], breaks=h1$breaks, plot=F)
h4 <- hist(tG_Diag$e$Distance[tG_Diag$e$Distance<0.05], breaks=h1$breaks, plot=F)

barplot(h3$counts / choose(n3,2), col=alpha(col3, 1), 
        border=rgb(0,0,0,0), space=0, xaxt='n',
        ylab='Frequency', xlab='TN93 distance', cex.lab=1.8, cex.axis = 1.5)

axis(side=1, at=seq(0, length(h1$counts), 5),
     labels=seq(0, 0.05, length.out=11), cex.axis=1.5)

barplot(h4$counts / choose(n4,2), add=T, col=alpha(col4,1), space=0, 
        border=rgb(0,0,0,0), cex.axis = 1.5, cex.names = 1.5)

barplot(h3$counts / choose(n3,2), add=T, col=alpha(col3,0.70), space=0, 
        border=rgb(0,0,0,0), cex.axis = 1.5, cex.names = 1.5)


legend('topleft', legend=c('Tennessee (Collection Date Set) ', 'Tennessee (Diagnostic Date Set)'),
       fill=c(col3, col4), bty='n', cex=1.8)


dev.off()
