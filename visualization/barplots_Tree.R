require(scales)


sT <- readRDS("Data/ThesisData/sT.rds")
nT <- readRDS("Data/ThesisData/nT.rds")
tT <- readRDS("Data/ThesisData/tT.rds")
tT_Diag <- readRDS("Data/ThesisData/tT_Diag.rds")

pdf(file="~/barplots_Tree.pdf", width=15, height=15)
par(mar=c(5,6,4,1)+.1, mfrow=c(2,1))
col1 <- "dodgerblue"
col2 <- "orange2"
col3 <- "indianred1"
col4 <- "indianred4"


temp <- sT$e$dist[1:nrow(sT$v), 1:nrow(sT$v)]
sT$el <- unlist(lapply(1:(length(temp[1,])-1), function(i){temp[i,(i+1):length(temp[1,])]}))

temp <- nT$e$dist[1:nrow(nT$v), 1:nrow(nT$v)]
nT$el <- unlist(lapply(1:(length(temp[1,])-1), function(i){temp[i,(i+1):length(temp[1,])]}))

temp <- tT$e$dist[1:nrow(tT$v), 1:nrow(tT$v)]
tT$el <- unlist(lapply(1:(length(temp[1,])-1), function(i){temp[i,(i+1):length(temp[1,])]}))

temp <- tT_Diag$e$dist[1:nrow(tT_Diag$v), 1:nrow(tT_Diag$v)]
tT_Diag$el <- unlist(lapply(1:(length(temp[1,])-1), function(i){temp[i,(i+1):length(temp[1,])]}))

h1 <- hist(sT$el[sT$el<0.15], breaks=50, plot=F)
h2 <- hist(nT$el[nT$el<0.15], breaks=h1$breaks, plot=F)

n1 <- length(sT$v$ID)
n2 <- length(nT$v$ID)

barplot((h1$counts / length(sT$el)), col=alpha("white", 1), 
        border=rgb(0,0,0,0), space=0, xaxt='n',
        ylab='Frequency', xlab='Patristic distance', cex.lab=1.8, cex.axis = 1.5)

axis(side=1, at=seq(0, length(h1$counts), 5),
     labels=seq(0, 0.15, 0.025), cex.axis=1.5)

barplot(h1$counts/length(sT$el), add=T, col=alpha(col1,1), space=0, 
        border=rgb(0,0,0,0), cex.axis = 1.5, cex.names = 1.5)

barplot(h2$counts/length(nT$el), add=T, col=alpha(col2,0.85), space=0, 
       border=rgb(0,0,0,0), cex.axis = 1.5, cex.names = 1.5)


legend('topleft', legend=c('Seattle', 'Alberta'),
       fill=c(col1, col2), bty='n', cex=1.8)


#Col V Diag
h3 <- hist(tT$el[sT$el<0.15], breaks=50, plot=F)
h4 <- hist(tT_Diag$el[nT$el<0.15], breaks=h3$breaks, plot=F)

n1 <- length(tT$v$ID)
n2 <- length(tT_Diag$v$ID)

barplot((h1$counts / length(sT$el)), col=alpha("white", 1), 
        border=rgb(0,0,0,0), space=0, xaxt='n',
        ylab='Frequency', xlab='Patristic distance', cex.lab=1.8, cex.axis = 1.5)

axis(side=1, at=seq(0, length(h1$counts), 5),
     labels=seq(0, 0.15, 0.025), cex.axis=1.5)

barplot(h4$counts/length(tT_Diag$el), add=T, col=alpha(col4,1), space=0, 
        border=rgb(0,0,0,0), cex.axis = 1.5, cex.names = 1.5)

barplot(h3$counts/length(tT$el), add=T, col=alpha(col3,0.85), space=0, 
        border=rgb(0,0,0,0), cex.axis = 1.5, cex.names = 1.5)


legend('topleft', legend=c('Tennessee (Collection Date Set) ', 'Tennessee (Diagnostic Date Set)'),
       fill=c(col3, col4), bty='n', cex=1.8)

dev.off()
