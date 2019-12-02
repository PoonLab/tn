tnD_met <- readRDS("Data/Paper1/tn93TnsubB_met_G.rds")
tnD_NM <- readRDS("Data/Paper1/tn93TnsubB_nomet_G.rds")

pdf(file="Data/Paper1/barplots_sup.pdf", width=15, height=15)

tab1 <- table(tnD_NM$v$Time)
tab2 <- table(tnD_met$v$Time)

col1 <- "indianred1"
col2 <- "indianred4"

par(mfrow=c(2,1), mar=c(5,5,4,2))
nmtot <- unique(c(names(tab2), names(tab1)))[17:35]
m <- matrix(dimnames=list(nmtot, c("Tennessee (Collection Date Set)", "Tennessee (Diagnostic Date Set)")),
            nrow=length(nmtot), ncol=2)

m[names(tab1),"Tennessee (Collection Date Set)"] <- as.numeric(tab1)
m[names(tab2)[17:31],"Tennessee (Diagnostic Date Set)"] <- as.numeric(tab2)[17:31]

barplot(t(m), beside=T, col=c(col1,col2),
        xlab='Time Points (year)', las=0,
        ylab='Number of cases', cex.lab=1.8, cex.axis = 1.5, cex.names = 1.5)

# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(50, 250, 50), col='white', lwd=3, lend=2)
abline(h=seq(25, 250, 50), col='white', lend=3)

legend('topleft', legend=c("Tennessee (Collection Date Set)", "Tennessee (Diagnostic Date Set)"),
       fill=c(col1, col2), bty='n', cex = 1.8)

#axis(side=2)
barplot(t(m), beside=T, col=c(col1,col2), add=T,
        xlab='Time Points (year)', 
        ylab='Number of cases', cex.lab=1.8, las=0, cex.axis = 1.5, cex.names = 1.5)

axis(1, at=seq(2,56,3), labels=NA)


h4 <- hist(tnD_met$e$Distance[tnD_met$e$Distance<0.05], breaks=h1$breaks, plot=F)
h5 <- hist(tnD_NM$e$Distance[tnD_NM$e$Distance<0.05], breaks=h1$breaks, plot=F)

n4 <- length(tnD_met$v$ID)
n5 <- length(tnD_NM$v$ID)

barplot(h4$counts / choose(n4,2), col=alpha(col1, 1), 
        border=rgb(0,0,0,0), space=0, xaxt='n',
        ylab='Frequency', xlab='TN93 distance', cex.lab=1.8, cex.axis = 1.5, cex.names = 1.5)

x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(0.005, 0.025, 0.005), col='white', lwd=3, lend=2)
abline(h=seq(0.0025, 0.025, 0.005), col='white', lend=3)


axis(side=1, at=seq(0, length(h4$counts), 5),
     labels=seq(0, 0.05, length.out=11), cex.axis=1.5)

barplot(h4$counts / choose(n4,2), add=T, col=alpha(col1,1), space=0, 
        border=rgb(0,0,0,0),  cex.axis = 1.5)

barplot(h5$counts / choose(n5,2), add=T, col=alpha(col2,0.75), space=0, 
        border=rgb(0,0,0,0),  cex.axis = 1.5)

legend('topleft', legend=c("Tennessee (Collection Date Set)", "Tennessee (Diagnostic Date Set)"),
       fill=c(col1, col2), bty='n', cex = 1.8)

dev.off()