stD <- readRDS("Data/Paper1/tn93StsubB_G.rds")
naD <- readRDS("Data/Paper1/tn93NAsubB_G.rds")

tnD <- readRDS("Data/Paper1/tn93TnsubB_G.rds")



pdf(file="Data/Paper1/barplots.pdf", width=15, height=15)

tab1 <- table(stD$v$Time)
tab2 <- table(naD$v$Time)
tab3 <- table(tnD$v$Time)

col1 <- "dodgerblue"
col2 <- "orange2"
col3 <- "orangered"

par(mfrow=c(2,1))
nmtot <- unique(c(names(tab1), names(tab2), names(tab3) ))
m <- matrix(dimnames=list(nmtot, c("Seattle", "North Alberta", "Tennessee")),
            nrow=length(nmtot), ncol=3, xlab=NA)

m[names(tab1),"Seattle"] <- as.numeric(tab1)
m[names(tab2),"North Alberta"] <- as.numeric(tab2)
m[names(tab3),"Tennessee"] <- as.numeric(tab3)

barplot(t(m), beside=T, col=c(col1,col2,col3),
        xlab='Sample collection date (year)', 
        ylab='Number of cases', cex.lab=1.2, las=1)

# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(50, 250, 50), col='white', lwd=3, lend=2)
abline(h=seq(25, 250, 50), col='white', lend=3)

legend('topleft', legend=c('Seattle', 'N.Alberta', 'Tennessee'),
       fill=c('dodgerblue', 'orange2', 'orangered'), bty='n')

#axis(side=2)

barplot(t(m), beside=T, col=c(col1,col2,col3), add=T,
        xlab='Sample collection date (year)', 
        ylab='Number of cases', cex.lab=1.2, las=1)


axis(1, at=seq(2.5,62.5,4), labels=nmtot)

h1 <- hist(stD$e$Distance[stD$e$Distance<0.05], breaks=50, plot=F)
h2 <- hist(naD$e$Distance[naD$e$Distance<0.05], breaks=h1$breaks, plot=F)
h3 <- hist(tnD$e$Distance[tnD$e$Distance<0.05], breaks=h1$breaks, plot=F)

n1 <- length(stD$v$ID)
n2 <- length(naD$v$ID)
n3 <- length(tnD$v$ID)

barplot(h1$counts / choose(n1,2), col=alpha(col1, 1), 
        border=rgb(0,0,0,0), space=0, xaxt='n',
        ylab='Frequency', xlab='TN93 distance', cex.lab=1.2)

x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(0.005, 0.025, 0.005), col='white', lwd=3, lend=2)
abline(h=seq(0.0025, 0.025, 0.005), col='white', lend=3)


axis(side=1, at=seq(0, length(h1$counts), 5),
     labels=seq(0, 0.05, length.out=11))


h2a <- c((h2$counts / choose(n2,2))[1:30] , rep(NA,20))
h2b <- c(rep(NA,30), (h2$counts / choose(n2,2))[31:50])

barplot(h1$counts / choose(n1,2), add=T, col=alpha(col1,1), space=0, 
        border=rgb(0,0,0,0))

barplot(h2a, add=T, col=alpha(col2,0.5), space=0, 
        border=rgb(0,0,0,0))
barplot(h2b, add=T, col=alpha(col2,1), space=0, 
        border=rgb(0,0,0,0))

barplot(h3$counts / choose(n3,2), add=T, col=alpha(col3,0.75), space=0, 
        border=rgb(0,0,0,0))

legend('topleft', legend=c('Seattle', 'N.Alberta', 'Tennessee'),
       fill=c('dodgerblue', 'orange2', 'orangered'), bty='n')
dev.off()
