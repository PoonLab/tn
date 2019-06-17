setwd('~/work/maup/data/pub1')

# sampling year distribution
naD <- read.csv('naD.txt', stringsAsFactors = FALSE)

# extract date info
naD.labels <- unique(c(naD$ID1, naD$ID2))

get.year <- function(x) {
   sapply(x, function(s) {
      as.integer(strsplit(s, "_")[[1]][2])
   })
}

naD.years <- get.year(naD.labels)

stD <- read.csv('stD.txt', stringsAsFactors = FALSE)
stD.labels <- unique(c(stD$ID1, stD$ID2))
stD.years <- get.year(stD.labels)

# prepare a data frame for barplot()
tab <- table(stD.years)
tab <- c(tab, '2013'=NA)
temp <- data.frame(year=names(tab), count.st=as.vector(tab))
tab <- table(naD.years)
temp$count.na <- as.vector(tab)[match(temp$year, names(tab))]


# prepare plot region
pdf(file='years.pdf', width=6, height=8)

par(mar=c(5,5,1,1), mfrow=c(2,1))
barplot(t(as.matrix(temp[,2:3])), beside=T, col='white', border=NA, 
        names.arg=temp$year, 
        xlab='Sample collection date (year)', 
        ylab='Number of cases', cex.lab=1.2, las=1)

# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=c(50, 100, 150), col='white', lwd=3, lend=2)
abline(h=seq(25, 200, 50), col='white', lend=3)

# draw axes (redraw y-axis)
axis(side=1, at=seq(2, 44, 6), labels=NA)
rug(x=seq(5, 41, 6), ticksize=-0.01, side=1)
axis(side=2)

# draw the bars in the foreground
barplot(t(as.matrix(temp[,2:3])), beside=T, add=T, axes=F,
        col=c('dodgerblue', 'orange2'), border='white')

legend(x=1, y=145, legend=c('Seattle', 'N.Alberta'), 
       fill=c('dodgerblue', 'orange2'), 
       y.intersp=1.2, bty='n')


# plot TN93 distributions


# collect histogram data
h1 <- hist(stD$Distance, breaks=50, plot=F)
h2 <- hist(naD$Distance, plot=F, breaks=h1$breaks)

n.st <- length(stD.labels)
n.na <- length(naD.labels)

barplot(h1$counts / choose(n.st,2), col=rgb(0.12,.56,1,.7), 
        border=rgb(0,0,0,0), space=0, xaxt='n',
        ylab='Frequency', xlab='TN93 distance', cex.lab=1.2)
axis(side=1, at=seq(0, length(h1$counts), 5),
     labels=seq(0, 0.05, length.out=11))
barplot(h2$counts / choose(n.na,2), add=T, col=rgb(.93,.604,.0,.7), space=0, 
        border=rgb(0,0,0,0))
#hist(stD$Distance, main=NA, border='white', col='dodgerblue')
legend(x=0.005, y=0.03, legend=c('Seattle', 'N.Alberta'),
       fill=c('dodgerblue', 'orange2'), bty='n')

dev.off()



#barplot(h1$counts / choose(n.st,2), col=rgb(0.12,.56,1,.7), 
#        border=rgb(0,0,0,0), space=0, xaxt='n', xlim=c(1,30),
#        ylab='Frequency', xlab='TN93 distance', cex.lab=1.2)
#axis(side=1, at=seq(0, length(h1$counts), 5),
#     labels=seq(0, 0.05, length.out=11))
#barplot(h2$counts / choose(n.na,2), add=T, col=rgb(.93,.604,.0,.7), space=0, 
#        border=rgb(0,0,0,0))



###########################