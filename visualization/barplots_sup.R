tG_Diag <- readRDS("Data/ThesisData/tG_Diag.rds")
tG <- readRDS("Data/ThesisData/tG.rds")
sG <- readRDS("Data/ThesisData/sG.rds")
nG <- readRDS("Data/ThesisData/nG.rds")
tG <- readRDS("Data/ThesisData/tG.rds")

pdf(file="~/barplots_sup.pdf", width=15, height=15)
par(mfrow=c(2,1), mar=c(5,5,4,2))



tab1 <- table(sG$v$Time)
tab2 <- table(nG$v$Time)
tab3 <- table(tG$v$Time)

col1 <- "dodgerblue"
col2 <- "orange2"
col3 <- "orangered"

nmtot <- unique(c(names(tab1), names(tab2), names(tab3) ))
m <- matrix(dimnames=list(nmtot, c("Seattle", "Alberta", "Tennessee")),
            nrow=length(nmtot), ncol=3)

m[names(tab1),"Seattle"] <- as.numeric(tab1)
m[names(tab2),"North Alberta"] <- as.numeric(tab2)
m[names(tab3),"Tennessee"] <- as.numeric(tab3)

barplot(t(m), beside=T, col=c(col1,col2, col3),
        xlab='Sample collection date (year)',
        ylab='Number of sequences', cex.lab=1.8, las=1, cex.axis = 1.5, cex.names =1.5, las=0)

# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(50, 250, 50), col='white', lwd=3, lend=2)
abline(h=seq(25, 250, 50), col='white', lend=3)

legend('topleft', legend=c('Seattle', 'Alberta', 'Tennessee'), 
       fill=c('dodgerblue', 'orange2', 'orangered'), bty='n',cex = 1.8)

barplot(t(m), beside=T, col=c(col1,col2, col3),
        xlab='Sample collection date (year)', add=T, 
        ylab='Number of cases', cex.lab=1.8, las=1, cex.axis = 1.5, cex.names =1.5, las=0)


col1 <- 'indianred4'

tab2 <- table(tG_Diag$v$Time)[17:31]

barplot(tab2, beside=T, col=c(col1),
        xlab='HIV diagnosis date (year)', las=0,  
        ylab='Number of sequences', cex.lab=1.8, cex.axis = 1.5, cex.names = 1.5)


# create a background
x <- par('usr')
rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
abline(h=seq(50, 250, 50), col='white', lwd=3, lend=2)
abline(h=seq(25, 250, 50), col='white', lend=3)

legend('topleft', legend=c("Tennessee (Diagnostic Date Set)"),
       fill=c(col1, col2), bty='n', cex = 1.8)

barplot(tab2, beside=T, col=c(col1),
        xlab='HIV diagnosis date (year)', las=0, add=T, 
        ylab='Number of sequences', cex.lab=1.8, cex.axis = 1.5, cex.names = 1.5)


if(F) {
        axis(1, at=seq(2,56,3), labels=NA)
        
        
        h4 <- hist(tG_Diag$e$Distance[tG_Diag$e$Distance<0.05], breaks=50, plot=F)
        h5 <- hist(tG$e$Distance[tG$e$Distance<0.05], breaks=h4$breaks, plot=F)
        
        n4 <- length(tG_Diag$v$ID)
        n5 <- length(tG$v$ID)
        
        barplot(h4$counts / choose(n4,2), col=alpha(col1, 1), 
                border=rgb(0,0,0,0), space=0, xaxt='n',
                ylab='Frequency', xlab='TN93 distance', cex.lab=1.8, cex.axis = 1.5, cex.names = 1.5)
        
        #x <- par('usr')
        #rect(xl=x[1], yb=x[3], xr=x[2], yt=x[4], col='ivory2', border=NA)
        #abline(h=seq(0.005, 0.025, 0.005), col='white', lwd=3, lend=2)
        #abline(h=seq(0.0025, 0.025, 0.005), col='white', lend=3)
        
        
        axis(side=1, at=seq(0, length(h4$counts), 10),
             labels=seq(0, 0.05, 0.01), cex.axis=1.5)
        
        barplot(h4$counts / choose(n4,2), add=T, col=alpha(col1,1), space=0, 
                border=rgb(0,0,0,0),  cex.axis = 1.5)
        
        barplot(h5$counts / choose(n5,2), add=T, col=alpha(col2,0.75), space=0, 
                border=rgb(0,0,0,0),  cex.axis = 1.5)
        
        legend('topleft', legend=c("Tennessee (Collection Date Set)", "Tennessee (Diagnostic Date Set)"),
               fill=c(col1, col2), bty='n', cex = 1.8)
        
}

dev.off()
