l1 <- readRDS("~/Data/Paper1/tn93StsubB_LD.rds")
l2 <- readRDS("~/Data/Paper1/tn93TnsubB_LD.rds")
l3 <- readRDS("~/Data/Paper1/tn93TnsubB_nomet_LD.rds")
l4 <- readRDS("~/Data/Paper1/tn93TnsubB_met_LD.rds")

cutoffs <- seq(0,0.04,0.0008)
#Plot Generation


plotLoss <- function(runs, colR, oFile, first=1) {
  #Obtain a list of vectors of GAICs for each filtered run
  gaics <- lapply(runs, function(run){sapply(run, function(x) {x$gaic})})
  
  #The step distance between cutoff points
  step <- max(cutoffs) / (length(cutoffs)-1)
  
  #The cutoff values which aquire the minimum GAIC. Also called the Minimum GAIC Estimator (MGAICE).
  minsLoc <- sapply(gaics, function(x){step*(which(x==min(x))[[1]]-1)}) 
  mins <- sapply(gaics, function(x){min(x)}) 
  maxs <- sapply(gaics, function(x){max(x)})
  
  #For defining range on a plot
  minmin <- -100
  maxmax <- 20
  maxT <- max((runs[[5]])[[1]]$v$Time)
  minT <- min((runs[[1]])[[1]]$v$Time)

  pdf(file = oFile, width=6, height=12)
  par(mfrow=c(5,1), mar = c(1,first,1,0) , oma=c(5,4,1,2))  
  
  #Make multiple plots for each run of GAICs with minimum labelled
  for (i in 1:length(gaics)) {
    
    #Initialize plot and background
    GAIC <- gaics[[i]]
    plot(cutoffs, GAIC, ylim = c(minmin+(0.2*minmin),maxmax), xlab="", ylab = "GAIC", cex.lab=1.8, cex.axis=1.5)
    bg <- par('usr')
    rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='linen', border=NA)
    abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
    abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
    abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
    abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
    
    #Plot GAIC
    lines(cutoffs, GAIC, lwd=2, col=colR)
    legend("bottomright", legend = paste0("Years ", minT, "-", maxT-(length(gaics)-i)), cex = 1.5,bg ="white")
    points(x=c(minsLoc[i]), y=c(mins[i]), cex=1)
    
    #Represents the location of the past run's MGAICE, Loss Ratio = minGAIC / Past minGAIC
    if (i>1){abline(v=minsLoc[i-1], lty=2)}
    
    #Draws an arrow to represent the follow through of the previous minimum
    if (i<5){arrows(minsLoc[i], mins[i], minsLoc[i], minmin+(0.2*minmin), length=0.05)}
    
    #To create the Cutoff label
    if (i==5){
      par(xpd=NA)
      title(xlab="Cutoffs", cex.lab=1.8)
    }
  }
  
  dev.off()
}

plotLoss(l1, "dodgerblue", "~/Data/Paper1/Loss_St.pdf", first=5)
plotLoss(l2, "orangered","~/Data/Paper1/Loss_Tn.pdf", first=5)
plotLoss(l3, "indianred1","~/Data/Paper1/Loss_Tn_Col.pdf" )
plotLoss(l4, "indianred4", "~/Data/Paper1/Loss_Tn_Diag.pdf")
