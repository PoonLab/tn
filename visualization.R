require(scales)

res <- resA
robF <- "Data/paperData/naFullTree/naFull_ROB.rds"
resComp <- NA
hmark <- numeric(0)
cCol <- "black"
cCol <- "grey"

plotGAIC <- function(res, resComp=NA, robF="", vmark=numeric(0), mCol="black", cCol="grey", randMod=F,
                     legLabs=numeric(0)) {
  
  #Searching for Robustness run info
  if(file.exists(robF)) {
    resRob <- readRDS(robF)
  } else {
    resRob <- res[numeric(0)]
  }
  
  #Searching for Comparison run info
  if(is.na(comp)) {
    resComp <- res[numeric(0)]
  } 

  x <- c(res$MaxDistance, resComp$MaxDistance, resRob$MaxDistance)
  y <- c(res$GAIC, resComp$GAIC, resRob$GAIC)
  plot(x, y, type='n', ylab="AIC Loss", xlab = "Bootstrap Threshold", xaxs = "i")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "antiquewhite")
  ticStep <- (par()$yaxp[2]-par()$yaxp[1])/(par()$yaxp[3])
  ticLoc <- seq(par()$yaxp[1],par()$yaxp[2], ticStep)
  abline(h=ticLoc, col="white", lwd=2.5)
  abline(h=(ticLoc+ticStep/2), col="white", lwd=1)
  abline(h=0, lty=2)
  lines(res$MaxDistance, res$GAIC, col=mCol, lwd=3)
  
  robLine <- sapply(res$MaxDistance, function(d) {
    dGAIC <- resRob[MaxDistance==d, (GAIC)]
    return(c(mean(dGAIC), sd(dGAIC)))
  })
  lines(res$MaxDistance, robLine[1,], col=cCol)
  polygon(c(res$MaxDistance, rev(res$MaxDistance)),
          c(robLine[1,]+robLine[2,],rev(robLine[1,]-robLine[2,])),
          col=alpha("black", 0.08), border = alpha("black", 0.65), lty=3)
}