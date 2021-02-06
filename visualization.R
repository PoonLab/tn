require(scales)


#Plot and compare AIC loss results
plotGAIC <- function(res, resComp=NA, robF="", mCol="black", cCol="grey", randMod=F,
                     legLabs=NA, xlim=NA, figLab="", ylim=NA) {
  #@param res: The result of a GAIC run. This is the "Main" result
  #@param resComp: The result of a second GAIC run for comparison
  #@param robF: A file containing a robust run (ie. repeated runs from the multiGAIC function)
  #@param mCol: The colour for the main result
  #@param cCol: The colour for the comparison result.
  #             Robustness info will also use this
  #@param randMod: If comparing to a random model
  #@param legLabs: The text to go into the legend describing res, and resComp/resRob
  #@param xlim: Manually set xrange for magnification
  #@param ylim: Manually set yrange for magnification
  #@param figLab: A label for the top left hand of the image. For multiple plotting
  
  #Searching for Robustness run info
  if(file.exists(robF)) {
    resRob <- readRDS(robF)
  } else {
    resRob <- res[numeric(0)]
  }
  
  #Calculating robustness info
  robLine <- sapply(res$MaxDistance, function(d) {
    dGAIC <- resRob[MaxDistance==d, (GAIC)]
    return(c(mean(dGAIC), sd(dGAIC)))
  })
  
  #Searching for Comparison run info
  if(is.na(resComp)) {
    resComp <- res[numeric(0)]
  } 

  #Create plot 
  x <- c(res$MaxDistance, resComp$MaxDistance, resRob$MaxDistance)
  y <- c(res$GAIC, resComp$GAIC, resRob$GAIC)
  
  par(mar=c(5,5,5.2,2)+0.1)
  if(is.na(xlim)) {xlim = range(x)}
  if(is.na(ylim)) {ylim = range(y)*1.10}
  plot(x, y, type='n', ylab="AIC Difference from Null Model", xlab = "Maximum Distance Threshold", xaxs = "i",
       xlim=xlim, ylim= ylim, cex.axis=1.4, cex.lab=1.65)
  
  #Background
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "antiquewhite")
  ticStep <- (par()$yaxp[2]-par()$yaxp[1])/(par()$yaxp[3])
  ticLoc <- seq(par()$yaxp[1],par()$yaxp[2], ticStep)
  abline(h=ticLoc, col="white", lwd=2.5)
  abline(h=(ticLoc+ticStep/2), col="white", lwd=1)
  abline(h=0, lty=2)
  
  abline(v=res$MaxDistance[which.min(res$GAIC)], col="grey", lty=2, lwd=2)
  abline(v=res$MaxDistance[which.min(resComp$GAIC)], col="grey", lty=2, lwd=2)
  abline(v=res$MaxDistance[which.min(robLine[1,])], col="grey", lty=2, lwd=2)
  
  #Plot comp
  lines(resComp$MaxDistance, resComp$GAIC, col=cCol, lwd=3)
  
  #Plot Rob
  lines(res$MaxDistance, robLine[1,], col=cCol, lwd=2)
  polygon(c(res$MaxDistance, rev(res$MaxDistance)),
          c(robLine[1,]+robLine[2,],rev(robLine[1,]-robLine[2,])),
          col=alpha("black", 0.08), border = alpha("black", 0.65), lty=1, lwd=0.5)
  
  #Plot main
  lines(res$MaxDistance, res$GAIC, col=mCol, lwd=3)
  
  #Define steps for text labels on plot
  vStep <- (par("usr")[4]-par("usr")[3])/50
  hStep <- (par("usr")[2]-par("usr")[1])/50
  if(nrow(resRob)>0) {
    if(which.min(res$GAIC)>which.min(robLine[1,])) {hStep <- -hStep}
  }
  if(nrow(resComp)>0) {
    if(which.min(res$GAIC)>which.min(resComp$GAIC)) {hStep <- -hStep}
  }
  
  #Highlight min loc for Main
  lines(rep(res$MaxDistance[which.min(res$GAIC)], 2), c(par("usr")[3], par("usr")[3]+vStep*3), col=mCol, lwd=4)
  lines(rep(res$MaxDistance[which.min(res$GAIC)], 2), c(par("usr")[3], par("usr")[3]+vStep*2), col="white", lwd=0.5)
  points(res$MaxDistance[which.min(res$GAIC)], par("usr")[3])
  text(res$MaxDistance[which.min(res$GAIC)]-hStep*2, par("usr")[3]+vStep*5, res$MaxDistance[which.min(res$GAIC)], cex=1.25)
  
  #Highlight min loc for Comp (if different)
  if(nrow(resComp)>0){
    if(which.min(res$GAIC)!=which.min(resComp$GAIC)) {
      lines(rep(res$MaxDistance[which.min(resComp$GAIC)], 2), c(par("usr")[3], par("usr")[3]+vStep*2), col=cCol, lwd=4)
      lines(rep(res$MaxDistance[which.min(resComp$GAIC)], 2), c(par("usr")[3], par("usr")[3]+vStep*1), col="white", lwd=0.5)
      points(res$MaxDistance[which.min(resComp$GAIC)], par("usr")[3])
      text(res$MaxDistance[which.min(resComp$GAIC)]+hStep*2, par("usr")[3]+vStep*4, res$MaxDistance[which.min(resComp$GAIC)], cex=1.25)
    }
  }
  
  if(nrow(resRob)>0){
    #Highlight min loc for Rob (if different)
    if(which.min(res$GAIC)!=which.min(robLine[1,])) {
      lines(rep(res$MaxDistance[which.min(robLine[1,])], 2), c(par("usr")[3], par("usr")[3]+vStep*2), col=cCol, lwd=4)
      lines(rep(res$MaxDistance[which.min(robLine[1,])], 2), c(par("usr")[3], par("usr")[3]+vStep*1), col="white", lwd=0.5)
      points(res$MaxDistance[which.min(robLine[1,])], par("usr")[3])
      text(res$MaxDistance[which.min(robLine[1,])]+hStep*2, par("usr")[3]+vStep*4, res$MaxDistance[which.min(robLine[1,])], cex=1.25)
    }
  }

  
  #Legend and figure labelling
  par(xpd=NA)
  if(!is.na(legLabs)){
    legend("bottomright", legend = legLabs, fill = c(mCol, cCol), bty='n',cex=1.5) 
  }
  shiftx <- par('usr')[1] - strwidth(figLab, cex=2)*2
  shifty <- par('usr')[4] + strheight(figLab, cex=2)*3
  text(shiftx, shifty, figLab, cex=2.4)
}

#Plot and compare AIC loss results
plotMinLoc <- function(robFs, cols, legLabs=NA, figLab="", xlim=NA, bw="nrd0") {
  #@param robF: A set of files containing robust runs (ie. repeated runs from the multiGAIC function)
  #@param cols: A set of colours for the results
  #@param legLabs: The text to go into the legend describing res, and resComp/resRob
  #@param figLab: A label for the top left hand of the image. For multiple plotting
  #@param xlim: Manually set xrange for magnification
  
  #Obtain all optima locations
  minLocs <- lapply(robFs, function(f){
    resRob <- readRDS(f)
    sapply(range(resRob$RunID)[1]:range(resRob$RunID)[2], function(i){
      x <- resRob[RunID==i]
      x[which.min(GAIC), (MaxDistance)]
    })
  })
  
  #Get density plots for optima locations
  ds <- lapply(minLocs, function(x){
    density(x, bw=bw)
  })
  
  #Create a plot
  par(mar=c(5,5,5.2,2)+0.1)
  x <- unlist(lapply(ds, function(d){d$x}))
  y <- unlist(lapply(ds, function(d){d$y}))
  if(is.na(xlim)) {xlim = range(x)}
  plot(x,y,type='n', xlim=xlim, xaxs="i", yaxs="i", xlab="Optimal Threshold", ylab="Density",
       cex.axis=1.4, cex.lab=1.65)
  
  #Background
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "seashell1")
  ticStep <- (par()$xaxp[2]-par()$xaxp[1])/(par()$xaxp[3])
  ticLoc <- seq(par()$xaxp[1],par()$xaxp[2], ticStep)
  abline(v=ticLoc, col="white", lwd=3)
  abline(v=(ticLoc+ticStep/4), col="white", lwd=1.8)
  abline(v=(ticLoc+2*ticStep/4), col="white", lwd=1.8)
  abline(v=(ticLoc+3*ticStep/4), col="white", lwd=1.8)
  abline(h=0, lty=2)
  
  
  #Add density lines
  for(i in 1:length(ds)){
    l <- ds[[i]]
    cl <- cols[[i]]
    lines(l, col=cl, lwd=3)
    polygon(l$x, l$y,col = alpha(cols[[i]], 0.20))
  }
  
  #Legend and figure labelling
  par(xpd=NA)
  if(!is.na(legLabs)){
    legend("topright", legend = legLabs, fill = cols, bty='n',cex=1.5) 
  }
  shiftx <- par('usr')[1] - strwidth(figLab, cex=2)*2
  shifty <- par('usr')[4] + strheight(figLab, cex=2)*3
  text(shiftx, shifty, figLab, cex=2.4)
}

plotEdgeLength <- function(treeL, cols, legLabs=NA, figLab="", xlim=NA) {
  
  blSpan <- sapply(treeL, function(t){max(t$edge.length)})
  breaks <- seq(0, max(blSpan), max(blSpan)/10)
  
  #Get hist plots for optima locations
  hs <- lapply(treeL, function(t) {
    x <- t$edge.length
    hist(x, breaks=breaks)
  })
  
  #Create a plot
  par(mar=c(5,5,5.2,2)+0.1)
  hts <- sapply(1:(length(breaks)-1), function(i){
    sapply(hs, function(h){h$counts[i]/sum(h$counts)})
  })
  
  if(is.na(xlim)) {xlim = range(x)}
  
  barplot(hts, beside=T)
  axis(side = 1, at=seq(0.5, (length(treeL)+1)*length(breaks), length(treeL)+1))

  
  #Background
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "seashell1")
  ticStep <- (par()$xaxp[2]-par()$xaxp[1])/(par()$xaxp[3])
  ticLoc <- seq(par()$xaxp[1],par()$xaxp[2], ticStep)
  abline(v=ticLoc, col="white", lwd=3)
  abline(v=(ticLoc+ticStep/4), col="white", lwd=1.8)
  abline(v=(ticLoc+2*ticStep/4), col="white", lwd=1.8)
  abline(v=(ticLoc+3*ticStep/4), col="white", lwd=1.8)
  abline(h=0, lty=2)
  
  
  #Add density lines
  for(i in 1:length(ds)){
    l <- hs[[i]]
    cl <- cols[[i]]
    
  }
  
  #Legend and figure labelling
  par(xpd=NA)
  if(!is.na(legLabs)){
    legend("topright", legend = legLabs, fill = cols, bty='n',cex=1.5) 
  }
  shiftx <- par('usr')[1] - strwidth(figLab, cex=2)*2
  shifty <- par('usr')[4] + strheight(figLab, cex=2)*3
  text(shiftx, shifty, figLab, cex=2.4)
  
  
}