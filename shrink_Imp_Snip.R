##TO-DO: Specify source-file Location
source("~/git/tn/comp_An.R")

g <- readRDS("~/Data/Paper1/tn93TnsubB_nomet_G.rds")
oFile <- "~/Data/Paper1/tn93TnsubB_nomet"

#remV <- rev(g$v$ID)[1:98]
#g$v <- subset(g$v, !ID%in%remV)
#g$e <- subset(g$e, !ID1%in%remV & !ID2%in%remV)
#g$f <- bpeFreq(g)


#Create a set of subgraphs with sequentially smaller new-years filtered subgraphs
gs  <- lapply(1:(unname(tail(table(g$v$Time),1))-1), function(x){
  iG <- g
  remV <- rev(iG$v$ID)[1:x]
  iG$v <- subset(iG$v, !ID%in%remV)
  iG$e <- subset(iG$e, !ID1%in%remV & !ID2%in%remV)
  iG$f <- likData(iG)
  
  print(nrow(subset(iG$v, Time==max(iG$v$Time))))
  
  return(iG)
}) 

#Run the GAIC algorithm for all subgraphs
runs <- lapply(1:length(gs) , function(i) {
  
  print(paste0(round(i/length(gs)*100),"%"))
  iG<-gs[[i]]
  
  gaicRun(iG)
})

#Save all growth data in accessable files
saveRDS(runs, file = paste0(oFile, "_SD.rds"))

#Obtain a list of vectors of GAICs for each filtered run
cutoffs <- as.numeric(names(runs[[1]])) 

minsLoc <- sapply(gaics, function(x){step*(which(x==min(x))[[1]]-1)}) 


plot(79:1, mins, xlab = "Number of New Cases", ylab="Absolute Minimum GAIC Obtained")

rMet <- readRDS("~/Data/Tennessee/analysis/tn93TnsubB_met_SD.rds")
rNM <- readRDS("~/Data/Tennessee/analysis/tn93TnsubB_nomet_SD.rds")

g1 <- lapply(rMet, function(run){sapply(run, function(x) {x$gaic})})
g2 <- lapply(rNM, function(run){sapply(run, function(x) {x$gaic})})

s1 <- lapply(rMet, function(run){nrow(subset(run[[1]]$v, Time==max(Time)))})
s2 <- lapply(rNM, function(run){nrow(subset(run[[1]]$v, Time==max(Time)))})

s1 <- lapply(rMet, function(run){table(run[[1]]$v$Time)})
s2 <- lapply(rNM, function(run){table(run[[1]]$v$Time)})


m1 <- sapply(g1, function(x){min(x)})
m2 <- sapply(g2, function(x){min(x)})

pdf(file="~/Data/Paper1/shrink.pdf", width=7.5, height=7.5)

plot(79:1, rev(abs(m1)), xlab='Number of New Cases in subset (|U|)',
     ylab='Magnitude of Minimum GAIC', type='n', ylim=c(0,20), cex.lab=1.5)

# draw background
bg <- par('usr')
rect(xl=bg[1], yb=bg[3], xr=bg[2], yt=bg[4], col='ivory2', border=NA)
abline(h=axTicks(side=2), col='white', lwd=3, lend=2)
abline(h=axTicks(side=2)+diff(axTicks(side=2))[1]/2, col='white', lend=2)
abline(v=axTicks(side=1), col='white', lwd=3, lend=2)
abline(v=axTicks(side=1)+diff(axTicks(side=1))[1]/2, col='white', lend=2)
box()

lines(rev(abs(m1)), col="indianred4", lwd=3)
lines(rev(abs(m2)), col="indianred1", lwd=3)

legend('topleft', legend=c('Tennessee: Diagnostic Date', 'Tennessee: Collection Date' ),
       fill=c('indianred4', 'indianred1'), bty='n')


dev.off()