## args <- "~/Seattle/Robust/RLout"

args <- commandArgs(trailingOnly = T)

runs <- lapply(list.files(args), function(x) {readRDS(file=paste0(args, "/", x))})
gaics <- lapply(rev(runs), function(run){sapply(run, function(x) {x$gaic})})
minsLoc <- sapply(gaics, function(x){0.0008*(which(x==min(x))[[1]]-1)}) 

par(mfrow=c(2,3))

Cutoffs <- seq(0,0.04, 0.0008)

for (i in 1:length(gaics)) {
  GAIC <- gaics[[i]]
  plot(Cutoffs, GAIC, main = paste0("2000-", 2006+i))
  lines(Cutoffs, GAIC)
  abline(v=minsLoc[i], lty=2, lwt=1.5)
  if (i>1){
    abline(v=minsLoc[i-1], lty=2, col="orangered")
    lossRat <- GAIC[as.character(minsLoc[i-1])]/GAIC[as.character(minsLoc[i])]
    legend("bottomright", legend = paste0("Loss Ratio: ", round(lossRat,2)))
  }
}
