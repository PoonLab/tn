#Import Libraries
library(ggplot2,verbose = FALSE)
library(easynls,verbose = FALSE)

## args <- "~/WWData/WWGD"

GDs <- lapply(list.files(args), function(x) {readRDS(file=paste0(args, "/", x))})

gaics <- lapply(GDs, function(run){sapply(run, function(x) {x$gaic})})

dev.off()
minsLoc <- sapply(gaics, function(x){0.0008*(which(x==min(x))[[1]]-1)}) 
mins <- sapply(gaics, function(x){min(x)}) 
maxs <- sapply(gaics, function(x){max(x)})
minmin <- min(mins)
maxmax <- max(maxs)

cutoffs <- seq(0,0.04,0.0008)

clrs <- colours(distinct=T)[101:(101+length(cutoffs))] 
countries <- c("Argentina", "Austrailia", "Brazil", "Camaroon", "Canada", "Germany", "Italy", "Kenya", "Mexico", "Poland", "South Africa", "Spain", "Thailand", "Uganda", "United States")
countriesShort <- c("Arg", "Aus", "Bzl", "Cam", "Can", "Grm", "Itl", "Ken", "Mex", "Pol", "S.A.", "Spn", "Thai", "Ugnd", "U.S.")

plot.new()
title(xlab= "Cutoff values used to construct models and measure growth", 
      ylab = "GAIC: A null model's AIC subtracted from a proposed model AIC")
plot.window(xlim = c(0,0.04), ylim = c(minmin, maxmax))
axis(2, at= seq(-210,50,10), labels = seq(-210, 50, 10), las=2, pos = 0)
axis(3, minsLoc, las=3, lwd=0, lines=-2)
axis(1, at=cutoffs, labels=cutoffs)

i<-0
for(gaic in gaics) {
  i<-i+1
  lines(x=cutoffs,y=unname(gaic), col=clrs[i], lwd=1.2)
  abline(v=minsLoc[i], col=clrs[i], lwd=1, lty=2)
}

legend("bottomright", legend=countries, fill=clrs, title = "Sources of Data")


inp <- read.csv("~/WWData/NationalEpiData.csv",header=T)
df <- data.frame(row.names = countries, MGAICE = minsLoc, MinGAIC = mins, 
           Incidence=sapply(1:nrow(inp), function(x){mean((inp[x,1:3])[!is.na(inp[x,1:3])])}), 
           Prevalence=sapply(1:nrow(inp), function(x){mean((inp[x,4:6])[!is.na(inp[x,4:6])])}), 
           Diagnosed=inp$Diag, Treatment=inp$Treat, Suppressed=inp$Sup)



par(mfrow=c(2,2))

plot(df$MGAICE,df$Prevalence, xlab = "MGAICE: The point where the minimum GAIC is obtained", ylab= "Mean Prevalence (age 15-45), from 2005, 2010 and 2017 Data")
fit1<-lm(Prevalence~MGAICE, data =df)
lines(df$MGAICE, fitted(fit1))
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit1)$r.squared, digits=4)))
text(df$MGAICE, df$Prevalence, labels=countriesShort, cex= 0.7, pos=4)

plot(df$MGAICE,df$Incidence, xlab = "MGAICE: The point where the minimum GAIC is obtained", ylab= "Mean Incidence (per 1000 individuals), from 2005, 2010 and 2017 Data")
fit2<-lm(Incidence~MGAICE, data =df)
lines(df$MGAICE, fitted(fit2))
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit2)$r.squared, digits=4)))
text(df$MGAICE, df$Incidence, labels=countriesShort, cex= 0.7, pos=4)

plot(df$MGAICE,log(df$Prevalence), xlab = "MGAICE: The point where the minimum GAIC is obtained", ylab= "Mean Prevalence (age 15-45), from 2005, 2010 and 2017 Data", main = "LOG-SCALE")
fit3<-lm(log(Prevalence)~MGAICE, data =df)
lines(df$MGAICE, fitted(fit3))
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit3)$r.squared, digits=4)))
text(df$MGAICE, log(df$Prevalence), labels=countriesShort, cex= 0.7, pos=4)


plot(df$MGAICE,log(df$Incidence), xlab = "MGAICE: The point where the minimum GAIC is obtained", ylab= "Mean Incidence (per 1000 individuals), from 2005, 2010 and 2017 Data", main = "LOG-SCALE")
fit4<-lm(log(Incidence)~MGAICE, data =df)
lines(df$MGAICE, fitted(fit4))
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit4)$r.squared, digits=4)))
text(df$MGAICE, log(df$Incidence), labels=countriesShort, cex= 0.7, pos=4)


