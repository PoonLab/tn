## args <- "~/WWData/WWGD"

GDs <- lapply(list.files(args), function(x) {readRDS(file=paste0(args, "/", x))})

gaics <- lapply(GDs, function(run){sapply(run, function(x) {x$gaic})})
size <- sapply(GDs, function(run){sum(run[[1]]$csize)})
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
df <- data.frame(row.names = countries, MGAICE = minsLoc, MinGAIC = mins, SampleS = size, 
           Incidence=sapply(1:nrow(inp), function(x){mean((inp[x,1:3])[!is.na(inp[x,1:3])])}), 
           Prevalence=sapply(1:nrow(inp), function(x){inp[x,11]*0.01*mean((inp[x,4:6])[!is.na(inp[x,4:6])])}), 
           Diagnosed=inp$Diag, Treatment=inp$Treat, Suppressed=inp$Sup)



par(mfrow=c(2,2))

plot(df$MGAICE,df$Prevalence, xlab = "MGAICE: The point where the minimum GAIC is obtained", ylab= "Mean Prevalence, from 2005, 2010 and 2017 Data")
fit1<-lm(Prevalence~MGAICE, data =df)
lines(df$MGAICE, fitted(fit1))
legend("topright", bty="n", legend=paste0("Pearson Correlation: ", format(cor(df$MGAICE, df$Prevalence, method="pearson"), digits=4)))
text(df$MGAICE, df$Prevalence, labels=countriesShort, cex= 0.7, pos=4)

plot(df$MGAICE,df$SampleS/df$Prevalence, xlab = "MGAICE: The point where the minimum GAIC is obtained", ylab= "Sample Size, Relative to Prevalence")
fit1<-lm(SampleS/Prevalence~MGAICE, data =df)
lines(df$MGAICE, fitted(fit1))
legend("topright", bty="n", legend=paste0("Pearson Correlation: ", format(cor(df$MGAICE, df$SampleS/df$Prevalence, method="pearson"), digits=4)))
text(df$MGAICE, df$SampleS/df$Prevalence, labels=countriesShort, cex= 0.7, pos=4)


plot(df$MGAICE,log(df$Prevalence), xlab = "MGAICE: The point where the minimum GAIC is obtained", ylab= "Mean Prevalence , from 2005, 2010 and 2017 Data", main = "LOG-SCALE")
fit3<-lm(log(Prevalence)~MGAICE, data =df)
lines(df$MGAICE, fitted(fit3))
legend("topright", bty="n", legend=paste0("Pearson Correlation: ", format(cor(df$MGAICE, log(df$Prevalence), method="pearson"), digits=4)))
text(df$MGAICE, log(df$Prevalence), labels=countriesShort, cex= 0.7, pos=4)

plot(df$MGAICE,log(df$SampleS/df$Prevalence), xlab = "MGAICE: The point where the minimum GAIC is obtained", ylab= "Sample Size, Relative to Prevalence", main = "LOG-SCALE")
fit3<-lm(log(df$SampleS/Prevalence)~MGAICE, data =df)
lines(df$MGAICE, fitted(fit3))
legend("topright", bty="n", legend=paste0("Pearson Correlation: ", format(cor(df$MGAICE, log(df$SampleS/df$Prevalence), method="pearson"), digits=4)))
text(df$MGAICE, log(df$SampleS/df$Prevalence), labels=countriesShort, cex= 0.7, pos=4)




