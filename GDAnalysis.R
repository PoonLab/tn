#Analyze Growth data from tn93GD.R

### USAGE: Rscript GDAnalysis.R tn93GDOutput.RData ###

##Stat Presentation functions
#__________________________________________________________________________________________#

#Obtain the GAIC, a measure of fit between predicted and actual growth
gaic <- function(res=res)  {
  stat <- sapply(1:ncol(res), function(x) {
    #Extract full and fit data
    fit <- res[[1,x]]
    full <- res[[2,x]]
    
    #Place growth and forecast data in dfs for fit and full growth
    #df1 <- data.frame(Growth = fit$growth, Pred = fit$forecast)
    #df2 <- data.frame(Growth = full$growth, Pred = full$forecast)
    
    #Place growth and forecast data in dfs for fit and full growth
    df1 <- data.frame(Growth = fit$growth, Pred = fit$forecast)
    df2 <- data.frame(Growth = fit$growth, Pred = fit$csize * (sum(fit$growth)/sum(fit$csize)))
    
    #Model growth as a function of forecast for fit and full growth models
    mod1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
    mod2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
    
    #Calculate GAIC
    mod1$aic-mod2$aic
  })
  
  #Present and return data
  plot(colnames(res), stat, ylab = "GAIC", xlab = "Cutoff", type= "l", cex.lab=1.5)
  return(stat)
}

#__________________________________________________________________________________________#

#Expecting the output from a run of tn93GD.Rdata
args = commandArgs(trailingOnly = T)

#Loads the output from tn93GD.RData
load(args)

## Skeleton Code For plotting
stat <- sapply(1:ncol(res), function(x) {
  fit <- res[[1,x]]
  print(fit)
})
plot(colnames(res), stat, ylab = "" , xlab= "")