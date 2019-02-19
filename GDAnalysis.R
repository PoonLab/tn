#Analyze Growth data from tn93GD.R

### USAGE: Rscript GDAnalysis.R tn93GDOutput.RData ###

##Stat Presentation functions
#__________________________________________________________________________________________#

#Obtain the VPC, a measure of variance due to the aggregate level over total variance
vpc <- function(res=res){
  stat <- sapply(colnames(res), function(x) {
    #Extract full and fit data
    fit <- res[[1,x]]
    full <- res[[2,x]]
    
    #Calculate VPC
    var(fit$growth/fit$csize) / var(full$growth/full$csize)
  })
  
  #Present and return data
  plot(colnames(res), stat, ylab = "VPC", xlab = colnames(res))
  return(stat)
}

#Obtain the GAIC, a measure of fit between predicted and actual growth
gaic <- function(res=res)  {
  stat <- sapply(colnames(res), function(x) {
    #Extract full and fit data
    fit <- res[[1,x]]
    full <- res[[2,x]]
    
    #Place growth and forecast data in dfs for fit and full growth
    df1 <- data.frame(Growth = fit$growth, Pred = fit$forecast)
    df2 <- data.frame(Growth = full$growth, Pred = full$forecast)
    
    #Model growth as a function of forecast for fit and full growth models
    mod1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
    mod2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
    
    #Calculate GAIC
    mod1$aic-mod2$aic
  })
  
  #Present and return data
  plot(colnames(res), stat, ylab = "GAIC", xlab = "Cutoff")
  return(stat)
}

#Obtain the Deviance, a classical fit measurement
dev <- function(res=res) {
  
  stat <- sapply(colnames(res), function(x) {
    #Extract full and fit data
    fit <- res[[1,x]]
    full <- res[[2,x]]
    
    #Place growth and forecast data in dfs for fit and full growth
    df1 <- data.frame(Growth = fit$growth, Pred = fit$forecast)
    df2 <- data.frame(Growth = full$growth, Pred = full$forecast)
    
    #Model growth as a function of forecast for fit and full growth models
    mod1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
    mod2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
    
    #Calculate Deviance
    1-mod1$deviance/mod2$deviance 
  })
  
  #Present and return data
  plot(colnames(res), stat, ylab = "Deviance", xlab = "Cutoff")
  return(stat)
}

#Obtain the Poisson Probability map of relative growth
ppmap <- function(res=res) {
  
  stat <- sapply(colnames(res), function(x) {
    #Extract fit data
    fit <- res[[1,x]]
    
    #Obtain an expectation based off of overall cluster growth
    exp <- fit$inc*(fit$csize)/(sum(fit$csize))
    
    #zscore can be calculated based off of this diff
    diff <- fit$growth - exp
    zscore <- diff / sqrt(exp)
    return(var(zscore))
  })
  
  #Present and return data
  plot(colnames(res), stat, ylab = "Poisson Probability", xlab = "Cutoff")
  return(stat)
}

#Obtain the proportion of growing clusters
gr <- function(res=res) {
  
  stat <- sapply(1:(ncol(res)-1), function(x) {
    #Extract fit data
    fit <- res[[1,x]]
    
    #Obtain the ratio of growing clusters over total clusters
    length(fit$growth[fit$growth>0])/fit$no
  })
  
  #Present and return data
  plot(stat, ylab = "Growth Ratio", xlab = "Cutoff")
  return(stat)
}

check <- function(res) {
  stat <- sapply(1:ncol(res), function(x) {
    fit <- res[[1,x]]
    g <- fit$growth/fit$csize
    print(g[g>0])
  })

}
  


#__________________________________________________________________________________________#

#Expecting the output from a run of tn93GD.Rdata
args = commandArgs(trailingOnly = T)

#Loads the output from tn93GD.RData
load(args)