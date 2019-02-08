#Analyze Growth data from tn93GD.R

### USAGE: Rscript GDAnalysis.R tn93GDOutput.RData ###

####- TO-DO: Add PPmap functionality to this and GD function -####

#__________________________________________________________________________________________#

#Obtain the VPC, a measure of variance due to the aggregate level over total variance
VPC <- function(res=res){
  stat <- sapply(colnames(res), function(x) {
    #Extract full and fit data
    fit <- res[[1,x]]
    full <- res[[2,x]]
    
    #Calculate VPC
    var(fit$growth/fit$csize) / var(full$growth/full$csize)
  })
  
  #Present and return data
  plot(stat, ylab = "VPC", xlab = colnames(res))
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
  plot(stat, ylab = "GAIC", xlab = colnames(res))
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
    mod2$deviance-mod1$deviance 
  })
  
  #Present and return data
  plot(stat, ylab = "Deviance", xlab = colnames(res))
  return(stat)
}

#__________________________________________________________________________________________#

#Expecting the output from a run of tn93GD.Rdata
args = commandArgs(trailingOnly = T)

#Loads the output from tn93GD.RData
load(args)

#Test Functions
VPC <- vpc()
GAIC <- gaic()
deviance <- dev()