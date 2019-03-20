#Plot Making Age Data Figure 1
fig1 <- function(input) {
  ageD <- input
  
  cut <- c(ageD[6], ageD[11], ageD[16])
  colours <- c("blue3", "orange2", "gray5")
  plot(1, xlab = "", ylab = "", ylim = c(0, 0.00057), xlim = c(1,12), cex.axis = 2.6, tck = 0.02, cex.lab = 2.5)
  l <- 0
  
  for (i in cut) {
    l <- l+1
    
    m <- sapply(levels(factor(i$Age)), function(x) {
      mean(i$Frequency[i$Age==x])
    })
    
    df <- data.frame(Age = as.numeric(names(m)), Frequency = unname(m))
    mod <- nls(Frequency ~ a*Age^b, data = df, start = list(a=1,b=1), control = list(maxiter=1000) )
    points(df,  cex = 3.0, pch=20, col=colours[[l]], add=T)
    a <- mod$m$getPars()[[1]]
    b <- mod$m$getPars()[[2]]
    curve(a*x^b, add = T, lwd = 2.5, col=colours[[l]])
  }

}

## Skeleton Code For plotting
stat <- sapply(1:(ncol(res)-1), function(x) {
  fit <- res[[1,x]]
  sum(fit$growth)
})
plot(head(colnames(res), -1), stat, ylab = "" , xlab= "", cex.axis = 2.6, tck= 0.02, lwd=6.0, type = "l", col = "blue3")

 
ADtn <- readRDS("Tennessee/tn93TnAD.rds")
ADst <- readRDS("Seattle/tn93StAD.rds")
ADna <- readRDS("NAlberta/tn93NAAD.rds")

GDtn <- readRDS("Tennessee/tn93TnGD.rds")
GDst <- readRDS("Seattle/tn93StGD.rds")
GDna <- readRDS("NAlberta/tn93NAGD.rds")

