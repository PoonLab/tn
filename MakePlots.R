library(ggplot2)
library(MASS)

#Plot Making Age Data Figure 1
fig1 <- function(ageD,letter) {
  cuts <- sapply(seq(6,16,5), function(x){
    i <- ageD[[x]]
    sapply(levels(factor(i$Age)), function(x) {
      mean(i$Frequency[i$Age==x])
    })
  })
  df <- data.frame(Age=as.numeric(levels(factor(ageD[[1]]$Age))), pt1=cuts[,1], pt2= cuts[,2], pt3=cuts[,3])
  
  pngTitle <- paste0("fig1", letter,".png")
  png(pngTitle, width=1500, height=1000)
  
  mod1 <- nls(pt1 ~ a*Age^b, data = df, start = list(a=1,b=1), control = list(maxiter=1000) )
  a1 <- mod1$m$getPars()[[1]]
  b1 <- mod1$m$getPars()[[2]]
  df$mod1 <- a1*df$Age^b1
  
  mod2 <- nls(pt2 ~ a*Age^b, data = df, start = list(a=1,b=1), control = list(maxiter=1000) )
  a2 <- mod2$m$getPars()[[1]]
  b2 <- mod2$m$getPars()[[2]]
  df$mod2 <- a2*df$Age^b2
  
  mod3 <- nls(pt3 ~ a*Age^b, data = df, start = list(a=1,b=1), control = list(maxiter=1000) )
  a3 <- mod3$m$getPars()[[1]]
  b3 <- mod3$m$getPars()[[2]]
  df$mod3 <- a3*df$Age^b3
  
  lines <- c("0.010" = "blue", "0.015" = "black", "0.020" = "orange")
  p <- ggplot(df, aes(x=Age)) +
    labs(title=letter, x="Time Between Case Sample Collection", y="Mean Frequency of Bipartite Edges") +
    theme(axis.title.x = element_text(size=20, margin=margin(t=20)),
          axis.title.y = element_text(size=20, margin=margin(r=20)), 
          axis.text.x = element_text(size=20), 
          axis.text.y = element_text(size=20),
          plot.title = element_text(size=35),
          legend.text = element_text(size=25),
          legend.title = element_text(size=30)) +
    geom_point(aes(y=pt1, colour="0.010")) +
    geom_point(aes(y=pt2, colour="0.015")) +
    geom_point(aes(y=pt3, colour="0.020")) +
    geom_smooth(aes(y=mod1, colour="0.010"), method="lm", formula=y~exp(-x), se=F) +
    geom_smooth(aes(y=mod2, colour="0.015"), method="lm", formula=y~exp(-x), se=F) +
    geom_smooth(aes(y=mod3, colour="0.020"), method="lm", formula=y~exp(-x), se=F) +
    scale_colour_manual(name="TN93 Cutoff Threshold", values=lines)
  
  print(p)
  dev.off()
}

fig2 <- function(res, letter) {C
  
  stat <- sapply(1:(ncol(res)), function(x) {
    fit <- res[[1,x]]
    c(sum(fit$growth), mean(fit$growth))
  })
  lines <- c("Mean Growth" = "blue", "Total Growth" = "orange")
  df <- data.frame(Threshold = seq(0.005,0.05,0.001), TotalGrowth = stat[1,], MeanGrowth = stat[2,])
  
  pngTitle <- paste0("fig2", letter,".png")
  png(pngTitle, width=1500, height=1000)
  
  p <- ggplot(df, aes(x=Threshold)) +
    labs(title=letter ,x= "TN93 Distance Cutoff Threshold", y="Growth (Number of New Cases Linked to Old Cluster Members)") +
    theme(axis.title.x = element_text(size=20, margin=margin(t=20)),
          axis.title.y = element_text(size=20, margin=margin(r=20)), 
          axis.text.x = element_text(size=20), 
          axis.text.y = element_text(size=20),
          plot.title = element_text(size=35),
          legend.text = element_text(size=25)) +
    geom_line(aes(y=MeanGrowth, colour = "Mean Growth"), size=2) +
    geom_point(aes(y=MeanGrowth, colour = "Mean Growth"), size=4, shape=21, stroke=1.5,fill="white") +
    geom_line(aes(y=TotalGrowth, colour = "Total Growth"), size=2) +
    geom_point(aes(y=TotalGrowth, colour = "Total Growth"), size=4, shape=21, stroke=1.5, fill="white") +
    scale_colour_manual(name="", values=lines)
  
  print(p)
  dev.off()
}

fig3 <- function(gaicD, letter) {
  
  df <- data.frame(Threshold = seq(0.005,0.05,0.001), GAIC = gaicD)
                   
  pngTitle <- paste0("fig3", letter,".png")
  png(pngTitle, width=1500, height=1000)
  
  p <- ggplot(df, aes(x=Threshold, y=GAIC)) +
    theme(axis.title.x = element_text(size=20, margin=margin(t=20)),
          axis.title.y = element_text(size=20, margin=margin(r=20)), 
          axis.text.x = element_text(size=20), 
          axis.text.y = element_text(size=20),
          plot.title = element_text(size=35),
          legend.text = element_text(size=25)) +
    geom_line()+
    geom_point()+
    labs(title=letter, x= "TN93 Distance Cutoff Threshold", y="GAIC") 
    
  print(p)
  dev.off()
}

ADfit2 <- function(ageD) {
  cuts <- sapply(ageD, function(i) {
    m <- sapply(levels(factor(i$Age)), function(x) {mean(i$Frequency[i$Age==x])}) 
    fit <- fitdistr(m, "exponential")
    test <- ks.test(m, "pexp", fit$estimate)
    return(test[[2]])
  })
  return(cuts)
}

ADfit1 <- function(ageD) {
  p <- sapply(ageD, function(i) {
    df <- data.frame(Age = as.numeric(levels(factor(i$Age))),  
                     Frequency = sapply(levels(factor(i$Age)), function(x) {mean(i$Frequency[i$Age==x])}), 
                     Null = rep(mean(i$Frequency), length(levels(factor(i$Age)))))
    lm(Frequency ~ Age, data=df)
  })
  return(p)
}

edgeFreq <- function(ageD){
  cuts <- sapply(ageD, function(i){
    sapply(levels(factor(i$Age)), function(x) {
      mean(i$Frequency[i$Age==x])
    })
  })
  return(meancuts)
}

justGrowth <- function(res){
  stat <- sapply(1:(ncol(res)), function(x) {
    fit <- res[[1,x]]
    mean(fit$growth)
  })
  return(stat)
}

ADtn <- readRDS("ColDateData/tnDAD.rds")
ADst <- readRDS("ColDateData/stDAD.rds")
ADna <- readRDS("ColDateData/naDAD.rds")
ADna_dates <- readRDS("ColDateData/naDDAD.rds")

GDtn <- readRDS("ColDateData/tnDGD.rds")
GDst <- readRDS("ColDateData/stDGD.rds")
GDna <- readRDS("ColDateData/naDGD.rds")
GDna_dates <- readRDS("ColDateData/naDDGD.rds")

DisAtn <- readRDS("ColDateData/tnDDisA.rds")
DisAst <- readRDS("ColDateData/stDDisA.rds")
DisAna <- readRDS("ColDateData/naDDisA.rds")
DisAna_dates <- readRDS("ColDateData/naDDUnW.rds")

UnWtn <- readRDS("ColDateData/tnDUnW.rds")
UnWst <- readRDS("ColDateData/stDUnW.rds")
UnWna <- readRDS("ColDateData/naDUnW.rds")
UnWna_dates <- readRDS("ColDateData/naDDDisA.rds")

fig1(ADtn, "A")
fig1(ADst, "B")
fig1(ADna, "C")
fig1(ADna_dates, "")

fig2(GDtn, "A")
fig2(GDst, "B")
fig2(GDna, "C")
fig2(GDna_dates, "")

fig3(DisAtn, "A")
fig3(DisAst, "B")
fig3(DisAna, "C")
fig3(DisAna_dates, "")

fig3(UnWtn, "A")
fig3(UnWst, "B")
fig3(UnWna, "C")
fig3(UnWna_dates, "")

###############
#Linear Update
####################################################################

GDst <- readRDS("pub1/stDGD2.rds")
GDna <- readRDS("pub1/naDGD2.rds")
GDna_dates <- readRDS("pub1/naDDGD2.rds")

UnWst <- readRDS("pub1/stDUnW2.rds")
UnWna <- readRDS("pub1/naDUnW2.rds")
UnWna_dates <- readRDS("pub1/naDDUnW2.rds")

ADst <- readRDS("pub1/stDAD.rds")
ADna <- readRDS("pub1/naDAD.rds")
ADna_dates <- readRDS("pub1/naDDAD.rds")

mean(ADfit1(ADst))
mean(ADfit1(ADna))
mean(ADfit1(ADna_dates))

sd(ADfit1(ADst))
sd(ADfit1(ADna))
sd(ADfit1(ADna_dates))



plot(ADfit2(ADst))
plot(ADfit2(ADna))
plot(ADfit2(ADna_dates))

GD <- justGrowth(GDst)
GD <- head(GD,-1)
df <- data.frame(Cutoff = seq(0.005,0.049, 0.001), Growth = GD)
fit <- nls(Growth ~ a*Cutoff^b, data = df, start = list(a=1,b=1))

justGrowth(GDna)
justGrowth(GDna_dates)