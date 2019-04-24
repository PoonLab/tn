library(ggplot2)
library(gghighlight)
library(MASS)
library(egg)

#Plot Making Age Data Figure 
expAge <- function(ageD,letter) {
  cuts <- sapply(seq(6,16,5), function(x){
    i <- ageD[[x]]
    sapply(levels(factor(i$Age)), function(x) {
      mean(i$Frequency[i$Age==x])
    })
  })
  df <- data.frame(Age=as.numeric(levels(factor(ageD[[1]]$Age))), pt1=cuts[,1], pt2= cuts[,2], pt3=cuts[,3])
  
  pngTitle <- paste0("fig1", letter,".png")
  png(pngTitle, width=1500, height=1000)
  
  mod1 <- lm(pt1 ~ a*Age^b, data = df, start = list(a=1,b=1), control = list(maxiter=1000) )
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
  return(cuts)
}

res <- GDna
sapply(1:(ncol(res)), function(x) {
  fit <- res[[1,x]]
  sum(fit$growth)
})


####Actual Figs

#Obtains a filtered subgraph of the full graph. Vertices are removed beyond a given year and edges are removed below a cutoff
graphPlot <- function(inG, y, d, col) {
  
  #Removes vertices beyond a current year
  outV <- V(inG)[V(inG)$year>y]
  outG <- inG - outV
  
  #Removes edges with distances above a certain cutoff
  outE <- E(outG)[E(outG)$Distance>=d]
  outG <- outG - outE
  
  #Plot option ignores clusters of size 1 and provides a graph (for ease of overview, not for calculations)
  outG <- subgraph.edges(outG, E(outG), delete.vertices = T)
  plot(outG, vertex.size = 2, vertex.label = NA, vertex.color= col,
       edge.width = 0.65, edge.color = 'black', 
       margin = c(0,0,0,0))
  #sub=paste0(title, " Network, at d=", d),
}

linAge <- function(ageD, letter) {
  cuts <- sapply(seq(6,16,5), function(x){
    i <- ageD[[x]]
    sapply(levels(factor(i$Age)), function(x) {
      mean(i$Frequency[i$Age==x])
    })
  })
  df <- data.frame(Age=as.numeric(levels(factor(ageD[[1]]$Age))), pt1=cuts[,1], pt2= cuts[,2], pt3=cuts[,3])
  
  
  lines <- c("0.010" = "blue", "0.015" = "dark blue", "0.020" = "black")
  p <- ggplot(df, aes(x=Age)) +
    labs(title=letter, x="Time Difference (collection year)", y="Mean of Edge Density in Bipartite Graph") +
    theme(axis.title.x = element_text(size=12, margin=margin(t=10)),
          axis.title.y = element_text(size=12, margin=margin(r=10)), 
          axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=20, hjust=0.5, vjust=-0.1, margin = margin(b=10, t=10)),
          legend.text = element_text(size=12),
          legend.title = element_text(size=15)) +
    geom_point(aes(y=pt1, colour="0.010")) +
    geom_point(aes(y=pt2, colour="0.015")) +
    geom_point(aes(y=pt3, colour="0.020")) +
    geom_smooth(aes(y=pt1, colour="0.010"), method = "lm", se=F) +
    geom_smooth(aes(y=pt2, colour="0.015"), method = "lm", se=F) +
    geom_smooth(aes(y=pt3, colour="0.020"), method = "lm", se=F) +
    scale_colour_manual(name="TN93 Cutoff Threshold", values=lines)
  
  return(p)
}

linGrowth <- function(growthD) {
  
  st <- sapply(1:(ncol(growthD[[1]])), function(x) {
    res <- growthD[[1]]
    fit <- res[[1,x]]
    c(sum(fit$growth), mean(fit$growth))
  })
  
  na <- sapply(1:(ncol(growthD[[2]])), function(x) {
    res <- growthD[[2]]
    fit <- res[[1,x]]
    c(sum(fit$growth), mean(fit$growth))
  }) 
  
  lines <- c("Seattle" = "blue", "North Alberta" = "gold")
  lines2<- c("Mean Growth" = "solid", "Sum Growth"="dashed")
  df <- data.frame(Threshold = seq(0.005,0.05,0.001), stTotalGrowth = st[1,], stMeanGrowth = st[2,],
                   naTotalGrowth = na[1,], naMeanGrowth = na[2,])
  
  ggplot(df, aes(x=Threshold)) +
    labs(title="" ,x= "TN93 Distance Cutoff Threshold", y="Growth of Clusters") +
    theme(axis.title.x = element_text(size=12, margin=margin(t=10)),
          axis.title.y = element_text(size=12), 
          axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=20, hjust=-0.05, vjust=-0.05),
          legend.text = element_text(size=15)) +
    geom_line(aes(y=stMeanGrowth, colour = "Seattle", linetype="Mean Growth"), size=1) +
    geom_line(aes(y=stTotalGrowth, colour = "Seattle", linetype="Sum Growth"), size=1) +
    geom_line(aes(y=naMeanGrowth, colour = "North Alberta", linetype="Mean Growth"), size=1.0) +
    geom_line(aes(y=naTotalGrowth, colour = "North Alberta", linetype="Sum Growth"), size=1.0) +
    scale_colour_manual(name="", values=lines)+
    scale_linetype_manual(name="", values =lines2 )
}

gaicPlot <- function(gaicD) {
  
  df <- data.frame(Threshold = seq(0.005,0.04,0.001), GAIC1 = head(gaicD[[1]], -10), GAIC2= head(gaicD[[2]], -10))
  lines <- c("Seattle" = "blue", "North Alberta" = "gold")
  
  ggplot(df, aes(x=Threshold)) +
    theme(axis.title.x = element_text(size=12, margin=margin(t=10)),
          axis.title.y = element_text(size=12), 
          axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=20, hjust=-0.05, vjust=-0.05),
          legend.text = element_text(size=15)) +
    geom_line(aes(y=GAIC1, colour="Seattle"), size=1.2)+
    geom_line(aes(y=GAIC2, colour="North Alberta"), size=1.2)+
    geom_vline(xintercept = c(df$Threshold[df$GAIC1==min(df$GAIC1)],df$Threshold[df$GAIC2==min(df$GAIC2)]),linetype=4, colour="black", alpha=0.5)+
    labs(title="", x= "TN93 Distance Cutoff Threshold", y="GAIC")+ 
    scale_colour_manual(name="", values=lines)
}

distPlot <- function(inG, letter, col) {
  ed <- data.frame(Distance=E(inG)$Distance)

  ggplot(ed, aes(x=Distance)) +
    labs(title=letter, x="Weight of Edge (TN93 Distance)", y="") +
    theme(axis.title.x = element_text(size=12, margin = margin(t=10)),
          axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=20, hjust=0.5, vjust=-0.1, margin = margin(b=10)))+
    stat_bin(binwidth = 0.001, colour=col, fill="grey") + 
    xlim(0.005,0.05)
}

yearPlot <- function(inG, col) {
  yd <- data.frame(Year=V(inG)$year)
  
  ggplot(yd, aes(x=Year)) +
    labs(title="", y="", x="Year of Vertex (Sequence Collection Year)") +
    theme(axis.title.x = element_text(size=12, margin = margin(t=15)),
          axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=10)+
    scale_x_continuous(breaks = seq(min(V(inG)$year), max(V(inG)$year), 1)))+
    stat_bin(binwidth = 1, colour=col, fill="grey")+
    stat_bin(binwidth = 1, geom="text", aes(label=..count..), vjust=-1)
}

###############
#Linear Update
####################################################################

ggarrange(linAge(ADst, "Seattle"), linAge(ADna, "North Alberta"), 
          nrow = 2, padding=10, labels = c("A", "B"), label.args = list(gp = grid::gpar(font = 1, cex =1.5)))
linGrowth(list(GDst,GDna))
gaicPlot(list(UnWst,UnWna))
ggarrange(distPlot(g1, "Seattle",  "blue"),
          distPlot(g2, "North Alberta", "gold"),
          yearPlot(g1, "blue"),
          yearPlot(g2, "gold"),
          nrow = 2, ncol=2, padding=10, labels = c("A", "B", "C", "D"), label.args = list(gp = grid::gpar(font = 1, cex =1.5)))

par(mfrow=c(1,2))
graphPlot(g1, 2011, 0.013, "blue")
title("Seatte at d=0.013", line=-3)
title("A", line=1, adj=0,cex.main=3)
graphPlot(g2, 2012, 0.011, "orange") 
title("North Alberta data at d=0.011",line=-3)
title("B", line=1, adj=0, cex.main=3)

GDst <- readRDS("pub1/stDGD2.rds")
GDna <- readRDS("pub1/naDGD2.rds")

UnWst <- readRDS("pub1/stDUnW2.rds")
UnWna <- readRDS("pub1/naDUnW2.rds")

ADst <- readRDS("pub1/stDAD.rds")
ADna <- readRDS("pub1/naDAD.rds")
