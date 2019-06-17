library(ggplot2)
library(gghighlight)
library(MASS)
library(egg)
library(igraph)

#Plot Making tDiff Data Figure 
exptDiff <- function(ageD,letter) {
  cuts <- sapply(seq(6,16,5), function(x){
    i <- ageD[[x]]
    sapply(levels(factor(i$tDiff)), function(x) {
      mean(i$Frequency[i$tDiff==x])
    })
  })
  df <- data.frame(tDiff=as.numeric(levels(factor(ageD[[1]]$tDiff))), pt1=cuts[,1], pt2= cuts[,2], pt3=cuts[,3])
  
  pngTitle <- paste0("fig1", letter,".png")
  png(pngTitle, width=1500, height=1000)
  
  mod1 <- lm(pt1 ~ a*tDiff^b, data = df, start = list(a=1,b=1), control = list(maxiter=1000) )
  a1 <- mod1$m$getPars()[[1]]
  b1 <- mod1$m$getPars()[[2]]
  df$mod1 <- a1*df$tDiff^b1
  
  mod2 <- nls(pt2 ~ a*tDiff^b, data = df, start = list(a=1,b=1), control = list(maxiter=1000) )
  a2 <- mod2$m$getPars()[[1]]
  b2 <- mod2$m$getPars()[[2]]
  df$mod2 <- a2*df$tDiff^b2
  
  mod3 <- nls(pt3 ~ a*tDiff^b, data = df, start = list(a=1,b=1), control = list(maxiter=1000) )
  a3 <- mod3$m$getPars()[[1]]
  b3 <- mod3$m$getPars()[[2]]
  df$mod3 <- a3*df$tDiff^b3
  
  lines <- c("0.010" = "blue", "0.015" = "black", "0.020" = "orange")
  p <- ggplot(df, aes(x=tDiff)) +
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
    m <- sapply(levels(factor(i$tDiff)), function(x) {mean(i$Frequency[i$tDiff==x])}) 
    fit <- fitdistr(m, "exponential")
    test <- ks.test(m, "pexp", fit$estimate)
    return(test[[2]])
  })
  return(cuts)
}

ADfitbn <- function(ageD) {
  sapply(ageD, function(ageDi) {
    mod <- 
    mod$coefficients[[2]]
  })
}

edgeFreq <- function(ageD){
  cuts <- sapply(ageD, function(i){
    sapply(levels(factor(i$tDiff)), function(x) {
      mean(i$Frequency[i$tDiff==x])
    })
  })
  return(cuts)
}

sapply(GDna, function(x) {
  mod <- x$ageD
  1-mod$deviance/mod$null.deviance
})

sapply(GDna, function(x) {
  x$gaic
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

lintDiff <- function(ageD, letter) {
  cuts <- sapply(seq(1,16,5), function(x){
    ageDi <- ageD[[x]]
    mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')
    predict(mod, data.frame(tDiff=seq(1,12,1)), type='response')
    })
  
  df <- data.frame(tDiff = seq(1,12,1), pt1 = cuts[,1], pt2 = cuts[,2], pt3 = cuts[,3], pt4 = cuts[,4] )
  lines <- c("0.005"= "royalblue", "0.010" = "blue", "0.015" = "dark blue", "0.020" = "black")
  
  p <- ggplot(df, aes(x=tDiff)) +
    labs(title=letter, x="Time Difference (collection year)", y="Mean of Edge Density in Bipartite Graph") +
    theme(axis.title.x = element_text(size=12, margin=margin(t=10)),
          axis.title.y = element_text(size=12, margin=margin(r=10)), 
          axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=20, hjust=0.5, vjust=-0.1, margin = margin(b=10, t=10)),
          legend.text = element_text(size=12),
          legend.title = element_text(size=15)) +
    geom_line(aes(y=pt1, colour="0.005")) +
    geom_line(aes(y=pt2, colour="0.010")) +
    geom_line(aes(y=pt3, colour="0.015")) +
    geom_line(aes(y=pt4, colour="0.020")) +
    scale_colour_manual(name="TN93 Cutoff Threshold", values=lines)
}

linGrowth <- function(growthD) {
  
  st <- sapply(growthD[[1]], function(x) {
    c(sum(x$growth), mean(x$growth))
  })
  
  na <- sapply(growthD[[2]], function(x) {
    c(sum(x$growth), mean(x$growth))
  })
  
  lines <- c("Seattle" = "blue", "North Alberta" = "orangered3")
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

gaicPlot <- function(growthD) {
  
  st <- sapply(growthD[[1]], function(x) {x$gaic})
  na <- sapply(growthD[[2]], function(x) {x$gaic})
  
  
  df <- data.frame(Threshold = seq(0.005,0.04,0.001), GAIC1 = head(st, -10), GAIC2= head(na, -10))
  lines <- c("Seattle" = "blue", "North Alberta" = "orangered")
  
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
    geom_text(aes(df$Threshold[df$GAIC1==min(df$GAIC1)],5,label = df$Threshold[df$GAIC1==min(df$GAIC1)],vjust =1.5))+
    geom_text(aes(df$Threshold[df$GAIC2==min(df$GAIC2)],5,label = df$Threshold[df$GAIC2==min(df$GAIC2)],vjust =1.5))+
    labs(title="", x= "TN93 Distance Cutoff Threshold", y="GAIC")+ 
    scale_colour_manual(name="", values=lines)
}

distPlot <- function(inG1,inG2) {
  h1 <- hist(E(inG1)$Distance, plot=F)
  h2 <- hist(E(inG2)$Distance, plot=F, breaks=h1$breaks)
  
  df <- data.frame(st = h1$counts/choose(length(V(inG1)),2), na = h2$counts/choose(length(V(inG2)),2))
  
  lines <- c("Seattle" = "blue", "North Alberta" = "orangered3")
  
  ggplot(df, aes(x=head(h1$breaks, -1))) +
    labs(title="", x="Weight of Edge (TN93 Distance)", y="Frequency") +
    theme(axis.title.x = element_text(size=12, margin = margin(t=10)),
          axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=20, hjust=0.5, vjust=-0.1, margin = margin(b=10)))+
    geom_bar(aes(y=st, fill="Seattle"), stat="identity", alpha=0.5) + 
    geom_bar(aes(y=na, fill="North Alberta"), stat="identity",  alpha=0.5 ) + 
    scale_fill_manual(name="", values=lines) 
}

yearPlot <- function(inG1, inG2) {
  h1 <- unname(table(V(inG1)$year))
  h2 <- unname(table(V(inG2)$year))
  
  df <- data.frame(st = c(h1,0), na = c(0,0,0,0,0,0,0, h2) )
  
  lines <- c("Seattle" = "blue", "North Alberta" = "orangered3")
  
  ggplot(df, aes(x=seq(2000,2013))) +
    labs(title="", y="Frequency", x="Year of Vertex (Sequence Collection Year)") +
    theme(axis.title.x = element_text(size=12, margin = margin(t=10)),
          axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=20, hjust=0.5, vjust=-0.1, margin = margin(b=10)))+
    geom_bar(aes(y=st, fill="Seattle"), stat="identity", alpha=0.5) + 
    geom_bar(aes(y=na, fill="North Alberta"), stat="identity",  alpha=0.5 ) + 
    scale_fill_manual(name="", values=lines) 
}

###############
#Linear Update
####################################################################

GDst <- readRDS("stDGD.rds")
GDna <- readRDS("naDGD.rds")

ADst <- readRDS("stDAD.rds")
ADna <- readRDS("naDAD.rds")

ggarrange(lintDiff(ADst, "Seattle"), lintDiff(ADna, "North Alberta"), 
          nrow = 2, padding=10, labels = c("A", "B"), 
          label.args = list(gp = grid::gpar(font = 1, cex =1.5)))
linGrowth(list(GDst,GDna))
gaicPlot(list(GDst,GDna))

# where do g1 and g2 come from?


ggarrange(distPlot(g1, g2),
          yearPlot(g1, g2),
          nrow = 2, padding=10, labels = c("A", "B"), 
          label.args = list(gp = grid::gpar(font = 1, cex =1.5)))

par(mfrow=c(1,2))
graphPlot(g1, 2011, 0.013, "dodgerblue")
title("Seattle", line=-3)
#title("A", line=1, adj=0,cex.main=3)
graphPlot(g2, 2012, 0.011, "orange2") 

title("North Alberta data at d=0.011",line=-3)
#title("B", line=1, adj=0, cex.main=3)

