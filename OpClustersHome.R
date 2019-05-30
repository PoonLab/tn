#Import Libraries
library(igraph,verbose = FALSE)
library(dplyr,verbose = FALSE)
library(parallel,verbose = FALSE)
library(ggplot2,verbose = FALSE)

#Expecting tn93 output as second param
## USAGE: Rscript ~/git/tn/OpClusters.R ___D.txt ##

#Obtain the frequency of edges in a bipartite Graph between two different years as a function of the difference between those years
bpeFreq <- function(iG) {
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A data frame of Number of positives (edges from one year to the newest year) 
  #         with total possible edges and time difference (in years) between the two years
  
  # Obtain the range of years
  maxY <- max(V(iG)$year)
  minY <- min(V(iG)$year)
  ys <- seq(minY, (maxY-1), 1)
  nV <- V(iG)[year==maxY]  # nodes in more recent year of subgraph
  
  # Obtain the frequency of new cases being connected to each year
  frequency <- sapply(ys, function(x) {
    pV <- V(iG)[year==x]  # nodes in older year
    bE <- E(iG)[pV%--%nV]  # bipartite edges
    pos <- length(bE)
    tot <- length(pV)*length(nV)
    return(c(pos,tot))
  })
  
  #Assign age to every case
  tDiff <- sapply(ys, function(x) maxY-x)
  
  #Create a data frame of case attachment frequency and case age
  df <- data.frame(tDiff = tDiff, Positive = frequency[1,], Total = frequency[2,])
  
  return(df)
}

#Filters the input graph such that all new cases are only linked to old cases by their closest edge to old cases
minFilt <- function(iG) {
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A filtered version of this same graph, with all new cases holding only one edge to old cases
  
  #Obtain the new year
  nY <- max(V(iG)$year)
  
  #Obtain edge id's of all of the shortest edge lengths from new cases (A new case can only be linked to 1 case)
  bE <- E(iG)[V(iG)[year==nY]%--%V(iG)[year<nY]]
  
  #To catch a case where no new cases link to old ones
  if (length(bE) > 0) {
    
    #Obtain the closest edges for each new case
    cE <-lapply(V(iG)[year==nY], function(x) {
      xE <- bE[inc(x)]
      
      #To catch a case that is new, but has no linkages to old cases
      ifelse(length(xE)==0, 0, (xE[Distance == min(Distance)])[1] ) 
    })#, mc.cores=8)
    
    #Remove the entries from new cases that dont connect to  old cases
    cE <- unname(unlist(cE[cE!=0]))
    
    #Filter out all edges except for the closest edges
    if(!is.null(cE)){
      iG <- iG - difference(bE, E(iG)[cE])
    }
  }
  
  return(iG)
}

#Obtains the growth of predefined old clusters based on the addition of new clusters
simGrow <- function(iG) { 
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: The cluster information for that subgraph, annotated with the growth of each cluster
  
  #Split the input graph into the new cases and present clusters
  nV <- V(iG)[year==max(V(iG)$year)]
  pG <- induced_subgraph(iG, V(iG)[year<max(V(iG)$year)])
  clu <- components(pG)
  
  #Assign cluster growth based on number of new cases linked to old cases in clusters 
  temp <- sapply(1:clu$no, function(x) {
    members <- names(clu$membership[unname(clu$membership)==x])
    memV <- V(iG)[name%in%members]
    bE <- E(iG)[memV%--%nV]
    forecast <- sum(memV$freq) 
    growth <- length(bE)
    return(c(growth,forecast))
  })
  
  #Assign growth, forecast (based on diagnostic date), and incidence
  clu$growth <- temp[1,]
  clu$forecast <- temp[2,]
  clu$inc <- length(nV)
  
  return(clu)
}

#Plot the GAIC between an informed and uninformed function over a set of thresholds
gaicPlot <- function(growthD,  thresh = cutoffs) {
  #@param growthD: A list of clustering information at various cutoffs, annotated with growth (simGrow, output)
  #@param thresh: A list of cutoff thresholds to representing the independant variable
  #@return: A visual graph of plotted GAIC between two models over the course of @thresh (a list of cutoffs)
  
  #Extract GAIC measurements
  gaicD <- sapply(growthD, function(x) {x$gaic})
  
  #PLace Data into frame
  df <- data.frame(Threshold = thresh, GAIC1 = gaicD)
  min <- df$Threshold[which(df$GAIC1==min(df$GAIC1))[[1]]]
  
  #Generate plot
  ggplot(df, aes(x=Threshold)) +
    theme(axis.title.x = element_text(size=12, margin=margin(t=10)),
          axis.title.y = element_text(size=12), 
          axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=20, hjust=-0.05, vjust=-0.05),
          legend.text = element_text(size=15)) +
    geom_line(aes(y=GAIC1), size=1.2)+
    geom_vline(xintercept = min, linetype=4, colour="black", alpha=0.5)+
    geom_text(aes(min, 5, label = min, vjust =1.5))+
    labs(title="", x= "TN93 Distance Cutoff Threshold", y="GAIC") 
}

## Importing Case data
#____________________________________________________________________________________________________________________________#

#Expecting the output from a tn93 run formatted to a csv file.
#Expecting patient information in the format ID_Date
args = commandArgs(trailingOnly = T)
input <- read.csv(args, stringsAsFactors = F)

#This script will give warnings due to the fact that there are low fit rates on the null model
options(warn=-1)

#Creates a graph based on the inputted data frame. The tn93 Distances become edge4 attributes
g <- graph_from_data_frame(input, directed=F, vertices=NULL)

#Adds the ID's and Sample collection years as different vertex attributes for each vertex
temp <- sapply(V(g)$name, function(x) strsplit(x, '_')[[1]])
V(g)$name <- temp[1,]
V(g)$year <- as.numeric(temp[2,])

#Obtain the range of years and the maximum input year
years <- as.integer(levels(factor(V(g)$year)))
nY <- max(years)
while (length(V(g)[year==nY])<63) {nY <- nY-1}

g <- induced_subgraph(g, V(g)[year<=nY])

years <- as.integer(levels(factor(V(g)$year)))
V(g)$tDiff <- sapply(V(g)$year, function(x) nY-x)

#Initialize a set of cutoffs to observe
steps <- head(hist(E(g)$Distance, plot=FALSE)$breaks,-5)
cutoffs <- seq(0 , max(steps), max(steps)/50)

g <- minFilt(g)

#Create a set of subgraphs at each cutoff
gs <- lapply(cutoffs, function(d) {
  subgraph.edges(g,E(g)[Distance<=d], delete.vertices = F)
})#, mc.cores=8) 
names(gs) <- cutoffs

## Generate Growth data
#__________________________________________________________________________________________________________________________#

res <- mclapply(cutoffs, function(d) {
  cat(paste0("\r", "Running Analysis ", d/max(cutoffs)*100, "%"))
  #Obtain a subGraph at the maximum year, removing edges above the distance cutoff and ensuring no merging by removing, non-closest edges to new cases
  subG <- gs[[as.character(d)]]
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case ag
  #This data may contain missing cases, hense the complete cases addition
  ageDi <- bind_rows(lapply(rev(tail(years,-2)), function(y){
    ssubG <- minFilt(induced_subgraph(subG, V(subG)[year<y]))
    bpeFreq(subG)
  }))
  
  mod <- glm(cbind(Positive, Total) ~ tDiff, data=ageDi, family='binomial')
  
  #Assign a predicted growth value to each member of the graph
  V(subG)$freq <- predict(mod, data.frame(tDiff=V(subG)$tDiff), type='response')
  
  #Obtain growth based on two models restricted model
  clu <- simGrow(subG)
  
  #Place growth and forecast data in dfs for fit and full growth
  df1 <- data.frame(Growth = clu$growth, Pred = clu$forecast)
  df2 <- data.frame(Growth = clu$growth, Pred = clu$csize * (sum(clu$growth)/sum(clu$csize)))
  
  #Model growth as a function of forecast for fit and full growth models
  mod1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
  mod2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
  
  #Save, gaic, model and age data
  clu$gaic <- mod1$aic-mod2$aic
  clu$mod <- mod
  clu$ageD <- ageDi
  
  return(clu)
}, mc.cores=8)

#Label data
names(res) <- cutoffs

## Generate Pictures and output
#__________________________________________________________________________________________________________________________#

#Obtain Minimum GAIC estemating cutoff threshold and the network associated with it
gaics <- sapply(res, function(x) {x$gaic})
do <- names(which(gaics==min(gaics))[1])
opt <- gs[[do]]

#Plot option ignores clusters of size 1 and provides a graph (for ease of overview, not for calculations)
optPG <- subgraph.edges(opt, E(opt), delete.vertices = T)

#Create output pdf
pdf(file = paste0(gsub("\\..*", "", args), "VS.pdf"))

#Plot GAIC
gaicPlot(res)

#Plot Network
plot(optPG, vertex.size = 2, vertex.label = NA, vertex.color= "orange",
     edge.width = 0.65, edge.color = 'black', 
     margin = c(0,0,0,0))

dev.off()

#Obtain the information from opt cluster and print it to stOut
optClu <- components(opt)
optClu$years <- table(V(g)$year)
optClu$no <- NULL
optClu$csize <- sort(table(optClu$membership)[table(optClu$membership)>1], decreasing =T)
print(optClu)

#Save all growth data in accessable files
saveRDS(res, file = paste0(gsub("\\..*", "", args), "GD.rds"))

cat(paste0("\n","Done" ))