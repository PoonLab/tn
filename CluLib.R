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
minFilt <- function(iG, Home=F) {
  #@param iG: A subGraph cut based on a threshold distance, with the latest cases representing New cases (ie. Upcoming cases)
  #@return: A filtered version of this same graph, with all new cases holding only one edge to old cases
  
  #Obtain the new year
  nY <- max(V(iG)$year)
  
  #Obtain edge id's of all of the shortest edge lengths from new cases (A new case can only be linked to 1 case)
  bE <- E(iG)[V(iG)[year==nY]%--%V(iG)[year<nY]]
  
  #To catch a case where no new cases link to old ones
  if (length(bE) > 0) {
    
    #Resolve Parallel problem
    if (Home) {
      #Obtain the closest edges for each new case
      cE <- lapply(V(iG)[year==nY], function(x) {
        xE <- bE[inc(x)]
        
        #To catch a case that is new, but has no linkages to old cases
        ifelse(length(xE)==0, 0, (xE[Distance == min(Distance)])[1] ) 
      }, mc.cores=8)
    
    } else {
      #Obtain the closest edges for each new case
      cE <-mclapply(V(iG)[year==nY], function(x) {
        xE <- bE[inc(x)]
        
        #To catch a case that is new, but has no linkages to old cases
        ifelse(length(xE)==0, 0, (xE[Distance == min(Distance)])[1] ) 
      }, mc.cores=8)
    }
    
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


clusterAnalyze <- function(subG) {
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
}

createGraphSet <- function(runArgs){
  infile <- runArgs[1]
  outfile <- ifelse(exists(runArgs[2]), runArgs[2], infile)
  inputFilter <- as.numeric(runArgs[3])
  
  input <- read.csv(infile, stringsAsFactors = F)
  
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
  while (length(V(g)[year==nY])<63 || inputFilter>0) {
    nY <- nY-1
    if (length(V(g)[year==nY])>63){inputFilter <- inputFilter-1}
  }
  g <- induced_subgraph(g, V(g)[year<=nY])
  
  years <- as.integer(levels(factor(V(g)$year)))
  V(g)$tDiff <- sapply(V(g)$year, function(x) nY-x)
  
  #Initialize a set of cutoffs to observe
  steps <- head(hist(E(g)$Distance, plot=FALSE)$breaks,-5)
  cutoffs <- seq(0 , max(steps), max(steps)/50)
  
  #Resolve potential for merging clusters
  g <- minFilt(g)
  
  #Create a set of subgraphs at each cutoff
  gs <- mclapply(cutoffs, function(d) {
    subgraph.edges(g,E(g)[Distance<=d], delete.vertices = F)
  }, mc.cores=8) 
  names(gs) <- cutoffs
  
  return(gs)
}