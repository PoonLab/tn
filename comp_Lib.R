#Supply functions necessary for efficient data storage and processing
#Possibly unnecessary with heavy restructuring (ie. if graph creation can be avoided)
library(dplyr)
library(data.table)
library(parallel)

#Creates a set of data-tables representing a graph of sequences, with the edges between those sequences representing the TN93 Distance.
#The time and location associated with the sequence can be taken either directly from the sequence header, or provided separately in a .csv file
#This set of data tables also includes the set of minimum retrospective edges from sequences at the newest time point.
impTN93 <- function(iFile, reVars="_", varInd=c(1, 2),  addVarN=NA, addVarT=NA, 
                    dateFormat="%Y", partQ=0.95){
  #@param iFile: The name/path of the input file (expecting tn93 output csv)
  #@param reVars: The regular expression used to extract variables from column headers. This is passed to strsplit, creating a vertex of values from the column header
  #@param varInd: A vector of numbers describing the order of variables in the split string. This should describe the index of the unique ID, the Timepoint and the location.
  #               ex. If the header would be split such that the 4th index is the Unique ID, then 4 should be the first number in this list
  #               ID and timepoint are currently required. If the location information is not available, it should be set as "0".
  #@param partQ: The proportion of the set that is to define "known" cases for the purposes of cluster growth. The remaining quantile is marked as new cases.
  #@oaram addvarN: The names of additional variables beyond the second.
  #@param addVarT: The variable types of additional variables beyond the second. Accepts several strings.
  #                "num" for numeric, "str" for characters, "bool" for logical, and "date" for dates (these must be formatted consistantly with other inputted dates) 
  #@return: A list of 3 Data frames. An edge list (weighted by TN93 genetic distance), a vertex list, 
  #         and a list of minimum edges, for the future establishment of a timepoint-based model
  
  #Obtain tn93 edge list from file
  idt <- fread(iFile)
  
  #Reformat edge list as data table object with predictors extracted from sequence header
  temp <- sapply(c(idt$ID1, idt$ID2), function(x) ((strsplit(x,reVars)[[1]]))[varInd])

  el <- data.table(ID1=as.character(temp[1,1:nrow(idt)]), t1=as.Date(temp[2,1:nrow(idt)], format=dateFormat),
                   ID2=as.character(temp[1,(nrow(idt)+1):(2*nrow(idt))]), t2=as.Date(temp[2,(nrow(idt)+1):(2*nrow(idt))], format=dateFormat),
                   Distance = as.numeric(idt$Distance), stringsAsFactors= F)
  
  #Obtain the maximum time and time difference between the head and tail of each edge
  el[,"tMax" := pmax(t1,t2)] 
  el[,"tDiff" := (t1-t2)]
  
  #Obtain list of unique sequences (also data table)
  vl <- data.table(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), stringsAsFactors=F)
  
  #In the event of additional variables
  if(length(varInd)>2){
    
    #A quick type-converter for the ability to input types
    typeConv <- function(x, type) {
      if(type%in%c("num", "numeric")) {return(as.numeric(x))}
      if(type%in%c("str", "string", "char", "character")) {return(as.character(x))}
      if(type%in%c("fac", "factor")) {return(as.factor(x))}
      if(type%in%c("log", "logic", "logical", "bool", "boolean")) {return(as.logical(x))}
      if(type%in%c("date", "time")) {return(as.Date(x, format=dateFormat))}
    }
    
    #Add additional variables
    for(i in 1:(length(varInd)-2)) {
      
      #Add variables named based on character strings in addVarN
      #These variables types are based on addVarT, converted with typeConv()
      v1 <- (paste0(addVarN[i],"1"))
      v2 <- (paste0(addVarN[i],"2"))
      el[,(v1):=typeConv(temp[i+2,1:nrow(idt)], addVarT[i])]
      el[,(v2):=typeConv(temp[i+2,(nrow(idt)+1):(2*nrow(idt))], addVarT[i])]
      
      #Add this variable to the vertex list
      vl[, (addVarN[i]) := unlist(list(el[, get(v1)], el[, get(v2)]))]
      
      #Calculate differences for edges if applicable
      if(addVarT[i]%in%c("num", "numeric", "date", "time")) {
        el[,paste0(addVarN[i],"Diff") := el[,..v1]-el[,..v2]]
      }
      
      #Calculate shared value if applicable
      if(addVarT[i]%in%c("fac", "factor", "char", "string", "str", "character")) {
        el[,paste0(addVarN[i],"Match") := F]
        el[get(v1) == get(v2), paste0(addVarN[i],"Match") := T]
      }
    }
  }
  
  #Collapse length of edge list into vertex list (1 row per sequence)
  vl <- vl[match(unique(vl$ID), vl$ID),]
  
  #Vertices and edges lists together start a graph object
  g <- list(v=vl[order(vl$Time),], e=el[order(el$tMax),])
  g$e[, "I" := .I]
  g$v[, "I" := .I]
  
  #Set a "New Point", such that cases after this point are considered new for the purposes of growth
  newPoint <- quantile(as.numeric(g$v$Time), partQ)
  
  #Split into testing and training partitions
  #Label new vertices as inherently New, new edges are those which contain a new vertex
  g$v[, "New" := F]
  g$e[, "New" := F]
  g$v[as.numeric(Time)>=newPoint, "New" := T]  
  g$e[as.numeric(tMax)>=newPoint, "New" := T]
  
  #This is the complete List of edges, however, edges are filtered by default if they lie between two new vertices
  #In the future, more edges may be filtered based on a distance cutoff requirement
  g$e[,"Filtered" := F]
  
  #Obtain minimum retrospective edges for new cases
  retE <- g$e[(New)&((t1<newPoint)|(t2<newPoint)),]
  iMRE <- sapply(g$v[(New), (ID)], function(id){
    iE <- retE[(ID1%in%id)|(ID2%in%id), ]
    iE[which.min(iE$Distance)[[1]], I] 
  })
  
  #Only the minimum retrospective edges for all new cases remain unfiltered
  g$e[(New), "Filtered" := T]
  g$e[iMRE, "Filtered" := F]
  
  return(g)
}

#Create clusters based on component clustering by some measure of genetic distance
#If a new case clusters with a known case, this is considered growth
compClu <- function(iG, maxD) {
  #@param iG: The inputted graph. Expecting all vertices, but some edges filtered by distance.
  #@param maxD: The maximum tn93 distance. Edges higher than this distance will be filtered out of the graph
  #@return: The inputted graph, annotated with a cluster size summary and case membership in the vertices section
  
  #Because data.tables are being used, this prevents original values being reassigned via pointer
  iG <- copy(iG)
  
  #Filter cases distant edges
  iG$e[(Distance>maxD), "Filtered" := T]
  
  #Simplify the list of unsorted vertices (just id's) and edges (just head and tail id's)
  #These are separated by whether or not they are new
  vid <- iG$v[!(New), (ID)]
  adj <- as.matrix(iG$e[!(Filtered)&!(New),c("ID1","ID2")])
  adjN <- as.matrix(iG$e[!(Filtered)&(New),c("ID1","ID2")])
  
  #Initialize the first cluster name and a column for cluster membership.
  iG$v$Cluster <- vector(mode="numeric", length=nrow(iG$v))

  #The search vertex becomes the first member of the first cluster and is removed from the searchable set of cluster names
  ci <- 1
  srchV <- vid[ci]
  memV <- srchV
  vid <- setdiff(vid, memV)
  growth <- integer(0)
  
  #Assigning Cluster Membership
  repeat {
    
    #Remove edges internal to search query and list outgoing edges
    adj <- adj[!((adj[,"ID1"]%in%srchV) & (adj[,"ID2"]%in%srchV)),,drop=F]
    exE <- adj[((adj[,"ID1"]%in%srchV) | (adj[,"ID2"]%in%srchV)),,drop=F]
    
    #Find all neighboring vertices to the search vertex (or search vertices) through external edges
    #These are then added to the list of member vertices and removed from the list of searchable vertices
    nbV <- setdiff(c(exE[,"ID1"],exE[,"ID2"]), srchV)
    memV <- c(memV, nbV) 
    vid <- setdiff(vid, nbV)
    
    #If there are no more neigbours to the search vertices, the cluster is completed and we reset the search parameters
    if (length(nbV)==0) {
      
      #Update the growth of this cluster
      memVN <- c(adjN[adjN[,"ID1"]%in%memV,"ID2"], adjN[adjN[,"ID2"]%in%memV,"ID1"])
      growth <- c(growth, length(memVN))
      
      #Update the vertices with their cluster membership and update the list of 
      iG$v[(ID)%in%memV, "Cluster" := ci]

      #The end condition, catching the event that there are no vertices to assign to clusters
      if (length(vid)==0) {break}
      
      #Reset search parameters
      ci <- ci+1
      srchV <- vid[1]
      memV <- srchV
      vid <- setdiff(vid, memV)
      
      next
    }
    
    #Remove all edges within the current cluster from the adjacency list
    adj <- adj[!((adj[,"ID1"]%in%srchV) | (adj[,"ID2"]%in%srchV)),,drop=F]
    srchV <- nbV
  }
  
  #Add the overall size of clusters (before calculating their connectivity to new cases).
  #Also add the connectivity to new clusters under the variable "growth"
  iG$c <- list()
  iG$c$Membership <- lapply(1:ci, function(x){iG$v[(Cluster)==x, (ID)]})
  iG$c$Info <- data.table(Old = as.numeric(table(iG$v[!(New), (Cluster)])), New = growth)
  
  return(iG)
}

#Run across a set of several subGraphs created at various filters, analyzing GAIC at each with clusterAnalyze
##- TO-DO: Monitor Thread-safety -##
GAICRun <- function(iG, maxDs=NA, runID=0, nCores=1, addVarInd=NA, plotGAIC=F) {
  #@param iG: Expecting the entire Graph, but in some cases may take a subset  
  #@param maxDs: A list of cutoff thresholds
  #@param runID: An identifier to stash this particular run and compare it to others
  #@param nCores: The number of cores for parallel functionality. -CURRENTLY SOMETIMES ERROR PRONE
  #@param plotGAIC: If true, plots a visual of the key result (GAIC)
  #@param addVarInd: The indices of additional variables to be used.
  #@return: A data table of each runs cluster information.
  #         Both null and proposed model AIC values, as well as the AIC loss ($nullAIC, $modAIC and $GAIC)
  #         The max size, average size and number of singletons ($SizeMax, $MeanSize and $Singletons)
  #         The total growth, largest growth and number of growing clusters ($GrowthTot, $GrowthMax, and $nGrowing)
  #         The ID of the largest cluster and the cluster with the highest growth ($SizeMaxID and $GrowthMaxID)
  #         The effect ratio of mean recency in growing clusters over mean recency in non-growing clusters ($xMag)
  
  #Initialize a set of cutoffs to observe (based on the genetic distance distribution)
  if  (is.na(maxDs)) {
    steps <- head(hist(subset(iG$e, Distance<0.05)$Distance, plot=FALSE)$breaks, -5)
    maxDs <- seq(0 , max(steps), max(steps)/50) 
  }
  
  #Initialize additional variables as the complete set of possible additional vars
  addVarN <- colnames(subG$v[,!c("ID", "Time", "I", "New", "Cluster")])
  if (!is.na(addVarInd)) {
    addVarN <- addVarN[addVarInd]
  } 
  
  #This function runs through severel comparisons of a model weighted by predictors, to a model without those variables
  df <- mclapply(maxDs, function(d) {
    
    #Obtain clusters
    subG <- compClu(iG, d)
    
    #Obtain recency (a sum value of all member tips collection date recency) for each cluster.
    subG$c$Info[, "Recency" := sapply(subG$c$Membership, function(x) {
      mean(as.numeric(subG$v[ID%in%x, (Time)]) - min(as.numeric(subG$v[,(Time)])))
    })]
    
    #Loop through additional variables 
    for(nm in addVarN){
      
      #For categorical data, a string is returned and factorized
      if(typeof(subG$v[, (nm)]) %in% c("integer", "character")) {
        #The most common value is used to represent the value of a given cluster
        subG$c$Info[, (nm) := as.factor(sapply(subG$c$Membership, function(x) {
          t <- table(subG$v[(ID)%in%x, get(nm)])
          names(t[which.max(t)])
        }))]
      }
      
      #For numeric data, an average is returned
      ## - UNTESTED - ##
      if(typeof(subG$v[, (nm)]) %in% c("double")) {
        subG$c$Info[, (nm) := sapply(subG$c$Membership, function(x) {mean(x)})]
      }
    }
    
    #Create proposed model formula dynamically based on number of new dependant variables
    propForm <- as.formula(paste0("New~", paste(colnames(subG$c$Info)[-2], collapse="+")))
    
    #Compares clusters with weights obtained through variables to clusters with even weights
    fit1 <- glm(formula = propForm, data = subG$c$Info, family = "poisson")
    fit2 <- glm(formula = New~Old, data = subG$c$Info, family = "poisson")
    
    #GAIC is the difference between the AIC of two models
    #Put another way, this is the AIC loss associated with predictive variables
    #Other descriptive data characteristics
    dfi <- data.frame(modAIC=fit1$aic, nullAIC=fit2$aic, GAIC=(fit1$aic-fit2$aic),
               GrowthTot=sum(subG$c$Info[,(New)]), Singletons=nrow(subG$c$Info[(Old)==1,]), MeanSize=mean(subG$c$Info[,(Old)]),
               GrowthMax=max(subG$c$Info[,(New)])[[1]], GrowthMaxID=which.max(subG$c$Info[,(New)]), nGrowing=nrow(subG$c$Info[(New)>0,]),
               SizeMax=max(subG$c$Info[,(Old)]), SizeMaxID=which.max(subG$c$Info[,(Old)]))
    
    #To track the potentially many variable effect sizes
    for(cf in names(fit1$coefficients)){
      dfi[,paste0(cf, "_Coefficient")]=fit1$coefficients[cf]
    }
    
    return(dfi)
  }, mc.cores=nCores)
  
  #Convert output to data table and add reference labels based on the run information
  dt <- as.data.table(bind_rows(df))
  dt[,"RunID" := runID]
  dt[,"MaxDistance" := maxDs]
  
  #For quick visual feedback
  if(plotGAIC) {
    plot(dt$MaxDistance, dt$GAIC, xlab="Threshold", ylab="AIC Loss")
    lines(dt$MaxDistance, dt$GAIC)
    abline(h=0, lty=2)
    abline(v=dt$MaxDistance[which.max(abs(dt$GAIC))], lty=2)
  }
  
  return(dt)
}

#Obtain several run results for a much larger set of sub-sampled graphs
multiGAICRun <- function(iG, n, maxDs=NA, prop=0.80) {
  #@param iG: Expecting the entire Graph  
  #@param maxDs: A list of cutoff thresholds
  #@param n: The number of subsample runs to determine
  #@param prop: The proportion of the data set sampled for each sub-sample.
  #@return: A data table of multiple runs worth of cluster information.
  #         Each run will be labelled with a particular run ID

  #Initialize a set of cutoffs to observe (based on the genetic distance distribution)
  if  (is.na(maxDs)) {
    steps <- head(hist(subset(iG$e, Distance<0.05)$Distance, plot=FALSE)$breaks, -5)
    maxDs <- seq(0 , max(steps), max(steps)/50) 
  }
  
  #Run n different runs, each labelled with their own runID
  dt <- bind_rows(lapply(1:n, function(i){
    
    #Copy a sample instance of the inserted graph
    sampG <- copy(iG)
    
    #Sample new and old IDs 
    sampIDs <- c(sample(iG$v[!(New),(ID)], round(prop*nrow(iG$v[!(New)]))),
                 iG$v[(New), (ID)])
    
    sampG$v <- iG$v[(ID)%in%sampIDs,]
    sampG$e <- iG$e[((ID1)%in%sampIDs)&((ID2)%in%sampIDs),]
    
    print(i)
    
    GAICRun(sampG, maxDs, runID=i)
  }))
  
  return(dt)
}
