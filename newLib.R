infile <- "stD.txt"
#infile <- "Data/Seattle/tn93StsubB.txt"
inputFilter <- 0
metData <- NA

#Creates a graph based on some inputted run arguments and potentially patient meta-data
createGraph <- function(infile, inputFilter, metData){
  #@param infile: The name/path of the input file (expecting tn93 output)
  #@param inputFilter: Will drop x of the most recent years from the total data set based on this input
  #@param metData: the filename for a dataframe of associated metadata (optional).
  #@return: A graph with each vertex having a name and an id as well as every edge representing some measure of distance.
  
  #From the input file, a tn93 output file. This will be an edgeList
  input <- read.csv(infile, stringsAsFactors = F)
  temp1 <- sapply(input$ID1, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp2 <- sapply(input$ID1, function(x) (strsplit(x,'_')[[1]])[[2]])
  temp3 <- sapply(input$ID2, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp4 <- sapply(input$ID2, function(x) (strsplit(x,'_')[[1]])[[2]])
  
  #Represents edge data and vertex data as two seperate dataframes
  el <- data.frame(ID1=as.character(temp1), t1=as.numeric(temp2), ID2=as.character(temp3), t2=as.numeric(temp4), 
                   Distance = as.numeric(input$Distance), stringsAsFactors= F)
  el$tDiff <- abs(el$t1 - el$t2)
  vl <- unique(data.frame(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), stringsAsFactors=F))
  
  ##TO-DO: Meta-Data (potentially broken in previous code before)
  
  #Arrange into a named list representing the graph
  g <- list(v=vl, e=el)
  
  return(g)
}

clusters <- function(inG) {
  
  #inG <- eFilt(g, 0.015)
  
  #To set up a column in vector info of cluster membership
  inG$v$Cluster <- vector(mode="numeric", length=nrow(inG$v))
  
  #Simplify the list of vertices (just id's) and edges (just head and tail id's)
  vid <- inG$v[,"ID"]
  adj <- inG$e[,c("ID1","ID2")]
  
  #Initialize the first cluster name and the list of clusters, with c0 reserved for all singletons
  i <- 1
  clu <- list("c0"=vector())
  
  #Initialize the search term, our first vertex to base the clustering off of. 
  #This also becomes the first member of this cluster
  search <- vid[1]
  member <- search
  
  #Because the first search vertex has been sorted into a cluster, we remove it from temp 
  vid <- setdiff(vid, member)
  repeat {

    #Find vertices from the search set and to the search set but not within the search set
    fromSearch <- which(adj$ID1%in%search)
    toSearch <- which(adj$ID2%in%search)
    withinSearch <- intersect(fromSearch, toSearch)
    fromSearch <- setdiff(fromSearch,withinSearch)
    toSearch <- setdiff(toSearch,withinSearch)

    #Find all neighbouring vertices to the search vertex (or search vertices)
    #These are added to the current cluster and removed from temp
    neighbours <- unique(c(adj$ID2[fromSearch],adj$ID1[toSearch]))
    member <- c(member, neighbours) 

    vid <- setdiff(vid, member)

    #If there are no more neigbours to the search vertex, the cluster is completed and we reset the search parameters
    if (length(neighbours)==0) {
      
      #To catch a singleton
      if (length(member)==1) {
        clu[["c0"]] <- c(clu[["c0"]], member)
      }
      else {
        clu[[paste0("c",i)]] <- member
        inG$v$Cluster[which(inG$v$ID%in%member)] <- i
        i <- i+1
      }
      
      #The end condition, catching the event that there are no vertices to assign to clusters
      if (length(vid)==0) {break}
      
      #Reset search parameters
      search <- vid[1]
      member <- search
      vid <- vid[-which(vid%in%search)]

      next
    }
    
    #Remove all edges within the current cluster from the adjacency list
    adj <- adj[-c(withinSearch,fromSearch,toSearch),]
    search <- neighbours
  }

  #Add some summary information regarding clusters
  inG$cSum <- sapply(tail(clu,-1), function(x){length(x)})
  inG$cNo <- length(clu)-1

  outG <- inG
  
  return(list(clu, outG))
}

simGrow <- function(inG, maxD, maxT) {
  
  inG <- tFilt(inG, maxT)
  inG <- dFilt(inG, maxD)
  inG <- clsFilt(inG)
  
  oldG <- tFilt(inG, (maxT-1))
  
  oClu <- clusters(oldG)[[2]]
  nClu <- clusters(inG)[[2]]
  
  oldG <- oClu[[2]]
  inG <- nClu[[2]]
  
  
}

#Remove edges from some graph that sit above a maximum reporting distance.
dFilt <- function(inG, maxD) {
  #@param inG: The input Graph with all edges present
  #@param maxD: The maximum distance, edges with distance above this are filtered out
  #@return: The input Graph with edges above a maximum reporting distance filtered
  
  inG$e <- inG$e[which(inG$e$Distance<maxD),]
  outG <- inG
  return(outG)
}

#Remove vertices from some graph that sit above a maximum time
tFilt <- function(inG, maxT) {
  #@param inG: The input Graph with all vertices present
  #@param maxT: The maximum distance, edges with Time above this are filtered out
  #@return: The input Graph with vertices above a maximum time filtered out
  inG$v <- inG$v[which(inG$v$Time<=maxT),]
  
  #To remove associated edges
  inG$e <- inG$e[-which(inG$e$t1>maxT),]
  inG$e <- inG$e[-which(inG$e$t2>maxT),]
  outG <- inG
  
  return(outG)
}

#Filter the edges coming from new cases such that the new cases have no edges to eachother and only one edge leading from them to old cases
#This is a simplification to help resolve the merging involved in cluster growth and to prevent growth by whole clusters.
clsFilt <- function(inG){
  nV <- inG$v$ID[which(inG$v$Time == max(inG$v$Time))]
  fromNew <- which(inG$e$ID1%in%nV)
  toNew <- which(inG$e$ID2%in%nV)
  withinNew <- intersect(fromNew,toNew)
  
  inG$e <- inG$e[-withinNew,]
  fromNew <- which(inG$e$ID1%in%nV)
  toNew <- which(inG$e$ID2%in%nV)
  
  minE <- sapply(nV, function(v) {
    fromV <- which(inG$e$ID1%in%v)
    toV <- which(inG$e$ID2%in%v)
    incV <- c(fromV,toV)
    
    if (length(incV) > 0){
      vE <- inG$e[incV,]
      minE <- which(vE$Distance == min(vE$Distance))[[1]]
      
      index <- c(fromV,toV)[[minE]]
    } else{
      index <- 0
    }
    
    return(index)
  })
    
  minE <- unname(minE[minE>0])
  remE <- setdiff(as.numeric(c(fromNew, toNew)), minE)
  
  outG <- inG$e[-remE,]

  return(outG)
} 

