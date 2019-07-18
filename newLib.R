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
      }
      
      #The end condition, catching the event that there are no vertices to assign to clusters
      if (length(vid)==0) {break}
      
      #Reset search parameters
      i <- i+1
      search <- vid[1]
      member <- search
      vid <- vid[-which(vid%in%search)]

      next
    }
    
    #Remove all edges within the current cluster from the adjacency list
    adj <- adj[-c(withinSearch,fromSearch,toSearch),]
    search <- neighbours
  }
  
  
  return(clu)
}

simGrow <- function(inG, maxD) {
  
  #Remove edges above the maximum reporting distance
  inG$e <- inG$e[which(inG$e["Distance"]<maxD),]
  
  #
  maxY <- max(c(inG$e[["t1"]],inG$e[["t2"]]))
  fromNew <- which(inG$e[["t1"]]==maxY)
  toNew <- which(inG$e[["t2"]]==maxY)
  
  growth <- inG$e[c(toNew,fromNew),]
  el <- inG$e[-c(toNew,fromNew),]
  vl <- inG$v[which(inG$v["Time"]<maxY),]
  cluG <- list(v=vl, e=el)

     
  

}

