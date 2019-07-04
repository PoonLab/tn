#Creates a graph based on some inputted run arguments and potentially patient meta-data
createGraph <- function(infile, inputFilter, metData){
  #@param infile: The name/path of the input file (expecting tn93 output)
  #@param inputFilter: Will drop x of the most recent years from the total data set based on this input
  #@param metData: the filename for a dataframe of associated metadata (optional).
  #@return: A graph with each vertex having a name and an id as well as every edge representing some measure of distance.
  
  #From the input file, a tn93 output file. This will be an edgeList
  input <- read.csv(infile, stringsAsFactors = F)
  temp1 <- sapply(eL$ID1, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp2 <- sapply(eL$ID1, function(x) (strsplit(x,'_')[[1]])[[2]])
  temp3 <- sapply(eL$ID2, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp4 <- sapply(eL$ID2, function(x) (strsplit(x,'_')[[1]])[[2]])
  
  #Represents edge data and vertex data as two seperate dataframes
  el <- data.frame(ID1=as.character(temp1), t1=as.numeric(temp2), ID2=as.character(temp3), t2=as.numeric(temp4), 
                   Distance = as.numeric(eL$Distance), stringsAsFactors= F)
  el$tDiff <- abs(el$t1 - el$t2)
  vl <- unique(data.frame(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), stringsAsFactors=F))
  
  ##TO-DO: Meta-Data (potentially broken in previous code before)
  
  #Arrange into a named list representing the graph
  g <- list(v=vl, e=el)
  
  return(g)
}

clusters <- function(inG) {
  
  temp <- inG$v[,"ID"]
  search <- temp[1]
  adj <- inG$e[,c("ID1","ID2")]
  i <- new <- 1
  clu <- list()
  
  while ( (nrow(temp)>0) && (nrow(adj)>0) ) {
    neighbours <- unique(c(adj$ID1[which(adj$ID2%in%search)], adj$ID2[which(adj$ID1%in%search)]))
    new <- neigbours - search
    search <- neighbours
  }
  
}