#' Create an implementation of a graph. This consists of an edge matrix (see ape::dist.dna()) and some sequence meta data
#' A large part of this process involves resolving growth, ensuring that new sequences are only added prospectively without merging clusters.
#' At this point, we also annotate the minimum retrospective edges of each edge. This is stored in sequence information
#' NOTE: A subset of these sequences are expected to be new .
create.graph <- function(seq.info, edge.info, growth.resolution = minimum.retrospective.edge) {
  #'@param seq.info: A set of sequence meta-data sorted by alignment header
  #'@param edge.info: A pairwise edge matrix of all associated headers in seq.info
  #'@param growth.resolution: The method by which growth is resolved. This ensures new cases don't merge clusters
  #'By default, each new sequence joins a cluster byonly it's minimum retrospective
  #'@return: A graph, with sequences and edge info. New sequences are only linked by their minimum retrospective edge
  
  #Check inputs
  if(!all(colnames(edge.info)%in%colnames(edge.info))){
    stop("The pairwise distance matrix does not contain the recognized headers")
  }
  if(!("New"%in%colnames(seq.info))){
    warning("No new sequences are specified by a New column in seq.info.")
    seq.info[,"New":=F]
  }

  #Assemble graph object
  g <- list()
  g$seq.info <- seq.info
  g$edge.info <- edge.info
  g$growth.resolved <- growth.resolution(g)
  
  return(g)
}


#' A growth resolution. New sequences only join old clusters through their minimum, retrospective edge
minimum.retrospective.edge <- function(g) {
  #' @param g: The input graph
  #' @return: A data table matching each new sequence with its closest retrospective neighbour

  #Find the minimum retrospective edge of each sequence
  new.seqs <- g$seq.info[(New), Header]
  old.seqs <- g$seq.info[!(New), Header]
  retro.edges <- g$edge.info[new.seqs,old.seqs]
  min.retro.edges <- sapply(1:nrow(retro.edges), function(i){
    names(which.min(retro.edges[i,]))
  })
  
  DT <- data.table::data.table(New=new.seqs, Neighbour=min.retro.edges)
  
  return(DT)
}


##-CURRENTLY UNUSED. APE'S DIST MATRIX FUNCTION IS MORE EFFECTIVE-##
#' A wrapper for tn93's basic run function to get an edgelist 
#' NOTE: The sequences referenced here will also be referenced in another data set (seq.info)
run.tn93 <- function(seqs.full){
  #'@param seqs.full: The full alignment. Including sequences excluded from the tree.
  #'@return: An edgelist of pairwise TN93 distances calculated using tn93 binaries
  
  #Prep temporary files
  seqs.file <- tempfile("seqs.full", fileext = ".fasta")
  edgelist.file <- tempfile("edgelist", fileext = ".csv")
  ape::write.FASTA(seqs.full, seqs.file)
  
  #Obtain TN93 distances as edgelist and clean up
  system(paste0("tn93 -t 1 -o ", edgelist.file, " ", seqs.file))
  edge.info <- data.table::fread(edgelist.file)
  
  unlink(seqs.file)
  unlink(edgelist.file)
  
  return(edge.info)      
}