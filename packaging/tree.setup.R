#Translates a set of sequence headers into a data.frame object for inputted.
#NOTE: This must contain, at minimum, a set of unique sequence id's (labelled ID) and one other variable

#'@param seqs: An inputted alignment using ape's sequence handling
#'@param var.names: The names of the variables represented in each header. This must contain "ID".
#'@param var.transformations: A list of transformation functions (such as as.character()) 
#'these transform each row into it's proper type. by default, each type is set to character.
#'The variable named "ID" will be automatically forced to a character object.
#'@param sep: The separator character that splits upthe headers in the fasta file
#'@return: A data.table object containing the information associated with each sequence
pull.headers <- function(seqs, var.names, var.transformations=list(), sep="_") {
  
  #Checking Inputs
  if(length(var.names)!=length(unique(var.names))){
    stop("var.names may not contain repeats")
  }
  if(length(var.transformations)==0){
    var.transformations <- lapply(1:length(var.names), function(x){as.character})
  } else {
    if(length(var.names)!=length(var.transformations)) {
      stop("var.names and var.transformations must be equal lengths")
      return(NULL)
    }
  }
  if(!("ID"%in%var.names)){
    stop("'ID' must be contained within var.names")
  } else{
    var.transformations[[which(var.names%in%"ID")]] <- as.character
  }
  
  #Split and transform data from headers
  split.headers <- sapply(names(seqs), function(x) {strsplit(x, sep)[[1]]})
  seq.info <- lapply(1:nrow(split.headers), function(i) {
    x <- unname(split.headers[i,])
    x <- var.transformations[[i]](x)
    data.table::data.table(x)
  })
  
  seq.info <- dplyr::bind_cols(seq.info)
  colnames(seq.info) <- var.names
  
  return(seq.info)
}


#Extends an ape tree object to include annotated node and path information. 
#By default this will include some basic node info like bootstraps. Other annotations are optional.

#'@param t: An inputted tree using ape's tree handling
#'@param seq.info: A data frame or data.table object containing the sequences
#'@param mc.cores: Passed to annotate.nodes as a parallel option
#'@param paths.annotated: An option to add path info for the tree. Useful for some tree-based clustering methods
#'@return: the tree annotated with node information and seq.info
extend.tree <- function(t, seq.info, paths.annotated=T, mc.cores=1) {
  
  #Root the tree (if unrooted) and resolve multichotomies
  if(!ape::is.rooted(t)){
    t <- phangorn::midpoint(t)
  }
  t <- ape::multi2di(t)
  
  #Check Sequence names inputs
  var.names <- colnames(seq.info)
  if(!("ID"%in%var.names)){
    stop("'ID' must be contained within var.names")
  } else{
    seq.info[,ID:=as.character(ID)]
    if(length(unique(seq.info[,"ID"]))!=length(seq.info[,"ID"])){
      seq.info <- seq.info[!duplicated(ID),]
      warning("Duplicate ID's have been arbitrarily removed")
    }
  }
  t$seq.info <- seq.info
  
  print("Annotating Nodes")
  t$node.info <- annotate.nodes(t, mc.cores)
  
  print("Annotating Paths")
  if(paths.annotated){
    t$path.info <- annotate.paths(t)
  }
  
  return(t)
}


#Called by extend.tree. Adds additional node info to the tree. Required for mono.cluster()
#NOTE: Unlabelled nodes are defaulted to a bootstrap value of 1

#'@param t: An inputted tree using ape's tree handling. This must be annotated with seq.info
#'@param mc.cores: A parallel option
#'@return: A slightly more detailed data table to replace "node.info" to a given tree.
annotate.nodes <- function(t, mc.cores=1) {
  
  #Unlabelled nodes are defaulted to a bootstrap value of 1
  t$node.label[which(t$node.label%in%"")] <- 1
  
  nodes <- 1:(2*length(t$tip.label)-1)
  
  #Store node info in data.table
  node.info <- data.table()
  node.info[,"ID" := c(t$seq.info$ID, nodes[(nrow(t$seq.info)+1):length(nodes)])] 
  node.info[,"Bootstrap" := c(rep(100, nrow(t$seq.info)), as.numeric(t$node.label))]
  
  node.info[is.na(Bootstrap), 
            "Bootstrap" := 10^ceiling(log10(max(node.info$Bootstrap[!is.na(node.info$Bootstrap)])))]
  node.info[,"Bootstrap" := (node.info$Bootstrap)/max(node.info$Bootstrap)]
  
  #Get descendant information for each node
  des <- phangorn::Descendants(t, type = "all")
  node.info$des <- des
  node.info$size <- sapply(des, function(x){length(x)})
  
  #Get pairwise patristic distance info
  patristic.dists <- ape::dist.nodes(t)
  dists.by.des <- lapply(des, function(x) {patristic.dists[x[x<=nrow(t$seq.info)],x[x<=nrow(t$seq.info)]]})
  node.info$max.patristic.dist <- sapply(dists.by.des, function(x){max(x)})
  node.info$mean.patristic.dist <- sapply(dists.by.des, function(x){mean(x[x>0])})
  
  #Get membership 
  for(nm in colnames(t$seq.info)[-1]){
    node.info[,(nm) := mclapply(des, function(x){
      t$seq.info[x[x<=nrow(t$seq.info)],get(nm)]
    }, mc.cores=mc.cores)]
  }
  node.info[, "Membership" := mclapply(des, function(x){
    t$seq.info[x[x<=nrow(t$seq.info)],(ID)]
  }, mc.cores=mc.cores)]
  
  return(node.info)
}


#Called by extend.tree. Adds path info to the tree. Required for step.cluster()

#'@param t: An inputted tree using ape's tree handling
#'@return: A matrix labelled "path.info" to attach to a given tree.
#'For each node in the path the branch lengths (below node) and bootstraps are given
#'For the terminal node, no branch length is given below the node and the bootstrap is 1
annotate.paths <- function(t) {
  
  #Get paths and length information from terminal nodes
  lens <- ape::node.depth.edgelength(t)
  paths <- ape::nodepath(t)
  names(paths) <- sapply(paths, function(p){p[length(p)]})
  
  #Extend apes nodepath() function to internal nodes
  i <- 1
  while(length(paths) < nrow(t$node.info)+1) {
    new.paths <- sapply(paths[i:length(paths)], function(x){x[-length(x)]})
    new.paths <- new.paths[!duplicated(new.paths)]  
    names(new.paths) <- sapply(new.paths, function(p){p[length(p)]})
    
    i <- length(paths)+1
    paths <- c(paths, new.paths[which(!(names(new.paths)%in%names(paths)))])
  }
  paths <- paths[which(!sapply(paths, function(p) {length(p)})==0)]
  paths <- paths[order(as.numeric(names(paths)))]
  lens <- lapply(paths, function(x){c((lens[x[-1]]-lens[x[-length(x)]]),NA) })
  
  #Obtain bootstrap branchlength and node number information for all paths
  boots <- lapply(paths, function(x){t$node.info$Bootstrap[x]})
  path.info <- lapply(1:nrow(t$node.info), function(i) {
    m <- matrix(ncol=length(paths[[i]]), nrow=3)
    rownames(m) <- c("Node", "Boot", "BranchLength")
    m[1,] <- rev(paths[[i]])
    m[2,] <- rev(boots[[i]])
    m[3,] <- rev(lens[[i]])
    return(m)
  })
  
  return(path.info)
}


#Add the growth information onto the known tree.
#This uses pplacer to ensure that previously defined clusters remain the same.

#'@param t: An inputted tree using ape's tree handling
#'@param t.growth: A set of trees from pplacer. Each with a newly added tip.
#'@param mc.cores: A parallel option
#'@return: The input tree annotated with growth information stored as growth.info.
annotate.growth <- function(t, t.grown, mc.cores=1) {
  
  print("Annotating Growth")
  
  #Obtain placement information from trees.
  #New IDs, bootstap + branch lengths of new node
  growth.info <- dplyr::bind_rows(
    parallel::mclapply(1:length(t.grown), function(i){
      
      x <- t.grown[[i]]
      new.ids <- setdiff(x$tip.label, t$tip.label)
      new.tips <- which(x$tip.label%in%new.ids)
      
      new.edges <- which(x$edge[,2]%in%new.tips)
      term.dists <- x$edge.length[new.edges]
      
      new.nodes <- x$edge[new.edges, 1]
      new.nodes.des <- which(x$edge[,1]%in%new.nodes)
      node.boots <- sapply(new.ids, function(x){
        l <- strsplit(x,"=")[[1]]
        as.numeric(l[length(l)])
      })
      
      pen.edges <- dplyr::setdiff(new.nodes.des, new.edges)
      pen.dists <- x$edge.length[pen.edges]
      
      neighbour.des <- lapply(new.nodes, function(n){phangorn::Descendants(x, n, "tips")[[1]]})
      old.tips <- lapply(neighbour.des, function(des){which(t$tip.label%in%x$tip.label[des])})
      terminal <- sapply(old.tips, function(x){length(x)==1})
      
      DT <- data.table::data.table("ID"=i, "NeighbourDes"=old.tips, "Bootstrap" = node.boots,
                                   "TermDistance"=term.dists, "PendantDistance"=pen.dists, "Terminal"=terminal)
      
      return(DT)
  }, mc.cores=mc.cores))
  
  #Collapse neighbour node descendant tips to their MRCA
  collapsed.neighbours <- sapply(growth.info[!(Terminal), NeighbourDes],  function(des){ape::getMRCA(t,des)})
  temp <- growth.info[, NeighbourDes]
  temp[which(!growth.info$Terminal)] <- collapsed.neighbours
  
  suppressWarnings(growth.info[, "NeighbourNode" := unlist(temp)])
  growth.info[, NeighbourDes := NULL]
  
  return(growth.info)
}
