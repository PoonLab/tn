#'Extends an ape tree object to include annotated node and path information. 
#'By default this will include some basic node info like bootstraps. Other annotations are optional.
extend.tree <- function(t, seq.info, paths.annotated=T, mc.cores=1) {
  #'@param t: An inputted tree using ape's tree handling
  #'@param seq.info: A data frame or data.table object containing the sequences. 
  #'By default, new sequences are assigned as those not in the tree.
  #'@param mc.cores: Passed to annotate.nodes as a parallel option
  #'@param paths.annotated: An option to add path info for the tree. Useful for some tree-based clustering methods
  #'@return: the tree annotated with node information and seq.info
  
  #Root the tree (if unrooted) and resolve multichotomies
  if(!ape::is.rooted(t)){
    t <- phangorn::midpoint(t)
  }
  t <- ape::multi2di(t)
  
  #Check Sequence names inputs
  var.names <- colnames(seq.info)
  if(!("Header"%in%var.names)){
    stop("Header must be contained within var.names")
  } else{
    if(length(unique(seq.info[,"Header"]))!=length(seq.info[,"Header"])){
      seq.info <- seq.info[!duplicated(Header),]
      warning("Duplicate Headers have been arbitrarily removed")
    }
    if(!all(t$tip.label%in%seq.info$Header)){
      stop("At least 1 tip labels in tree is not represent seq.labels")
    }
  }
  seq.info <- annotate.new(seq.info, which(!(seq.info$Header%in%t$tip.label)))
  t$seq.info <- seq.info
  
  t$node.info <- annotate.nodes(t, mc.cores)
  
  if(paths.annotated){
    t$path.info <- annotate.paths(t)
  }
  
  return(t)
}


#'Called by extend.tree. Adds additional node info to the tree. Required for mono.cluster()
#'NOTE: Unlabeled nodes are defaulted to a bootstrap value of 1
annotate.nodes <- function(t, mc.cores=1) {
  #'@param t: An inputted tree using ape's tree handling. This must be annotated with seq.info
  #'@param mc.cores: A parallel option
  #'@return: A slightly more detailed data table to replace "node.info" to a given tree.
  
  #Unlabelled nodes are defaulted to a bootstrap value of 1
  t$node.label[which(t$node.label%in%"")] <- 1
  
  nodes <- 1:(2*length(t$tip.label)-1)
  
  #Store node info in data.table
  node.info <- data.table::data.table()
  node.info[,"ID" := nodes] 
  node.info[,"Bootstrap" := c(rep(100, length(t$tip.label)), as.numeric(t$node.label))]
  
  node.info[is.na(Bootstrap), 
            "Bootstrap" := 10^ceiling(log10(max(node.info$Bootstrap[!is.na(node.info$Bootstrap)])))]
  node.info[,"Bootstrap" := (node.info$Bootstrap)/max(node.info$Bootstrap)]
  
  #Get descendant information for each node
  des <- phangorn::Descendants(t, type = "all")
  node.info[, "Descendants" := des] 
  
  #Get pairwise patristic distance info
  patristic.dists <- ape::dist.nodes(t)
  dists.by.des <- lapply(des, function(x) {patristic.dists[x[x<= length(t$tip.label)],x[x<=length(t$tip.label)]]})
  node.info$max.patristic.dist <- sapply(dists.by.des, function(x){max(x)})
  node.info$mean.patristic.dist <- sapply(dists.by.des, function(x){mean(x[x>0])})
  
  return(node.info)
}


#'Called by extend.tree. Adds path info to the tree. Required for step.cluster()
annotate.paths <- function(t) {
  #'@param t: An inputted tree using ape's tree handling
  #'@return: A matrix labelled "path.info" to attach to a given tree.
  #'For each node in the path the branch lengths (below node) and bootstraps are given
  #'For the terminal node, no branch length is given below the node and the bootstrap is 1
  
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


#'Add the growth information onto the known tree.
#'This uses pplacer to ensure that previously defined clusters remain the same.
annotate.growth <- function(t, t.grown, mc.cores=1) {
  #'@param t: An inputted tree using ape's tree handling
  #'@param t.growth: A set of trees from pplacer. Each with a newly added tip.
  #'@param mc.cores: A parallel option
  #'@return: The input tree annotated with growth information stored as growth.info.

  #Obtain placement information from trees.
  #New IDs, bootstap + branch lengths of new node
  growth.info <- dplyr::bind_rows(
    parallel::mclapply(t.grown, function(x){
      
      new.ids <- setdiff(x$tip.label, t$tip.label)
      new.tips <- which(x$tip.label%in%new.ids)
      
      h <- gsub("_#.*","",new.ids[1])
      
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
      
      DT <- data.table::data.table("Header"=h, "NeighbourDes"=old.tips, "Bootstrap" = node.boots,
                                   "TermDistance"=term.dists, "PendantDistance"=pen.dists, "Terminal"=terminal)
      
      return(DT)
  }, mc.cores=mc.cores))
  
  #Collapse neighbour node descendant tips to their MRCA
  collapsed.neighbours <- sapply(growth.info[!(Terminal), NeighbourDes],  function(des){ape::getMRCA(t,des)})
  temp <- growth.info[, NeighbourDes]
  temp[which(!growth.info$Terminal)] <- collapsed.neighbours
  
  suppressWarnings(growth.info[, "NeighbourNode" := unlist(temp)])
  growth.info[, NeighbourDes := NULL]
  
  if(!all(growth.info$Header%in%t$seq.info[(New),Header])){
    warning("Not all newly added sequences are noted in the seq.info of the tree")
  }
  
  return(growth.info)
}
