#' Parses the logfile output of tree building software to find information relevant to pplacer.
#' Prints stats to a temporary ref.pkg file and returns the path to that file. 
#' 
#' NOTE: Currently compatible with GTR substitution model and either FastTree, RAxML or IQ-TREE logfiles.
translate.log <- function(log.file, program, substitution.model="GTR") {
  #'@param log.file: A path to the logfile from a tree construction run
  #'@param program: The software used to build the tree.
  #'Can be "FastTree", "IQ-TREE" or "RAxML".
  #'@param substitution.model: The substitution model used. Currently only accepts "GTR"
  #'@return: a json output to be written to a given stats.file for pplacer
  
  #Open connection to log.file
  con <- file(log.file)
  lns <- readLines(con)
  close(con)
  
  #Extracts and normalizes list of frequencies
  if(program%in%"FastTree") {
    p <- "GTRRates"
    s <- lns[grep(p, lns)]
    s <- strsplit(s, '\t')[[1]]
    s <- as.numeric(s[c(2,3,4,5,6,7)])
  }
  if(program%in%"IQ-TREE") {
    p <- "Rate parameters"
    s <- lns[grep(p, lns)]
    s <- strsplit(s[length(s)], " ")[[1]]
    s <- as.numeric(s[c(5,20,11,8,14,17)])
  }
  
  #Write stats information to .json
  stats.json <- jsonlite::toJSON(list(
    "empirical_frequencies"=TRUE,
    "datatype"="DNA",
    "subs_model"=substitution.model,
    "program"=program,
    ##-TO-DO: Test Correctness of gamma assumption -##
    "ras_model"="gamma",
    "gamma"=list(
      "alpha"=1.0,
      "n_cats"=as.integer(20)
    ),
    "subs_rates"=list(
      "ac"=s[1],
      "gt"=s[2],
      "at"=s[3],
      "ag"=s[4],
      "cg"=s[5],
      "ct"=s[6]
    )
  ), pretty=T, always_decimal=T, auto_unbox=T)

  return(stats.json)
}


#' A wrapper for the taxit create function used by pplacer. This will generate a summary json
#' See pplacer's basic function regarding alignments and tree function
#' NOTE: Creates temporary files. These are only deleted with the end of the session 
taxit.create <- function(t, seqs.full, stats.json, locus="LOCUS") {
  #'@param t: The tree (made on a subset of the full alignment)
  #'@param seqs.full: The full alignment. Including sequences excluded from the tree.
  #'@param stats.json: Path to the full alignment file (new + old seqs)
  #'@param locus: Extra information required for the summary json
  #'@return: Path to a temporary refpkg directory
  
  #Set up and populate temporary file system
  temp.dir <- tempdir()
  
  ali.file <- tempfile("ali", temp.dir,".fasta")
  ali.file.name <- tail(strsplit(ali.file, "/")[[1]],1)
  ape::write.FASTA(seqs.full, ali.file)
  
  stats.file <- tempfile("stats", temp.dir, ".json")
  stats.file.name <- tail(strsplit(stats.file, "/")[[1]],1)
  jsonlite::write_json(stats.json, stats.file)
  
  tree.file <- tempfile("tree", temp.dir,".nwk")
  tree.file.name <- tail(strsplit(tree.file, "/")[[1]],1)
  ape::write.tree(t, tree.file)
  
  log.file <- tempfile("log", temp.dir, ".txt")
  log.file.name <- tail(strsplit(log.file, "/")[[1]],1)
  write(log.file, "Sample log file. Created Using taxit.create() wrapper")
  
  #Create generic JSON summary for refpkg.
  summary.json <- jsonlite::toJSON(list(
    "files" = list(
      "aln_fasta" = ali.file.name,
      "phylo_model" = stats.file.name,
      "tree" = tree.file.name,
      "tree_stats" = log.file.name
    ),
    "rollback"=NULL,
    "log"= c("Stripped refpkg (removed 0 files)",
             "Loaded initial files into empty refpkg"),
    "metadata" = list(
      "create_date" = as.character(Sys.Date()),
      "format_version" = "1.1",
      "locus" = locus 
    ),
    "rollforward"=NULL,
    "md5"= list(
      "aln_fasta" = digest::digest(ali.file.name, algo = "md5"),
      "phylo_model" = digest::digest(stats.file.name, algo = "md5"),
      "tree" = digest::digest(tree.file.name, algo = "md5"),
      "tree_stats" = digest::digest(log.file.name, algo = "md5")
    )
  ), pretty=T, always_decimal=T, auto_unbox=T)
  summary.file <- tempfile("summary", temp.dir,".json")
  jsonlite::write_json(summary.json, summary.file)
  
  return(temp.dir)
}

#' A wrapper for pplacer's basic run function coupled with guppy's sing function. 
#' Together, these extend fixed trees with most likely placement locations.
#' TODO: Include package binaries such that it is not a requirement to install both pplacer and guppy
run.pplacer_guppy <- function(refpkg){
  #'@param: A reference package to use as input for pplacer
  #'@return: A set of trees

}
