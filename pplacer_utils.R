require(ape)
require(rjson)
require(digest)

#Obtain n samples from a given fasta file, save this to a directory
sampleFasta <- function(iFile, sampsDir, n=100){
  
  seqs <- read.FASTA(iFile)
  
  for(i in 1:n) {
    samp <- sample(seqs, 0.8*(length(seqs)), replace = F)
    s <- paste0(sampsDir, "/samp", i, ".fasta")
    write.FASTA(samp, s)
  }
}

#Run FastTree on a set of sequences generated through sample fasta
multiRun <- function(sampsDir) {
  
  n <- length(list.files(sampsDir, pattern = "\\.fasta$"))
  
  for(i in 1:n) {
    print(i)
    system(paste0("FastTree -gtr -log ", sampsDir, "/log", i, ".txt ",
                  "-nt ", sampsDir, "/samp", i, ".fasta ",
                  "> ", sampsDir, "/tree", i, ".nwk"))
  }
}

#A wrapper for the python translator script
#This translates the FastTree output into a .json file which can be used in .jplace file for 
translator <- function(sampsDir) {
  
  files <- list.files(sampsDir, patter="\\.txt$")
  
  for(i in 1:length(files)) {
    system(paste0("python3 translator.py -i ", sampsDir, "/", files[i],
                  " -o ", sampsDir, "/ststats", i, ".json"))
  }
}

#Effectively the reverse of makeNSFile
makeFiltFile <-function(iFile, new, oFile) {
  seqs <- read.FASTA(iFile)
  seqs <- seqs[!grepl(new, names(seqs))]
  
  write.FASTA(seqs, oFile)
}

#Find new sequences in a complete input file and separate them in a fast
#New sequences are based on a series of strings in the 
getNS <- function(iFile, oFile, new){
  
  seqs <- read.FASTA(iFile)
  seqs <- seqs[grepl(new, names(seqs))]
  write.FASTA(seqs, oFile)
}

#A similar function to the taxit create function found in pplacer
taxitCreate <- function(treeF, logF, statsF, locus="pol_SubB", refpkg) {
  
  #Copy files and update filenames to trim complete paths
  file.copy(logF, refpkg)
  file.copy(treeF, refpkg)
  file.copy(statsF, refpkg)
  logF <- tail(strsplit(logF, "/")[[1]],1)
  treeF <- tail(strsplit(treeF, "/")[[1]],1)
  statsF <- tail(strsplit(statsF, "/")[[1]],1)
  
  #Create contents JSON file
  conF <- file(paste0(refpkg, "/CONTENTS.json"))
  conText <- toJSON(list(
    "files" = list(
      "aln_fasta" = "fullAln.fasta" ,
      "phylo_model" = statsF,
      "tree" = treeF,
      "tree_stats" = logF
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
      "aln_fasta" = digest("fullAln.fasta", algo = "md5"),
      "phylo_model" = digest(logF, algo = "md5"),
      "tree" = digest(treeF, algo = "md5"),
      "tree_stats" = digest(logF, algo = "md5")
    )
  ), indent=4)
  write(conText, conF)
  close(conF)
}

#A function to call taxit create many times on a directory of repeated runs
multiTaxit <- function(sampsDir, n, short="", nSeqs, oDir) {
  
  oDir <- paste0(oDir, "/", short, "refpackages")
  dir.create(oDir)
  
  #
  seqs <- lapply(1:n, function(i) {
    samp <- paste0(sampsDir, "/", short, "samp", i, ".fasta")
    seqs <- c(read.FASTA(samp), read.FASTA(nSeqs))
    return(seqs)
  })
  
  for (i in 1:n) {
    #File Manipulation and directory setup 
    refpkg <- paste0(oDir, "/", short, i, ".refpkg")
    dir.create(refpkg)
    
    write.FASTA(seqs[[i]], paste0(refpkg, "/fullAlign.fasta"))
    
    treeF <- paste0(sampsDir, "/", short, "tree", i, ".nwk")
    logF <- paste0(sampsDir, "/", short, "log", i, ".txt")
    statsF <- paste0(sampsDir, "/", "st", "stats", i, ".json")
    
    taxitCreate(treeF=treeF, logF=logF, statsF=statsF, refpkg = refpkg)
  }
}

#Runs pplacer on all files in a directory to obtain placement files
multiPplacer <- function(refDir, oDir, short, n=100) {
  for (i in 1:n) {
    f <- paste0(refDir, "/", short, i, ".refpkg")
    system(paste0("pplacer -c ", f, " -o ", oDir, "/jplaceFile.jplace", i, " ", f, "/fullAlign.fasta"))
    
  }
}

#Runs guppy on all files in a directory to obtain growth files (inputs for tree growth method)
multiGuppy <- function(refDir, oDir, n=100) {
  for (i in 1:n) {
    f <- paste0(refDir, "/jplaceFile.jplace", i)
    system(paste0("guppy sing ", f, " -o ", oDir, "/growthFile", i, ".tree"))
  }
}