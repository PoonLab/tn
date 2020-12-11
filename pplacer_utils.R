require(ape)
require(rjson)
require(digest)

#Split new cases from old cases based on a Marker
newSplit <- function(iFile, newMark) {
  #@param iFile: Input file (a full set of fasta sequences)
  #@param newMark: The marker defining new sequences from a set of time-defined seqs
  #                This must be something that exists in the headers of ONLY new sequences
  
  #Take fasta file input and identify new sequences
  seqs <- read.FASTA(iFile)
  iNewSeq <- grepl(newMark, names(seqs))
  
  #Two files are made to delimit old and new sequences
  oFileO <- paste0(gsub(".fasta$", "_Old.fasta", iFile))
  oFileN <- paste0(gsub(".fasta$", "_New.fasta", iFile))
  write.FASTA(seqs[!iNewSeq], oFileO)
  write.FASTA(seqs[iNewSeq], oFileN)
}

#Obtain n random samples from a given fasta file, save this to a directory. 
#Each file in the directory will simply be labelled as "samp<i>.fasta"
#Alternatively - If a marker for "New" sequences is added, this can be used 
sampleFasta <- function(iFile, sampsDir, n=100, prop=0.8, short=NA, nFile=NA) {
  #@param iFile: Input file (For robustness tests, this is a filtered set of fasta sequences)
  #@param sampsDir: The filepath to the output directory for samples
  #@param n: The number of random samples to take.
  #@param nFile: The file containing specifically the newest sequences
  #@param prop: The proportion of the sequences to be taken as a random sample
  #@param short: The short form label to be used for these sequences.
  #              This is defaulted to the name of the fasta file if not supplied
  #              Any "_" characters in short names will be replaced with "-"
  
  #Take file input
  seqs <- read.FASTA(iFile)
  
  #Standardize sampsDir to include terminal / and create dir if it doesn't exist
  sampsDir <- gsub("$|/$", "/", sampsDir)
  if(!file.exists(sampsDir)) {dir.create(sampsDir)}
  
  #Default nickname is that of the fasta file
  if(is.na(short)){
    temp <- strsplit(iFile, "/|\\.")[[1]]
    short <- temp[length(temp)-1]
    short <- gsub("_", "-", short)
  }
  
  #Sample and write to sample file 
  for(i in 1:n) {
    samp <- sample(seqs, prop*(length(seqs)), replace = F)
    oFile <- paste0(sampsDir, short, "_samp", i, ".fasta")
    write.FASTA(samp, oFile)
  }
  
  #Attach new sequences to the set of samples to obtain sets of "complete" fastas
  if(is.na(nFile)) {nFile <- gsub("_Old", "_New", iFile)}
  if(file.exists(nFile)) {
    
    seqsN <- read.FASTA(oFileN)
    seqFiles <- list.files(sampsDir)
    
    #Loop through samples
    for(i in 1:length(seqFiles)) {
      iSeqs <- read.FASTA(paste0(sampsDir, seqFiles[i]))
      oSeqs <- c(iSeqs, seqsN)
      oFile <- paste0(sampsDir, short, "_full", i, ".fasta")
      write.FASTA(oSeqs, oFile)
    }
  }
}

#Run FastTree on a set of sequences generated through sample fasta
#The output is also saved to the samps directory
multiTree <- function(sampsDir) {
  #@param sampsDir: The filepath to the output directory for samples
  
  #Pull sequence files from sample directory
  sampsDir <- gsub("$|/$", "/", sampsDir)
  seqFiles <- list.files(sampsDir, pattern = "\\_samp[0-9]+.fasta$")
  n <- length(seqFiles)
  
  #Standardize sampsDir to include terminal / and standardize nickname
  short <- (strsplit(seqFiles[1], "_")[[1]])[1]
  
  #Create FastTree for each. Logfiles are recorded for future extraction of stats.
  #This uses GTR methods, and creates a Newick format tree
  for(i in 1:n) {
    system(paste0("FastTree -gtr -log ", sampsDir, short, "_log", i, ".txt ",
                  "-nt ", sampsDir, short, "_samp", i, ".fasta ",
                  "> ", sampsDir, short, "_tree", i, ".nwk"))
  }
}

##POINT OF REVIEW##

#A wrapper for the python translator script
#This translates the FastTree output into a .json "stats" file 
multiTranslator <- function(sampsDir, fPath) {
  #@param sampsDir: The filepath to the output directory for samples
  #@param fPath: Path to the python translator function

  #Standardize sampsDir to include terminal / and standardize nickname
  sampsDir <- gsub("$|/$", "/", sampsDir)
  short <- (strsplit(logFiles[1], "_")[[1]])[1]
    
  #Pull log files from sample directory
  logFiles <- list.files(sampsDir, pattern = "\\.txt$")
  n <- length(logFiles)
  
  #For each log file, create a stats 
  for(i in 1:n) {
    system(paste0("python3", fPath, "-i ", sampsDir, short, "_log",
                  " -o ", sampsDir, short, "_stats", i, ".json"))
  }
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