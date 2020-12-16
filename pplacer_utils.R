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
    
    seqsN <- read.FASTA(nFile)
    
    #Loop through samples
    for(i in 1:n) {
      iSeqs <- read.FASTA(paste0(sampsDir, short, "_samp", i, ".fasta"))
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

#A wrapper for the python translator script
#This translates the FastTree output into a .json "stats" file 
multiTranslator <- function(sampsDir, fPath) {
  #@param sampsDir: The filepath to the output directory for samples
  #@param fPath: Path to the python translator function

  #Standardize sampsDir to include terminal / and standardize nickname
  #In addition, pull log files from sample directory
  sampsDir <- gsub("$|/$", "/", sampsDir)
  logFiles <- list.files(sampsDir, pattern = "\\_log[0-9]+.txt$")
  short <- (strsplit(logFiles[1], "_")[[1]])[1]
  n <- length(logFiles)
  
  #For each log file, create a stats json file
  for(i in 1:n) {
    system(paste0("python3 ", fPath, " -i ", sampsDir, short, "_log", i, ".txt",
                  " -o ", sampsDir, short, "_stats", i, ".json"))
  }
}

#A wrapper for the taxit create function found in pplacer
#Can be run in Isolation
taxitCreate <- function(treeF, logF, statsF, fullF, oDir, locus="LOCUS") {
  #@param treeF: Path to the newick tree file
  #@param fullF: Path to the full alignment file (new + old seqs)
  #@param logF: Path to the log file from the tree building process
  #@param statsF: Path to the stats .json file created by translator
  #@param locus: Extra information required for the summary json
  #@param oDir: The output refpackage directory name
  
  #Standardize oDir to include terminal / and create dir if it doesn't exist
  oDir <- gsub("$|/$", "/", oDir)
  if(!file.exists(oDir)) {dir.create(oDir)}
  
  #Copy files and update filenames
  file.copy(treeF, oDir)  
  file.copy(logF, oDir)
  file.copy(statsF, oDir)
  file.copy(fullF, oDir)
  
  #Trim complete filenames
  logFT <- tail(strsplit(logF, "/")[[1]],1)
  treeFT <- tail(strsplit(treeF, "/")[[1]],1)
  statsFT <- tail(strsplit(statsF, "/")[[1]],1)
  fullFT <- tail(strsplit(fullF, "/")[[1]],1)
  
  #Create contents JSON file
  conF <- file(paste0(oDir, "/CONTENTS.json"))
  conText <- toJSON(list(
    "files" = list(
      "aln_fasta" = fullFT,
      "phylo_model" = statsFT,
      "tree" = treeFT,
      "tree_stats" = logFT
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
      "aln_fasta" = digest(fullFT, algo = "md5"),
      "phylo_model" = digest(logFT, algo = "md5"),
      "tree" = digest(treeFT, algo = "md5"),
      "tree_stats" = digest(statsFT, algo = "md5")
    )
  ), indent=4)
  write(conText, conF)
  close(conF)
}

#A function to call taxit create many times on a directory of repeated runs
#This is expected to sort through hundreds of files in a single folder, sorting sets into refpackages
multiTaxit <- function(sampsDir, oDir=NA, locus="LOCUS") {
  #@param sampsDir: The filepath to the output directory for samples
  #@param oDir: The output directory for the sample directory
  #@param locus: Passed to taxit create for the Contents json. 
  
  #Standardize sampsDir and oDir to include terminal / and create dir if it doesn't exist
  sampsDir <- gsub("$|/$", "/", sampsDir)
  
  #Obtain nickname and sample number from sample directory
  seqFiles <- list.files(sampsDir, pattern = "\\_samp[0-9]+.fasta$")
  n <- length(seqFiles)
  short <- (strsplit(seqFiles[1], "_")[[1]])[1]
  
  #Initialize output directory
  if(is.na(oDir)) {oDir <- paste0(sampsDir, short, "_refpackages")}
  oDir <- gsub("$|/$", "/", oDir)
  if(!file.exists(oDir)) {dir.create(oDir)}
  
  #Iterate through listed files and create refpackages for each.
  for (i in 1:n) {
    
    treeF <- paste0(sampsDir, short, "_tree", i, ".nwk")
    logF <- paste0(sampsDir, short, "_log", i, ".txt")
    statsF <- paste0(sampsDir, short, "_stats", i, ".json")
    fullF <- paste0(sampsDir, short, "_full", i, ".fasta")
    
    refDir <- paste0(oDir, short, "_refpkg", i)
    
    taxitCreate(treeF, logF, statsF, fullF, refDir, locus, noCopies = T)
    
    #Tidy up original files
    file.remove(treeF)
    file.remove(logF)
    file.remove(statsF)
    file.remove(fullF)
  }
}

#Runs pplacer on all files in a directory to obtain placement files
multiPplacer <- function(iDir, oDir) {
  #@param iDir: Input directory of reference packages
  #@param oDir: Output directory of .jplace files
  
  #Standardize iDir and oDir to include terminal / and create dir if it doesn't exist
  iDir <- gsub("$|/$", "/", iDir)
  oDir <- gsub("$|/$", "/", oDir)
  if(!file.exists(oDir)) {dir.create(oDir)}
  
  #Obtain nickname and sample number from sample directory
  refDirs <- list.files(iDir, pattern = "\\_refpkg[0-9]+$")
  n <- length(refDirs)
  short <- (strsplit(refDirs[1], "_")[[1]])[1]
  
  for (i in 1:n) {
    refDir <- paste0(iDir, short, "_refpkg", i)
    system(paste0("pplacer -c ", refDir, " -o ", oDir, short, i, ".jplace", " ", 
                  refDir, "/", short, "_full", i,  ".fasta"))
    
  }
}

#Runs guppy on all files in a directory to obtain growth files (inputs for tree growth method)
multiGuppy <- function(iDir, oDir) {
  #@param iDir: Input directory of .jplace files
  #@param oDir: Output directory of growth files
  
  #Standardize iDir and oDir to include terminal / and create dir if it doesn't exist
  iDir <- gsub("$|/$", "/", iDir)
  oDir <- gsub("$|/$", "/", oDir)
  if(!file.exists(oDir)) {dir.create(oDir)}
  
  #Obtain nickname and sample number from sample directory
  jpFiles <- list.files(iDir, pattern = "\\.jplace$")
  n <- length(jpFiles)
  short <- (strsplit(jpFiles[1], "_")[[1]])[1]
  
  for (i in 1:n) {
    jpFile <- paste0(iDir, short, i, ".jplace")
    system(paste0("guppy sing ", jpFile, " -o ", oDir, short, "_growth", i, ".tree"))
  }
}