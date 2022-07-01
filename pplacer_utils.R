 require("ape")
require("rjson")
require("digest")
require("jsonlite")

##TO-DO: Too Bulky - use multi-tree files or temp files or compression to save space -##

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
  t <- read.tree(paste0(sampsDir, short, "_tree", i, ".nwk"))
  t <- multi2di(t)
  write.tree(t, paste0(sampsDir, short, "_tree", i, ".nwk"))
}

#A wrapper for the python translator script
#This translates the FastTree output into a .json "stats" file 
translator <- function(logF, program, oFile, dataType="DNA", subsModel="GTR") {
  #@param sampsDir: The filepath to the output directory for samples
  #@param program: Changeshow log files are parsed depending on program used
  #@param oFile: The output file for this stats jason
  #@param dataType: The type of data of interest. Can alternatively be AA
  #                 CURRENTLY ONLY BUILT TO HANDLE DNA
  #@param subsModel: The substitution model used
  #                  CURRENTLY ONLY BUILT TO HANDLE GTR
  
  #Open connection to 
  con <- file(logF)
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
    subsModel
    s <- lns[grep(p, lns)]
    s <- strsplit(s[length(s)], " ")[[1]]
    s <- as.numeric(s[c(5,20,11,8,14,17)])
  }
  
  #Write stats information to .json
  con <- file(oFile)
  j <- toJSON(list(
    "empirical_frequencies"=TRUE,
    "datatype"=dataType,
    "subs_model"=subsModel,
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
  write(j, con)
  close(con)
}

#A wrapper for the taxit create function found in pplacer
#Can be run in Isolation
taxitCreate <- function(treeF, logF, fullF, oDir, program="FastTree", locus="LOCUS") {
  #@param treeF: Path to the newick tree file
  #@param program: Passed to translator
  #@param fullF: Path to the full alignment file (new + old seqs)
  #@param logF: Path to the log file from the tree building process
  #@param locus: Extra information required for the summary json
  #@param oDir: The output refpackage directory name
  
  #Standardize oDir to include terminal / and create dir if it doesn't exist
  oDir <- gsub("$|/$", "/", oDir)
  if(!file.exists(oDir)) {dir.create(oDir)}
  
  #Copy files and update filenames
  file.copy(treeF, oDir)
  file.copy(logF, oDir)
  file.copy(fullF, oDir)
    
  #Trim complete filenames
  logFT <- tail(strsplit(logF, "/")[[1]],1)
  treeFT <- tail(strsplit(treeF, "/")[[1]],1)
  fullFT <- tail(strsplit(fullF, "/")[[1]],1)
  
  #Create stats file
  statsFT <- paste0(gsub(".log|_log.txt$", "_stats.json", logFT))
  translator(logF, program, paste0(oDir, statsFT))
  
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
  ), pretty=T, always_decimal=T, auto_unbox=T)
  write(conText, conF)
  close(conF)
}

#A function to call taxit create many times on a directory of repeated runs
#This is expected to sort through hundreds of files in a single folder, sorting sets into refpackages
multiTaxit <- function(sampsDir) {
  #@param sampsDir: The filepath to the output directory for samples
  #                 This should have standardized names
  
  #Standardize sampsDir and oDir to include terminal / and create dir if it doesn't exist
  sampsDir <- gsub("$|/$", "/", sampsDir)
  
  #Obtain nickname and sample number from sample directory
  seqFiles <- list.files(sampsDir, pattern = "\\_samp[0-9]+.fasta$")
  n <- length(seqFiles)
  short <- (strsplit(seqFiles[1], "_")[[1]])[1]
  
  #Initialize output directory
  oDir <- paste0(sampsDir, short, "_refpackages")
  oDir <- gsub("$|/$", "/", oDir)
  if(!file.exists(oDir)) {dir.create(oDir)}
  
  #Iterate through listed files and create refpackages for each.
  for (i in 1:n) {
    
    treeF <- paste0(sampsDir, short, "_tree", i, ".nwk")
    logF <- paste0(sampsDir, short, "_log", i, ".txt")
    fullF <- paste0(sampsDir, short, "_full", i, ".fasta")
    refDir <- paste0(oDir, short, "_refpkg", i)
    
    taxitCreate(treeF, logF, fullF, refDir)
    
    #Tidy up original files
    file.remove(treeF)
    file.remove(logF)
    file.remove(fullF)
  }
}

# Simple wrapper for running pplacer on a refpackage
pplacer_guppy <- function(refpkg, oDir) {
  refpkgnm <- (strsplit(refpkg, '/')[[1]])[-1]
  short <- (strsplit(refpkgnm, "[.]")[[1]])[1]
  jpFile <- paste0(oDir, short, ".jplace")

  seqFile <- list.files(refpkg)[grepl(".fasta", list.files(refpkg))]
  system(paste0("pplacer -c ", refpkg, " -o ", jpFile, " ", refpkg, "/", seqFile))
  system(paste0("guppy sing ", jpFile, " -o ", oDir, short, "_growth.tree"))
}

#Runs pplacer on all files in a directory to obtain placement files
multiPplacer <- function(sampsDir) {
  #@param sampsDir: The filepath to the output directory for samples
  #                 This should have standardized names
  
  #Standardize sampsDir and oDir to include terminal / and create dir if it doesn't exist
  sampsDir <- gsub("$|/$", "/", sampsDir)
  
  #Obtain nickname and sample number from sample directory
  refDir <- list.files(sampsDir, pattern = "\\_refpackages$")
  short <- (strsplit(refDir, "_")[[1]])[1]
  
  #Initialize output directory
  oDir <- paste0(sampsDir, short, "_jplaceFiles")
  oDir <- gsub("$|/$", "/", oDir)
  if(!file.exists(oDir)) {dir.create(oDir)}

  n <- length(list.files(paste0(sampsDir, refDir)))
  for (i in 1:n) {
    refPkg <- paste0(sampsDir, refDir, "/", short, "_refpkg", i)
    system(paste0("pplacer -c ", refPkg, " -o ", oDir, short, i, ".jplace", " ", 
                  refPkg, "/", short, "_full", i,  ".fasta"))
    
  }
}

#Runs guppy on all files in a directory to obtain growth files (inputs for tree growth method)
multiGuppy <- function(sampsDir) {
  #@param sampsDir: The filepath to the output directory for samples
  #                 This should have standardized names
  
  #Standardize sampsDir and oDir to include terminal / and create dir if it doesn't exist
  sampsDir <- gsub("$|/$", "/", sampsDir)
  
  #Obtain nickname and sample number from sample directory
  jplaceDir <- list.files(sampsDir, pattern = "\\_jplaceFiles$")
  short <- (strsplit(jplaceDir, "_")[[1]])[1]
  
  #Initialize output directory
  oDir <- paste0(sampsDir, short, "_GrowthFiles")
  oDir <- gsub("$|/$", "/", oDir)
  if(!file.exists(oDir)) {dir.create(oDir)}
  
  n <- length(list.files(paste0(sampsDir, jplaceDir)))
  for (i in 1:n) {
    jpFile <- paste0(sampsDir, jplaceDir, "/", short, i, ".jplace")
    system(paste0("guppy sing ", jpFile, " -o ", oDir, short, "_growth", i, ".tree"))
  }
}