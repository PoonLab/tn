require(ape)
require(rjson)

#Obtain n samples from a given fasta file, save this to a directory
sampleFasta <- function(iFile, n=100){
  
  seqs <- read.FASTA(iFile)
  
  for(i in 1:n) {
    samp <- sample(seqs, 0.8*(length(seqs)), replace = F)
    dir.create(paste0(iFile,"Samps"))
    setwd(paste0(iFile,"Samps"))
    s <- paste0("samp", i, ".fasta")
    write.FASTA(samp, s)
  }
}

#Run FastTree on a set of sequences generated through sample fasta
multiRun <- function(sampsDir, n) {
  
  for(i in 1:n) {
    system(paste0("FastTree -gtr -log ", sampsDir, "/log", i, ".txt ",
           "-nt ", sampsDir, "/samp", i, ".fasta ",
           "> ", sampsDir, "/tree", i, ".nwk"))
  }
}

#A wrapper for the python translator script
#This translates the FastTree output into a .json file which can be used in .jplace file for 
translator <- function(sampsDir, n) {
  
  files <- list.files(sampsDir, patter="\\.txt$")
  
  for(i in 1:length(files)) {
    system(paste0("python3 translator.py -i ", files[i],
                  "-o ", sampsDir, "/ststats", i, ".json"))
  }
}

#Effectively the reverse of makeNSFile
##-UNTESTED Not Used for first results -##
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
taxitCreate <- function(samp, tree, json, stats, log, oDir, nSeqs, short="") {
  
  #File Manipulation and directory setup 
  refpkg <- paste0(oDir, "/", short, ".refpkg")
  dir.create(refpkg)
  
  #add components together
  file.copy(json, refpkg)
  file.copy(log, refpkg)
  file.copy(tree, refpkg)
  
  seqs <- c(read.FASTA(samp), read.FASTA(nSeqs))
  write.FASTA(seqs, paste0(refpkgN, "/fullSeqs.fasta"))
  
  conFile <- file(paste0(refpkg, "/CONTENTS.json"))
  writeLines(paste0(
    "{\n
    \t\"files\": {\n
    \t\t\"aln_fasta\": \"",  
    
    "))


}