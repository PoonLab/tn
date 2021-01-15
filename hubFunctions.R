source("git/tn/pplacer_utils.R")
source("git/tn/subT_Lib.R")
require(ape)

#Testing / Use
sampsDir <- "~/Data/paperData/tn_DiagSamps"
iFile <- "~/Data/paperData/tn_DiagTreeData/tn.refpkg/fullAln.fasta"
short <- "tn_Diag"
maxDs <- seq(0,0.04,0.001)
dateFormat <- "%Y"

#Creates Tree-based resamples and re-runs the tree based clustering method on each
runFullMulti <- function(sampsDir, iFile, short, dateFormat, maxDs, oFile) {
  #@param sampsDir: A directory to be populated with random 80% resamples of a full sequence
  #                 If none are specified, one is created
  #@param iFile: An input full fasta file for re-sampling
  #@param short: A short name for the naming of directories within the sample directory
  #@param dateFormat: Passed to multi run 
  #@param maxDs: Passed to multi run
  #@param oFile: The output RDS file for run info
  
  ##-TODO: Make Generalizeable
  source("~/git/tn/subT_Lib.R")
  
  #Sample Folder Set Up
  sampleFasta(iFile, sampsDir)
  multiTree(sampsDir)
  multiTranslator(sampsDir, "~/git/tn/translator.py")
  multiTaxit(sampsDir)
  
  #Run Pplacer and guppy on all created refpackages
  multiPplacer(paste0(sampsDir, "/", short, "_refpackages"),
               paste0(sampsDir, "/", short, "_jplaceFiles"))
  multiGuppy(paste0(sampsDir, "/", short, "_refpackages"),
             paste0(sampsDir, "/", short, "_GrowthFiles"))

  #Save multiple parallel run info to output file
  save(multiGAICRun(sampsDir, maxDs=maxDs, dateFormat = dateFormat), oFile)
}

runFullMulti(sampsDir, iFile, short, dateFormat, maxDs)


