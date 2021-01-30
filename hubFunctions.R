require(ape)

#Creates Tree-based resamples and re-runs the tree based clustering method on each
runFullMulti <- function(iFile, dateFormat, maxDs, newMark, prop=0.8, n=100, varInd=c(1,2)) {
  #@param iFile: An input full fasta file for splitting and resampling
  #@param dateFormat: Passed to multi run 
  #@param maxDs: Passed to multiGAICRun()
  #@param newMark: The marker defining new sequences from a set of time-defined seqs
  #                This must be something that exists in the headers of ONLY new sequences
  #@param prop: passed to sampleFasta() by default matches sampleFasta() default
  #@param n: passed to sampleFasta() by default matches sampleFasta() default
  
  ##-TODO: Make Generalizeable
  source("git/tn/pplacer_utils.R")
  source("~/git/tn/subT_Lib.R")
  
  #Sample Folder Set Up
  newSplit(iFile, newMark)
  oFileO <- paste0(gsub(".fasta$|.fas$", "_Old.fasta", iFile))
  nFile <- paste0(gsub(c(".fasta$|.fas$"), "_New.fasta", iFile))
  sampsDir <- paste0(gsub(".fasta$", "Samps", iFile))
  sampleFasta(oFileO, sampsDir, n=n, prop=prop, nFile =nFile)
  multiTree(sampsDir)
  multiTaxit(sampsDir)
  
  #Run Pplacer and guppy on all created refpackages
  multiPplacer(sampsDir)
  multiGuppy(sampsDir)

  #Save multiple parallel run info to output file
  oFile <- paste0(gsub(".fasta$|.fas$", "_ROB.rds", iFile))
  res <- multiGAICRun(sampsDir, maxDs=maxDs, dateFormat = dateFormat, varInd = varInd)
  saveRDS(res, oFile)
}

#Run GAIC test for tree-based data
runTreeGAIC <- function(tFile, reVars="_", varInd=c(1,2), dateFormat="%Y", addVarN=NA,
                        gFile, maxDs=NA, minB=0, nCores=1) {
  
  source("git/tn/subT_Lib.R")
  
  oT <- impTree(tFile, reVars, varInd, dateFormat, addVarN)
  oT <- growthSim(oT, gFile)
  res <- GAICRun(oT, maxDs, minB, nCores)
  return(res)
}

#Test scripts
if(F){
  
  #Set inputs for test
  reVars <- '_'
  varInd <- c(1,7,2)
  dateFormat <- "%Y"
  addVarN <- "SubType"
  
  tFile <- "~/Data/paperData/beijFullTree/beijFull_PRO_Old.fasta.treefile"
  gFile <- "~/Data/paperData/beijFullTree/beijFull_PRO_Growth.tre"
    
  oT <- impTree(tFile, reVars, varInd, dateFormat, addVarN)
  oT <- growthSim(oT, gFile)
  
  maxDs <- seq(0, 0.04, 0.001)
  res <- GAICRun(oT, maxDs, minB=0.90, nCores=8)
  
  plot(maxDs, res$GAIC, xlab="Thresholds", ylab="AIC Loss", ylim =c(-100, 3))
  lines(maxDs, res$GAIC)
  lines(maxDs, res$randAIC-res$nullAIC, col="red")
  
  #tFile <- "~/Data/Seattle/IqTree_Bootstrap/SeattleB_PRO_Filt.fasta.treefile"
  #gFile <- "~/Data/Seattle/IqTree_Bootstrap/st.tre"
  
  #tFile <- "~/Data/Tennessee/tn_ColTreeData/tn.refpkg/tn_Coltree.nwk"
  #gFile <- "~/Data/Tennessee/tn_ColTreeData/tn_ColGrowth.tre"
  
  #tFile <- "~/Data/Tennessee/tn_DiagTreeData/tn.refpkg/tn.tre"
  #gFile <- "~/Data/Tennessee/tn_DiagTreeData/tnTGrowth.tre"
  
  #tFile <- "~/Data/Seattle/stKingTrees/stKing_PRO_H_Filt.fasta.treefile"
  #gFile <- "~/Data/Seattle/stKingTrees/stKing.tre"
  
  #tFile <- "~/Data/NAlberta/IqTree_Bootstrap/na.refpkg/NorthAlbertaB_PRO_Filt.fasta.treefile"
  #gFile <- "~/Data/NAlberta/IqTree_Bootstrap/na.tre"
  
  oT <- impTree(tFile, reVars='_', varInd = c(2,3), dateFormat = "%Y-%m-%d")
  oT <- growthSim(oT, gFile)
  
  maxDs <- seq(0, 0.04, 0.001)
  res <- GAICRun(oT, maxDs, minB=0.90, nCores=8, modFormula=New~Old+Recency+DominantSubType)
  
  plot(maxDs, res$GAIC, xlab="Thresholds", ylab="AIC Loss")
  lines(maxDs, res$GAIC)
  #lines(maxDs, res$randAIC-res$nullAIC, col="red")
}
