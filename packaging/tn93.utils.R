#' A wrapper for tn93's basic run function to get an edgelist 
#' NOTE: The sequences referenced here will also be referenced in another data set (seq.info)
run.tn93 <- function(seqs.full, seq.info){
  #'@param seqs.full: The full alignment. Including sequences excluded from the tree.
  #'@param seqs.info: A data.table of sequence information with unique sequence ID's separate each sequence (see pull.headers())
  #'These ID's should match the order of and contain all sequences within the .fasta (seqs.full), to insure consistency between sequences and data.
  #'@return: An edgelist of pairwise TN93 distances calculated using tn93 binaries

  #Check inputs
  if(nrow(seq.info)!= length(seqs.full)){
    stop("The data used to match fasta sequences does not contain the same number fasta sequences")
  }
  if(!("ID"%in%colnames(seq.info))){
    stop("'ID' must be contained within the names of seq.info")
  } else {
    seq.info[, ID := as.character(seq.info$ID)]
  }

  #Overwrite labels for seqs.full to make them consistent with seq.info
  names(seqs.full) <- seq.info[, ID]
  
  #Prep temporary files
  seqs.file <- tempfile("seqs.full", fileext = ".fasta")
  edgelist.file <- tempfile("edgelist", fileext = ".csv")
  ape::write.FASTA(seqs.full, seqs.file)
  
  #Obtain TN93 distances as edgelist and clean up
  system(paste0("tn93 -t 1 -o ", edgelist.file, " ", seqs.file))
  edge.info <- data.table::fread(edgelist.file)
  
  unlink(seqs.file)
  unlink(edgelist.file)
  
  return(edge.info)      
}
