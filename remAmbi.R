#Automated removal of ambiguous .fasta sequences

require("ape")

## USAGE: Rscript remAmbi.R [INSERT FASTA FILE HERE] ##

#Take input      
iFile <- commandArgs(trailingOnly = T)
seqs <- read.dna(iFile, format = "fasta", as.character = T)

#Check highly gapped regions for trimming
positGaps <- sapply(1:ncol(seqs), function(i) {
  iSeq <- seqs[,i]
  lGaps <- length(iSeq[iSeq%in%"-"])
  lGaps/length(iSeq)
})

#Collect Frequencies of ambiguity
fs <- sapply(1:nrow(seqs), function(i) {
  iSeq <- seqs[i,]
  lGaps <- length(iSeq[iSeq%in%"-"])
  lAmbi <- length(iSeq[iSeq%in%c("y","r","w","s","k","m","d","v","h","b","n")])
  lAmbi/(length(iSeq)-lGaps)
})

#Collect Lengths of sequences
lngs <- sapply(1:nrow(seqs), function(i) {
  iSeq <- seqs[i,]
  lGaps <- length(iSeq[iSeq%in%"-"])
  (length(iSeq) - lGaps)/length(iSeq)
})

badSeq <- union(which(fs>0.05),which(lngs<0.85))

if(length(badSeq)==0) {
  print("No unnacceptable Sequences")
} else {
  #Write only those sequences with ambiguity below 1.5% and sequence length above 85%
  write.dna(seqs[-badSeq,], colsep = "", 
            gsub(".fasta$", "_PRO.fasta", iFile), "fasta")
}

